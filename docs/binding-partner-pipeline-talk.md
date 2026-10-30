# Binding-partner classification from solved structures

A talk summary of the design for a distance-only method that, given one residue and a list of solved structure files, reports how exposed the residue is and what it contacts: another protein, the same protein (another subunit), or a ligand.

---

## 1. Problem

**Inputs**

- One residue in a protein structure (chain, residue name, residue number in a reference file)
- A list of solved structure files (PDB / mmCIF)

**Outputs**

- The residue's relative solvent accessibility (a reported number; it does not gate the partner call, see 5.3)
- What the residue contacts, if anything:
  - another protein (different subunit / different chain, different sequence entity)
  - the same protein (homomer: another chain of the same entity)
  - a ligand (non-polymer hetero group)
  - optionally: ion, nucleic acid, same-chain packing
- If it contacts nothing: `none` (exposed) or `buried` (not exposed)

Structures entering this pipeline are already solved. The method does not dock, score energies, or invent partners. It only reads coordinates that are already present.

---

## 2. Relation to the existing electrostatics pipeline

The current repo (`electrostatics-residue-pipeline`) already assumes a known split:

- `structure.protein_chains`
- `structure.partner_chains`

Binding electrostatics is scored as:

```text
ΔG_binding ≈ E(complex) − E(protein) − E(partner)
```

That split is currently written by hand in config (demo: barnase–barstar, PDB `1BRS`, chains A and D). The method discussed here is an upstream classifier: discover the partner type (and chain / hetero identity) from structure files, so that a later step can write those fields back into config.

The classifier is built standalone. It only emits labels and evidence; it does not touch the electrostatics code and computes no energies.

---

## 3. Clarification: distance only

The partner call uses **distance constraints on deposited coordinates only**.

1. **Accessibility** — How exposed is the residue on its own chain? (Reported, not a gate.)
2. **Partner identity** — What atoms are within a short shell, and what do those atoms belong to?

No docking, no binding free-energy scoring, and no dynamics enter this classification step.

---

## 4. Labels

Chain ID is only a handle. The class comes from what the contacting atoms belong to. A site may support more than one label in the same structure (e.g. protein partner and ligand).

| Label | Meaning |
| --- | --- |
| `other_protein` | Standard amino acids on another polymer chain with a different sequence entity |
| `same_protein` | Standard amino acids on another chain of the same entity (homomer) |
| `same_chain` | Contacts on the query chain far enough in sequence to be outside local packing |
| `ligand` | Non-polymer hetero group after solvent, crystallization-agent and glycan filtering |
| `ion` | Metal or monoatomic ion (kept separate; AIMNet2 already refuses unsupported metals) |
| `nucleic_acid` | DNA or RNA polymer chain (its own label, not forced into protein or ligand) |
| `none` | Nothing in the contact shell, and the residue is exposed (relative SASA at or above the cutoff) |
| `buried` | Nothing in the contact shell, and the residue is below the SASA cutoff |

`same_chain` is intramolecular packing or domain–domain contact. It should **not** be written into `partner_chains`.

Glycans (mmCIF branched entities, and components such as `NAG`, `NDG`, `BMA`, `MAN`, `FUC`) are ignored entirely. They never count as a partner and never produce a `ligand` label.

---

## 5. Method outline

### 5.1 Place the residue across files

Fix the residue in the reference structure (chain, name, number). Read that chain's sequence. In every other file, place the site by local sequence alignment to each protein chain. Author residue numbers can drift between depositions; the alignment index is the identity carried forward.

Every chain that aligns gives one placed copy of the site (an asymmetric unit with six HA1 chains gives six). Drop a file when the position is absent from every chain, when the aligned residue type differs from the reference residue, or when a chain's identity to the reference chain is below 30%.

### 5.2 Prefer the biological unit

Prefer the first biological assembly (mmCIF assembly operators, or PDB `REMARK 350` / `BIOMT`). Contacts that appear only through a crystallographic symmetry mate outside that assembly are recorded as crystal-only for that file and do not count as binding partners. Files that are already a biological unit (for example the `_3mer` extracts, which carry no assembly metadata) are flagged by a filename glob in config; symmetry is not applied to them again. Use the first model of an NMR ensemble and the highest-occupancy altloc.

### 5.3 Accessibility (SASA)

Compute relative solvent-accessible surface area for the query residue on the query chain alone (protein atoms of that chain only: no partners, ligands, or waters), for each placed copy. The value is reported.

It never blocks a partner call. A residue with a qualifying contact in 5.4 is at an interface whatever its accessibility. On `1BRS`, A:GLU73 is only about 10% accessible on chain A alone (it sits in the barnase active-site pocket) and is still the textbook protein-contact residue, so a hard gate would reject the demo.

SASA decides only the no-contact case: relative SASA at or above the cutoff (a typical starting point is ~20% of the residue's maximum area) gives `none`; below it gives `buried`.

### 5.4 Contact shell (distance)

From the heavy atoms of the query residue, collect every other heavy atom within **4.5 Å**. That shell is an interface contact distance. (The existing 8 Å cluster cutoff in `cluster.py` is a second shell for quantum neighborhood construction and is too wide for partner calling.)

Require a minimum number of atom pairs per partner so a single grazing approach does not qualify: **3** for protein, nucleic-acid, ligand and same-chain partners, **1** for ions (a coordinated ion touches only one or two atoms). A partner is one chain for polymers and one hetero residue for ligands and ions. Ignore hydrogens.

### 5.5 Exclude non-partners

Remove from the shell:

- the query residue itself
- sequence neighbors on the same chain with \(|i - j| \le 5\)
- waters
- glycans (see section 4)

Residues inside a polymer chain, including modified ones (`MSE`, `SEP`, `TPO`, `PTR`, etc.), stay amino acids. Polymer membership comes from the entity type, not from the `HETATM` record.

### 5.6 Classify what remains

- Amino-acid polymer, different chain: in mmCIF the entity decides. Same entity → `same_protein`; different entity → `other_protein`. PDB files carry no entity table, so entities are assigned by identical sequence, and a different entity with sequence identity **≥ 90%** to the query chain also counts as `same_protein`.
- Same chain and \(|i - j| > 5\) → `same_chain`
- Nucleic-acid polymer → `nucleic_acid`
- Non-polymer hetero group, after removing glycans and a fixed crystallization list (`HOH`, `SO4`, `PO4`, `GOL`, `EDO`, `PEG`, `ACT`, `DMS`, `MPD`, `TRS`, `MES`, `CIT`, and usual buffers):
  - metal or monoatomic ion (`ZN`, `MG`, `CA`, `MN`, `FE`, `NA`, `CL`, …) → `ion`
  - anything else → `ligand`

Record partner id, minimum distance, atom-pair count, and whether the contact sits in the biological assembly.

### 5.7 Combine placed copies and files

One crystal is weak evidence, and one PDB entry can arrive as several files, so evidence is counted per PDB entry, not per file.

- Within a file, a label is present when it appears in at least half of the placed copies.
- Files from the same PDB entry (for example an asymmetric-unit file and its `_3mer` extract) form one entry. A label is present for the entry if it is present in any of its files.
- A label is **supported** when it is present in at least **two entries** and in at least **30%** of the entries that could contain that partner (both conditions must hold).

An entry could contain a partner when its assembly has:

- `other_protein`: another protein chain of a different entity
- `same_protein`: another chain of the query chain's entity
- `ligand` or `ion`: that same chemical component
- `nucleic_acid`: a nucleic-acid chain
- `same_chain`, `none`, `buried`: every entry counts

Ligand and ion identity is the chemical component, so `SO4` in one crystal never supports `FMN` in another.

If a single entry covers the site, report what it shows with confidence `low`. Otherwise confidence is `ok` for supported labels. Emit one evidence row per file and placed copy, and one consensus row per label (a `none` or `buried` row when no label is supported). The consensus accessibility is the median relative SASA across entries.

---

## 6. How labels would feed the electrostatics run

Out of scope for the classifier build, which only writes labels. This is the intended use of them in a later step:

| Consensus label | Config / scoring consequence |
| --- | --- |
| `other_protein` or `same_protein` | Set `protein_chains` to the query chain and `partner_chains` to the contacting chain(s). Keep complex − protein − partner APBS. Homomers use the same sequence on both sides. |
| `ligand` | Keep that hetero group through preparation. Today `prepare.py` calls `removeHeterogens`, so ligand partners never reach APBS; preparation must retain the ligand and score protein+ligand, protein alone, and ligand alone. |
| `ion` | Record it; leave AIMNet2 refused when unsupported elements are present (existing metal refusal). |
| `nucleic_acid` | Record it; no protein–protein split applies. |
| `same_chain`, `buried` or `none` | Do **not** report a binding \(\Delta\Delta G\). The site is intramolecular, buried, or solvent-only; the binding split is the wrong quantity. |

The thresholds (4.5 Å, minimum atom pairs, sequence gap of 5, 90% identity, 30% reference identity, two-entry and 30% vote minimums, ~20% relative SASA, and the crystallization, glycan and ion component lists) live in a `partner:` block of the config.

---

## 7. Sanity check on the demo

The hand-written config the classifier should reproduce for `1BRS` is:

```yaml
structure:
  pdb: data/demo/1BRS.pdb
  protein_chains: ["A"]
  partner_chains: ["D"]
```

The first draft of this plan expected residue **A:GLU73** to classify as `other_protein` against chain **D** (barstar). Measured on the file, at the settings above it does not. Within 4.5 Å of A:GLU73 there is one non-water atom pair to barstar (`OE1` to D:ASP39 `OD1`, 4.49 Å) plus three contacts to a water (`HOH91`), which is below the three-pair minimum. It only reaches barstar with a shell of about 5 Å (3 pairs). Its relative SASA on chain A alone is about 0.10.

Other residues at the same interface classify as expected, `other_protein` against chain D: A:ARG83, A:ARG87, A:LYS27 and A:HIS102. A:ARG87 is only about 2% accessible on chain A alone and is still called, which is why SASA does not gate the partner call.

Open: keep GLU73 as the demo with a wider shell, or switch the demo residue to one of the four above. With the single `1BRS` file the consensus confidence is `low`, because one entry cannot meet the two-entry minimum.

---

## 8. Agreed scope (short form)

On already-solved structures:

1. Report how exposed the residue is (relative SASA); this never blocks the partner call.
2. Use a short heavy-atom distance shell to say whether the partner is another protein, the same protein (another subunit), or a ligand (or ion / nucleic acid), and `none` or `buried` when nothing is in contact.
3. Combine evidence per PDB entry across the file list and report a consensus with a confidence.

That classifier is upstream of the existing electrostatics residue-substitution pipeline and exists to populate partner identity rather than to compute binding energies.

---

## 9. Changes after the first review

- SASA reports a number instead of gating partner calls; `buried` and `none` are the no-contact outcomes.
- Glycans are ignored; ligand status comes from the entity type, so `MSE` inside a chain stays an amino acid.
- Minimum atom pairs are per class (3, but 1 for ions).
- Entity ID decides `same_protein` versus `other_protein` in mmCIF; the 90% identity rule is the PDB fallback.
- The vote uses "and" over PDB entries, with a per-partner denominator, and gives one consensus row per label.
- Placement requires the same residue type and at least 30% identity to the reference chain.
- The `1BRS` demo residue GLU73 does not meet the 4.5 Å / 3-pair rule (see section 7).
- PDB entries are identified by `_entry.id`, because extracts such as `H1_1RU7_3mer.cif` carry their file name as the data-block name.
