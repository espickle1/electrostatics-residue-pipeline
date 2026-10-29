# Binding-partner classification from solved structures

A talk summary of the design for a distance-only method that, given one residue and a list of solved structure files, decides whether that residue is on the surface and, if so, what it contacts: another protein, the same protein (another subunit), or a ligand.

---

## 1. Problem

**Inputs**

- One residue in a protein structure (chain, residue name, residue number in a reference file)
- A list of solved structure files (PDB / mmCIF)

**Outputs**

- Whether the residue is on the surface
- If it sits at a residue–ligand or residue–protein contact surface, what the partner is:
  - another protein (different subunit / different chain, different sequence entity)
  - the same protein (homomer: another chain of the same entity)
  - a ligand (non-polymer hetero group)
  - optionally: ion, same-chain packing, or none (solvent only)

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

That split is currently written by hand in config (demo: barnase–barstar, PDB `1BRS`, chains A and D). The method discussed here is an upstream classifier: discover the partner type (and chain / hetero identity) from structure files, then write those fields back into config when a protein–protein or protein–ligand split is warranted.

---

## 3. Clarification: distance only

The partner call uses **distance constraints on deposited coordinates only**.

1. **Surface** — Is the residue solvent-exposed?
2. **Partner identity** — If it is at a contact surface, what atoms are within a short shell, and what do those atoms belong to?

No docking, no binding free-energy scoring, and no dynamics enter this classification step. Energies remain the job of the downstream APBS / AIMNet2 pipeline once chains and ligands are known.

---

## 4. Labels

Chain ID is only a handle. The class comes from what the contacting atoms belong to. A site may support more than one label in the same structure (e.g. protein partner and ligand).

| Label | Meaning |
| --- | --- |
| `other_protein` | Standard amino acids on another polymer chain with a different sequence entity |
| `same_protein` | Standard amino acids on another chain of the same entity (homomer) |
| `same_chain` | Contacts on the query chain far enough in sequence to be outside local packing |
| `ligand` | Non-polymer hetero group after solvent / crystallization-agent filtering |
| `ion` | Metal or monoatomic ion (kept separate; AIMNet2 already refuses unsupported metals) |
| `none` | Surface residue with nothing left in the contact shell (solvent-facing) |

`same_chain` is intramolecular packing or domain–domain contact. It should **not** be written into `partner_chains`. Nucleic acid, if present, should get its own label rather than being forced into protein or ligand.

---

## 5. Method outline

### 5.1 Place the residue across files

Fix the residue in the reference structure (chain, name, number). Read that chain’s sequence. In every other file, place the site by sequence alignment to each polymer chain. Author residue numbers can drift between depositions; the alignment index is the identity carried forward. Drop a file only when that position is absent.

### 5.2 Prefer the biological unit

Prefer the first biological assembly (mmCIF assembly operators, or PDB `REMARK 350` / `BIOMT`). Contacts that appear only through a crystallographic symmetry mate outside that assembly are tagged as crystal contacts and do not count as binding partners. If a file is already a biological unit, do not apply symmetry again. Use the first model of an NMR ensemble and the highest-occupancy altloc.

### 5.3 Surface check (SASA)

Compute solvent-accessible surface area for the query residue on the query chain alone. The residue is on the surface when its relative accessibility is above a cutoff (a typical starting point is ~20% of the residue’s maximum area). Residues below the cutoff are buried and are not treated as binding-surface sites for this classifier.

### 5.4 Contact shell (distance)

From the heavy atoms of the query residue, collect every other heavy atom within **4.5 Å**. That shell is an interface contact distance. (The existing 8 Å cluster cutoff in `cluster.py` is a second shell for quantum neighborhood construction and is too wide for partner calling.)

Require at least **3** atom pairs so a single grazing approach does not qualify. Ignore hydrogens.

### 5.5 Exclude non-partners

Remove from the shell:

- the query residue itself
- sequence neighbors on the same chain with \(|i - j| \le 5\)
- waters
- modified forms of the residue itself when they are the query site (`MSE`, `SEP`, `TPO`, `PTR`, etc.)

### 5.6 Classify what remains

- Amino-acid polymer, different chain, sequence identity with the query chain **&lt; 90%** (or a different mmCIF entity id) → `other_protein`
- Amino-acid polymer, different chain, same entity or identity **≥ 90%** → `same_protein`
- Same chain and \(|i - j| > 5\) → `same_chain`
- Non-polymer hetero group → `ligand`, after removing a fixed crystallization list (`HOH`, `SO4`, `PO4`, `GOL`, `EDO`, `PEG`, `ACT`, `DMS`, `MPD`, `TRS`, `MES`, `CIT`, and usual buffers)
- Metal or monoatomic ion (`ZN`, `MG`, `CA`, `MN`, `FE`, `NA`, `CL`, …) → `ion`

Record partner id, minimum distance, atom-pair count, and whether the contact sits in the biological assembly.

### 5.7 Vote across the file list

One crystal is weak evidence. Prefer a class that appears in the biological assembly in at least two files, or in a substantial fraction of files where the residue is resolved (e.g. ≥ 30%). Ligand identity should recur as the same chemical component, not a different buffer in each crystal. Emit one evidence row per file and one consensus row per residue. Keep confidence low when only one structure covers the site.

---

## 6. How labels feed the electrostatics run

| Consensus label | Config / scoring consequence |
| --- | --- |
| `other_protein` or `same_protein` | Set `protein_chains` to the query chain and `partner_chains` to the contacting chain(s). Keep complex − protein − partner APBS. Homomers use the same sequence on both sides. |
| `ligand` | Keep that hetero group through preparation. Today `prepare.py` calls `removeHeterogens`, so ligand partners never reach APBS; preparation must retain the ligand and score protein+ligand, protein alone, and ligand alone. |
| `ion` | Record it; leave AIMNet2 refused when unsupported elements are present (existing metal refusal). |
| `same_chain` or `none` | Do **not** report a binding \(\Delta\Delta G\). The site is intramolecular or solvent-only; the binding split is the wrong quantity. |

Suggested thresholds (4.5 Å, ≥ 3 atom pairs, sequence gap of 5, 90% identity, two-structure minimum, ~20% relative SASA) should live in config next to the current chain lists.

---

## 7. Sanity check on the demo

On PDB `1BRS`, residue **A:GLU73** should classify as `other_protein` against chain **D** (barstar), with no ligand partner, matching the hand-written config:

```yaml
structure:
  pdb: data/demo/1BRS.pdb
  protein_chains: ["A"]
  partner_chains: ["D"]
```

---

## 8. Agreed scope (short form)

On already-solved structures:

1. Decide if the residue is on the surface (SASA).
2. If it is at a contact surface, use a short heavy-atom distance shell to say whether the partner is another protein, the same protein (another subunit), or a ligand (or ion / none).

That classifier is upstream of the existing electrostatics residue-substitution pipeline and exists to populate partner identity rather than to compute binding energies.
