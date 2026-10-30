"""Read a structure file and reduce it to one model with one conformer per residue."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi

PROTEIN, NUCLEIC, OTHER = "protein", "nucleic", "other"
WATER_NAMES = frozenset({"HOH", "DOD"})
STANDARD_CODES = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class ChainInfo:
    name: str
    kind: str
    residues: list[gemmi.Residue]  # polymer residues only
    entity_id: str  # empty for PDB files, which carry no entity table
    seq: str  # one-letter sequence of the polymer residues (protein chains only)


def load_structure(path: str | Path) -> gemmi.Structure:
    """First model, hydrogens removed, highest-occupancy altloc per residue."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_hydrogens()
    if len(st) > 1:
        del st[1:]
    _keep_best_altloc(st[0])
    return st


def _keep_best_altloc(model: gemmi.Model) -> None:
    for chain in model:
        for res in chain:
            occupancy: dict[str, float] = {}
            for atom in res:
                if atom.altloc != "\0":
                    occupancy[atom.altloc] = occupancy.get(atom.altloc, 0.0) + atom.occ
            if len(occupancy) < 2:
                continue
            # Sorted first so an occupancy tie keeps the earliest altloc.
            best = max(sorted(occupancy), key=lambda alt: occupancy[alt])
            drop = [(a.name, a.altloc) for a in res if a.altloc not in ("\0", best)]
            for name, altloc in drop:
                res.remove_atom(name, altloc)


def entry_id(st: gemmi.Structure, path: str | Path) -> str:
    """PDB id: _entry.id when present (kept by extracts such as H1_1RU7_3mer.cif), else the name."""
    info = dict(st.info.items())  # gemmi.InfoMap has no .get
    return (info.get("_entry.id") or st.name or Path(path).stem).strip().upper()


def build_model(st: gemmi.Structure, already_biounit: bool) -> tuple[gemmi.Model, str]:
    """The first biological assembly when the file defines one, else the model as read."""
    if not already_biounit and len(st.assemblies) > 0:
        model = gemmi.make_assembly(
            st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.Short
        )
        return model, "assembly"
    return st[0], "as_is"


def one_letter(res: gemmi.Residue) -> str:
    info = gemmi.find_tabulated_residue(res.name)
    code = info.one_letter_code.strip().upper() if info and info.is_amino_acid() else ""
    return code if code in STANDARD_CODES else "X"


def seq_pos(res: gemmi.Residue) -> int:
    """Sequence position: mmCIF label_seq_id when present, else the residue number."""
    return res.label_seq if res.label_seq is not None else res.seqid.num


def chain_infos(model: gemmi.Model) -> dict[str, ChainInfo]:
    """Polymer chains of a model. Polymer membership comes from the entity type, not HETATM."""
    out: dict[str, ChainInfo] = {}
    for chain in model:
        residues = [r for r in chain if r.entity_type == gemmi.EntityType.Polymer]
        if not residues:
            continue
        infos = [gemmi.find_tabulated_residue(r.name) for r in residues]
        n_aa = sum(1 for i in infos if i and i.is_amino_acid())
        n_na = sum(1 for i in infos if i and i.is_nucleic_acid())
        if n_aa and n_aa >= n_na:
            kind = PROTEIN
        elif n_na:
            kind = NUCLEIC
        else:
            kind = OTHER
        seq = "".join(one_letter(r) for r in residues) if kind == PROTEIN else ""
        out[chain.name] = ChainInfo(chain.name, kind, residues, residues[0].entity_id, seq)
    return out
