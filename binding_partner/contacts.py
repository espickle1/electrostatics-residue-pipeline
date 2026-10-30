"""Distance shell around the query residue and classification of what is in it."""

from __future__ import annotations

from dataclasses import dataclass

import gemmi

from .align import pair_identity
from .config import Config
from .structure import NUCLEIC, PROTEIN, WATER_NAMES, ChainInfo, seq_pos


@dataclass(frozen=True)
class Partner:
    label: str
    partner_id: str  # chain name for polymers; chain:resnameNUM for hetero groups
    component: str  # chemical component for ligand / ion, empty otherwise
    n_pairs: int
    min_dist: float
    entity_id: str


def protein_label(query: ChainInfo, other: ChainInfo, mmcif: bool, cfg: Config) -> str:
    """mmCIF: the entity decides. PDB has no entity table, so sequence identity does."""
    if mmcif and query.entity_id and other.entity_id:
        return "same_protein" if query.entity_id == other.entity_id else "other_protein"
    identity = pair_identity(query.seq, other.seq)
    return "same_protein" if identity >= cfg.same_protein_identity else "other_protein"


def hetero_label(res: gemmi.Residue, cfg: Config) -> str | None:
    """ligand / ion for a non-polymer residue, None when it never counts as a partner."""
    if res.entity_type in (
        gemmi.EntityType.Polymer,
        gemmi.EntityType.Water,
        gemmi.EntityType.Branched,
    ):
        return None
    name = res.name
    if name in WATER_NAMES or name in cfg.glycan_components or name in cfg.excluded_components:
        return None
    if name in cfg.ion_components or len(res) == 1:
        return "ion"
    return "ligand"


def find_partners(
    model: gemmi.Model,
    ns: gemmi.NeighborSearch,
    chains: dict[str, ChainInfo],
    query_chain: str,
    query_index: int,
    cfg: Config,
    mmcif: bool,
    symmetry: bool = False,
) -> list[Partner]:
    """Partners of chains[query_chain].residues[query_index] within the contact shell.

    symmetry=False keeps contacts inside the model. symmetry=True is for a neighbor
    search built with the real unit cell and keeps only contacts through symmetry mates.
    """
    query = chains[query_chain]
    qres = query.residues[query_index]
    pairs: dict[tuple[str, str, str], list[float]] = {}
    entities: dict[tuple[str, str, str], str] = {}
    for atom in qres:
        for mark in ns.find_atoms(atom.pos, "\0", radius=cfg.contact_distance):
            if (mark.image_idx != 0) != symmetry:
                continue
            dist = atom.pos.dist(mark.pos)
            if dist > cfg.contact_distance:
                continue
            cra = mark.to_cra(model)
            res, chain_name = cra.residue, cra.chain.name
            same_copy = chain_name == query_chain and mark.image_idx == 0
            if same_copy and res.seqid == qres.seqid and res.name == qres.name:
                continue
            if res.entity_type == gemmi.EntityType.Polymer:
                info = chains.get(chain_name)
                if info is None:
                    continue
                if info.kind == NUCLEIC:
                    label = "nucleic_acid"
                elif info.kind != PROTEIN:
                    continue
                elif same_copy:
                    if abs(seq_pos(res) - seq_pos(qres)) <= cfg.same_chain_gap:
                        continue
                    label = "same_chain"
                else:
                    label = protein_label(query, info, mmcif, cfg)
                partner_id, component, entity_id = chain_name, "", info.entity_id
            else:
                label = hetero_label(res, cfg)
                if label is None:
                    continue
                partner_id = f"{chain_name}:{res.name}{res.seqid.num}{res.seqid.icode.strip()}"
                component, entity_id = res.name, ""
            if mark.image_idx:
                partner_id += f"@sym{mark.image_idx}"
            key = (label, partner_id, component)
            pairs.setdefault(key, []).append(dist)
            entities[key] = entity_id
    partners = [
        Partner(label, pid, comp, len(d), min(d), entities[(label, pid, comp)])
        for (label, pid, comp), d in pairs.items()
        if len(d) >= cfg.min_pairs_for(label)
    ]
    partners.sort(key=lambda p: (p.label, -p.n_pairs, p.partner_id))
    return partners
