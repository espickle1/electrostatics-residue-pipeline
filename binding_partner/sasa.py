"""Relative solvent accessibility of one residue on its own chain."""

from __future__ import annotations

import math

import freesasa

from .config import Config
from .structure import ChainInfo

freesasa.setVerbosity(freesasa.nowarnings)


def relative_sasa(chain: ChainInfo, index: int, cfg: Config) -> float | None:
    """Relative SASA of chain.residues[index] with only that chain's polymer atoms present.

    None when freesasa has no reference area for the residue type (e.g. MSE).
    """
    structure = freesasa.Structure()
    target = ""
    for i, res in enumerate(chain.residues):
        number = f"{res.seqid.num}{res.seqid.icode.strip()}"
        if i == index:
            target = number
        for atom in res:
            structure.addAtom(atom.name, res.name, number, "A", atom.pos.x, atom.pos.y, atom.pos.z)
    params = freesasa.Parameters({"probe-radius": cfg.probe_radius})
    area = freesasa.calc(structure, params).residueAreas()["A"][target]
    return None if math.isnan(area.relativeTotal) else area.relativeTotal
