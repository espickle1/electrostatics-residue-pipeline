"""Classify the reference residue in each structure file."""

from __future__ import annotations

from dataclasses import dataclass, field
from fnmatch import fnmatch
from pathlib import Path

import gemmi

from .align import Placement, Placer
from .config import Config
from .contacts import Partner, find_partners, hetero_label, protein_label
from .sasa import relative_sasa
from .structure import (
    NUCLEIC,
    PROTEIN,
    ChainInfo,
    build_model,
    chain_infos,
    entry_id,
    load_structure,
)

RESOLVED = "resolved"
Key = tuple[str, str]  # (label, component); component is empty except for ligand / ion


@dataclass(frozen=True)
class Reference:
    path: str
    chain: str
    resname: str
    resnum: int
    icode: str
    seq: str
    index: int


@dataclass
class Site:
    """One placed copy of the reference residue (one aligned chain)."""

    chain: str
    resname: str
    resnum: int
    identity: float
    rsa: float | None
    partners: list[Partner]
    baseline: str  # none / buried when no partner qualifies, else empty

    def keys(self) -> set[Key]:
        return {(p.label, p.component) for p in self.partners}


@dataclass
class FileResult:
    path: str
    entry_id: str = ""
    status: str = ""
    detail: str = ""
    model_kind: str = ""  # assembly or as_is
    sites: list[Site] = field(default_factory=list)
    available: set[Key] = field(default_factory=set)  # partners that could have been found
    crystal_only: list[Partner] = field(default_factory=list)


def make_reference(
    path: str | Path, chain: str, resname: str, resnum: int, icode: str = ""
) -> Reference:
    st = load_structure(path)
    chains = chain_infos(st[0])
    info = chains.get(chain)
    if info is None or info.kind != PROTEIN:
        raise ValueError(f"{path}: no protein chain {chain!r} (polymer chains: {sorted(chains)})")
    for index, res in enumerate(info.residues):
        if res.seqid.num == resnum and res.seqid.icode.strip() == icode.strip():
            if res.name != resname.upper():
                raise ValueError(
                    f"{path}: chain {chain} residue {resnum}{icode} is {res.name}, "
                    f"not {resname.upper()}"
                )
            return Reference(str(path), chain, res.name, resnum, icode.strip(), info.seq, index)
    raise ValueError(f"{path}: chain {chain} has no residue {resnum}{icode}")


def _place(
    chains: dict[str, ChainInfo], placer: Placer
) -> tuple[list[tuple[str, Placement]], list[str]]:
    placed: list[tuple[str, Placement]] = []
    reasons: list[str] = []
    for name, info in chains.items():
        if info.kind != PROTEIN:
            continue
        out = placer.place(info.seq)
        if isinstance(out, Placement):
            placed.append((name, out))
        else:
            reasons.append(out)
    return placed, reasons


def _available(
    model: gemmi.Model,
    chains: dict[str, ChainInfo],
    placed: list[tuple[str, Placement]],
    mmcif: bool,
    cfg: Config,
) -> set[Key]:
    """Partners the model could show: the denominator for the vote."""
    keys: set[Key] = {("same_chain", "")}
    for name, _ in placed:
        for other_name, other in chains.items():
            if other_name == name:
                continue
            if other.kind == PROTEIN:
                keys.add((protein_label(chains[name], other, mmcif, cfg), ""))
            elif other.kind == NUCLEIC:
                keys.add(("nucleic_acid", ""))
    for chain in model:
        for res in chain:
            label = hetero_label(res, cfg)
            if label:
                keys.add((label, res.name))
    return keys


def _crystal_only(
    st: gemmi.Structure,
    placer: Placer,
    sites: list[Site],
    mmcif: bool,
    cfg: Config,
) -> list[Partner]:
    """Partners reached only through symmetry mates of the asymmetric unit."""
    asu = st[0]
    chains = chain_infos(asu)
    placed, _ = _place(chains, placer)
    if not placed:
        return []
    ns = gemmi.NeighborSearch(asu, st.cell, 5).populate()
    in_assembly = {key for site in sites for key in site.keys()}
    found: dict[Key, Partner] = {}
    for name, placement in placed:
        for p in find_partners(asu, ns, chains, name, placement.index, cfg, mmcif, symmetry=True):
            key = (p.label, p.component)
            if key in in_assembly:
                continue
            if key not in found or p.n_pairs > found[key].n_pairs:
                found[key] = p
    return list(found.values())


def analyze_file(path: str | Path, ref: Reference, cfg: Config) -> FileResult:
    path = Path(path)
    st = load_structure(path)
    result = FileResult(path=str(path), entry_id=entry_id(st, path))
    mmcif = st.input_format == gemmi.CoorFormat.Mmcif
    already_biounit = any(fnmatch(path.name, g) for g in cfg.biounit_globs)
    model, result.model_kind = build_model(st, already_biounit)
    chains = chain_infos(model)
    placer = Placer(ref.seq, ref.index, cfg)
    placed, reasons = _place(chains, placer)
    if not placed:
        for reason in ("residue_mismatch", "low_identity", "position_absent"):
            if reason in reasons:
                result.status = reason
                break
        else:
            result.status = "no_protein_chain"
        result.detail = f"{len(reasons)} protein chain(s) tried"
        return result

    ns = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate()
    for name, placement in placed:
        info = chains[name]
        qres = info.residues[placement.index]
        partners = find_partners(model, ns, chains, name, placement.index, cfg, mmcif)
        rsa = relative_sasa(info, placement.index, cfg)
        if partners:
            baseline = ""
        else:
            baseline = "none" if rsa is not None and rsa >= cfg.surface_rsa_cutoff else "buried"
        result.sites.append(
            Site(name, qres.name, qres.seqid.num, placement.identity, rsa, partners, baseline)
        )
    result.available = _available(model, chains, placed, mmcif, cfg)
    if not already_biounit and st.cell.is_crystal():
        result.crystal_only = _crystal_only(st, placer, result.sites, mmcif, cfg)
    result.status = RESOLVED
    return result


def analyze_files(paths: list[str | Path], ref: Reference, cfg: Config) -> list[FileResult]:
    results: list[FileResult] = []
    for path in paths:
        try:
            results.append(analyze_file(path, ref, cfg))
        except Exception as exc:  # reported per file in files.csv, never silent
            results.append(
                FileResult(path=str(path), status="error", detail=f"{type(exc).__name__}: {exc}")
            )
    return results
