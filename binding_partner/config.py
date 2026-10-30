"""Thresholds and component lists for the binding-partner classifier."""

from __future__ import annotations

from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any

import yaml


def _names(text: str) -> frozenset[str]:
    return frozenset(text.split())


@dataclass(frozen=True)
class Config:
    # Accessibility (reported; only decides none vs buried when nothing is in contact)
    surface_rsa_cutoff: float = 0.20
    probe_radius: float = 1.4
    # Contact shell
    contact_distance: float = 4.5
    min_pairs: int = 3
    min_pairs_ion: int = 1
    same_chain_gap: int = 5
    same_protein_identity: float = 0.90
    # Placement across files
    min_ref_identity: float = 0.30
    require_same_residue: bool = True
    biounit_globs: tuple[str, ...] = ()
    # Combining copies, files and entries
    copy_fraction: float = 0.5
    min_entries: int = 2
    min_entry_fraction: float = 0.30
    # Component lists (residue names)
    excluded_components: frozenset[str] = _names(
        "HOH DOD SO4 PO4 GOL EDO PEG PGE PG4 1PE P6G ACT ACY FMT DMS MPD TRS MES CIT "
        "EPE IMD NO3 BME DTT TAR MLI"
    )
    glycan_components: frozenset[str] = _names(
        "NAG NDG NGA A2G BMA MAN FUC FUL GAL GLA BGC GLC SIA XYP BDP"
    )
    ion_components: frozenset[str] = _names(
        "ZN MG CA MN FE FE2 CU CU1 NI CO CD HG K NA CL BR IOD CS RB LI SR BA"
    )

    def min_pairs_for(self, label: str) -> int:
        return self.min_pairs_ion if label == "ion" else self.min_pairs


def load_config(path: str | Path | None = None) -> Config:
    """Read the ``partner:`` block of a YAML file (or the whole file) over the defaults."""
    if path is None:
        return Config()
    with open(path, encoding="utf-8") as fh:
        raw: Any = yaml.safe_load(fh) or {}
    block = raw.get("partner", raw) if isinstance(raw, dict) else raw
    if not isinstance(block, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    unknown = sorted(set(block) - {f.name for f in fields(Config)})
    if unknown:
        raise ValueError(f"Unknown partner config keys in {path}: {unknown}")
    kwargs: dict[str, Any] = {}
    for key, value in block.items():
        if key.endswith("_components"):
            kwargs[key] = frozenset(str(v).upper() for v in value)
        elif key == "biounit_globs":
            kwargs[key] = tuple(value)
        else:
            kwargs[key] = value
    return Config(**kwargs)
