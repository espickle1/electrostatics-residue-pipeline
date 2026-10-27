"""Path layout on the Modal volume and local results directory."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

# Mount point inside Modal containers
VOLUME_ROOT = Path("/results")
RUN_NAME_DEFAULT = "vertical_slice"


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    return cfg


def run_root(cfg: dict[str, Any], volume_root: Path = VOLUME_ROOT) -> Path:
    name = cfg.get("run_name", RUN_NAME_DEFAULT)
    return volume_root / name


def ensure_tree(root: Path) -> None:
    for sub in (
        "structures",
        "frames",
        "prot",
        "apbs",
        "clusters",
        "aimnet",
        "results",
        "manifests",
    ):
        (root / sub).mkdir(parents=True, exist_ok=True)


def structure_pdb(root: Path, variant: str) -> Path:
    return root / "structures" / f"{variant}.pdb"


def frame_pdb(root: Path, variant: str, frame: str = "0") -> Path:
    return root / "frames" / variant / f"frame_{frame}.pdb"


def prot_dir(root: Path, variant: str, frame: str = "0") -> Path:
    d = root / "prot" / variant / f"frame_{frame}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def apbs_dir(root: Path, variant: str, frame: str = "0") -> Path:
    d = root / "apbs" / variant / f"frame_{frame}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cluster_dir(root: Path, variant: str, frame: str = "0") -> Path:
    d = root / "clusters" / variant / f"frame_{frame}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def aimnet_dir(root: Path, variant: str, frame: str = "0") -> Path:
    d = root / "aimnet" / variant / f"frame_{frame}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def manifest_path(root: Path, variant: str, frame: str = "0") -> Path:
    return root / "manifests" / f"{variant}_frame_{frame}.yaml"


def results_dir(root: Path) -> Path:
    d = root / "results"
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_results_dir() -> Path:
    d = Path(__file__).resolve().parents[1] / "results"
    d.mkdir(parents=True, exist_ok=True)
    return d


def variant_ids(cfg: dict[str, Any]) -> list[str]:
    ids = ["WT"]
    for v in cfg.get("variants", []):
        ids.append(v["id"])
    return ids
