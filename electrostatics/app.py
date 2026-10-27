"""Modal app: classical CPU + AIMNet2 GPU vertical slice."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Optional

import modal
import yaml

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
SRC = HERE / "src"
MINUTES = 60

# Pinned stack notes (also written next to results on each run).
CLASSICAL_PIN = {
    "python": "3.11",
    "pdb2pqr": "conda-forge",
    "apbs": "conda-forge",
    "openmm": "conda-forge",
    "pdbfixer": "conda-forge",
    "mdtraj": "conda-forge",
    "griddataformats": "conda-forge",
}
AIMNET_PIN = {
    "python": "3.11",
    "aimnet": "pip aimnet[ase]",
    "torch": "cu126",
}

app = modal.App(name="electrostatics-residue-pipeline")

results_volume = modal.Volume.from_name("electrostatics-results", create_if_missing=True)
RESULTS = Path("/results")

classical_image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        "pdb2pqr",
        "apbs",
        "openmm",
        "pdbfixer",
        "mdtraj",
        "numpy",
        "pandas",
        "griddataformats",
        "pyyaml",
        "pip",
        channels=["conda-forge"],
    )
    .env({"PYTHONPATH": "/pkg"})
    .add_local_dir(str(SRC), remote_path="/pkg")
)

aimnet_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch",
        index_url="https://download.pytorch.org/whl/cu126",
    )
    .pip_install("aimnet[ase]", "pyyaml", "numpy", "pandas")
    .env({"PYTHONPATH": "/pkg"})
    .add_local_dir(str(SRC), remote_path="/pkg")
)


def _load_cfg(config_yaml: str) -> dict[str, Any]:
    return yaml.safe_load(config_yaml)


def _variant_map(cfg: dict[str, Any]) -> dict[str, dict[str, Any] | None]:
    out: dict[str, dict[str, Any] | None] = {"WT": None}
    for v in cfg.get("variants", []):
        out[v["id"]] = v
    return out


@app.function(
    image=classical_image,
    volumes={RESULTS: results_volume},
    timeout=60 * MINUTES,
    cpu=4,
    memory=16384,
)
def prepare_frame(config_yaml: str, pdb_bytes: bytes, variant_id: str) -> dict[str, Any]:
    from io_layout import (
        ensure_tree,
        frame_pdb,
        manifest_path,
        run_root,
        structure_pdb,
    )
    from prepare import prepare_variant

    cfg = _load_cfg(config_yaml)
    root = run_root(cfg, RESULTS)
    ensure_tree(root)

    src_pdb = root / "structures" / "_input.pdb"
    if not src_pdb.exists():
        src_pdb.write_bytes(pdb_bytes)

    vmap = _variant_map(cfg)
    if variant_id not in vmap:
        raise KeyError(variant_id)
    out_struct = structure_pdb(root, variant_id)
    manifest = prepare_variant(src_pdb, out_struct, cfg, vmap[variant_id])
    # Step 1: single frame copy (no MD).
    frame = frame_pdb(root, variant_id, "0")
    frame.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out_struct, frame)
    man = {
        **manifest,
        "frame": "0",
        "versions": CLASSICAL_PIN,
    }
    mp = manifest_path(root, variant_id, "0")
    mp.parent.mkdir(parents=True, exist_ok=True)
    mp.write_text(yaml.dump(man), encoding="utf-8")
    results_volume.commit()
    return man


@app.function(
    image=classical_image,
    volumes={RESULTS: results_volume},
    timeout=30 * MINUTES,
    cpu=2,
    memory=8192,
)
def protonate(config_yaml: str, variant_id: str, frame: str = "0") -> dict[str, Any]:
    from io_layout import frame_pdb, prot_dir, run_root
    from protonate import run_pdb2pqr

    cfg = _load_cfg(config_yaml)
    root = run_root(cfg, RESULTS)
    pdb_in = frame_pdb(root, variant_id, frame)
    out = prot_dir(root, variant_id, frame)
    result = run_pdb2pqr(pdb_in, out, cfg)
    results_volume.commit()
    return result


@app.function(
    image=classical_image,
    volumes={RESULTS: results_volume},
    timeout=60 * MINUTES,
    cpu=4,
    memory=16384,
)
def apbs_triple(config_yaml: str, variant_id: str, frame: str = "0") -> dict[str, Any]:
    from io_layout import apbs_dir, prot_dir, run_root
    from apbs_run import run_apbs_triple

    cfg = _load_cfg(config_yaml)
    root = run_root(cfg, RESULTS)
    pqr = Path(prot_dir(root, variant_id, frame) / "structure.pqr")
    out = apbs_dir(root, variant_id, frame)
    result = run_apbs_triple(pqr, out, cfg)
    results_volume.commit()
    return result


@app.function(
    image=classical_image,
    volumes={RESULTS: results_volume},
    timeout=15 * MINUTES,
    cpu=2,
    memory=4096,
)
def cluster_cut(config_yaml: str, variant_id: str, frame: str = "0") -> dict[str, Any]:
    from io_layout import cluster_dir, prot_dir, run_root
    from cluster import cut_cluster

    cfg = _load_cfg(config_yaml)
    root = run_root(cfg, RESULTS)
    vmap = _variant_map(cfg)
    prot = prot_dir(root, variant_id, frame)
    # Prefer protonated PDB; fall back to frame PDB path recorded by pdb2pqr.
    pdb_path = Path(prot / "structure_protonated.pdb")
    if not pdb_path.exists():
        from io_layout import frame_pdb

        pdb_path = frame_pdb(root, variant_id, frame)
    pqr_path = prot / "structure.pqr"
    out = cluster_dir(root, variant_id, frame)
    sub = vmap[variant_id]
    substituted = None
    if sub is not None:
        substituted = {
            "chain": sub["chain"],
            "resid": sub["resid"],
            "from": sub["from_resname"],
            "to": sub["to_resname"],
        }
    result = cut_cluster(pdb_path, pqr_path, out, cfg, variant_id, substituted)
    results_volume.commit()
    return result


@app.function(
    image=aimnet_image,
    volumes={RESULTS: results_volume},
    gpu="T4",
    timeout=30 * MINUTES,
    memory=8192,
)
def aimnet_energy(config_yaml: str, variant_id: str, frame: str = "0") -> dict[str, Any]:
    from io_layout import aimnet_dir, cluster_dir, run_root
    from aimnet_run import run_aimnet

    cfg = _load_cfg(config_yaml)
    root = run_root(cfg, RESULTS)
    cdir = cluster_dir(root, variant_id, frame)
    manifest = json.loads((cdir / "cluster_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("metal_check", {}).get("quantum_refused"):
        out = {
            "status": "refused",
            "metal_check": manifest["metal_check"],
        }
        adir = aimnet_dir(root, variant_id, frame)
        (adir / "aimnet_result.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
        results_volume.commit()
        return out

    result = run_aimnet(
        Path(manifest["cluster_pdb"]),
        aimnet_dir(root, variant_id, frame),
        cfg,
        int(manifest["net_charge_int"]),
        protein_chains=set(cfg["structure"]["protein_chains"]),
        partner_chains=set(cfg["structure"]["partner_chains"]),
    )
    results_volume.commit()
    return result


@app.function(
    image=classical_image,
    volumes={RESULTS: results_volume},
    timeout=15 * MINUTES,
    cpu=2,
    memory=4096,
)
def aggregate_results(config_yaml: str, stage: dict[str, Any]) -> dict[str, Any]:
    from io_layout import apbs_dir, results_dir, run_root
    from aggregate import (
        make_difference_dx,
        rows_from_stage_results,
        two_row_summary,
        write_long_table,
        write_versions,
    )

    cfg = _load_cfg(config_yaml)
    root = run_root(cfg, RESULTS)
    rows: list[dict[str, Any]] = []
    for variant_id, payload in stage["variants"].items():
        rows.extend(
            rows_from_stage_results(
                variant_id,
                "0",
                payload["apbs"],
                payload.get("aimnet"),
                payload.get("manifest", {}),
            )
        )
    rdir = results_dir(root)
    long_path = rdir / "long_table.csv"
    summary_path = rdir / "summary_two_row.csv"
    df = write_long_table(rows, long_path)
    summary = two_row_summary(df, wt_id="WT", out_path=summary_path)

    # WT − variant complex dx for the first non-WT variant.
    non_wt = [v for v in stage["variants"] if v != "WT"]
    dx_out = None
    if non_wt:
        vid = non_wt[0]
        wt_dx = Path(apbs_dir(root, "WT", "0") / "complex.dx")
        var_dx = Path(apbs_dir(root, vid, "0") / "complex.dx")
        dx_out = rdir / f"WT_minus_{vid}_complex.dx"
        make_difference_dx(wt_dx, var_dx, dx_out)

    write_versions(
        rdir / "versions.json",
        {"classical": CLASSICAL_PIN, "aimnet": AIMNET_PIN},
    )
    # Protonation mismatch flags if provided.
    flags = stage.get("protonation_flags", [])
    (rdir / "protonation_flags.json").write_text(
        json.dumps(flags, indent=2), encoding="utf-8"
    )
    results_volume.commit()

    return {
        "long_table": str(long_path),
        "summary": str(summary_path),
        "summary_records": summary.to_dict(orient="records"),
        "dx_diff": str(dx_out) if dx_out else None,
        "long_preview": df.head(20).to_dict(orient="records"),
    }


@app.function(
    image=classical_image,
    volumes={RESULTS: results_volume},
    timeout=5 * MINUTES,
)
def download_results(config_yaml: str) -> dict[str, bytes]:
    cfg = _load_cfg(config_yaml)
    from io_layout import results_dir, run_root

    root = run_root(cfg, RESULTS)
    rdir = results_dir(root)
    out: dict[str, bytes] = {}
    for path in rdir.rglob("*"):
        if path.is_file():
            rel = str(path.relative_to(rdir)).replace("\\", "/")
            out[rel] = path.read_bytes()
    return out


@app.local_entrypoint()
def main(config: Optional[str] = None):
    """Vertical slice: WT + one variant, one frame, no MD."""
    import sys

    sys.path.insert(0, str(SRC))
    from protonate import flag_protonation_mismatches

    config_path = Path(config) if config else HERE / "config.example.yaml"
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config_yaml = config_path.read_text(encoding="utf-8")

    pdb_path = Path(cfg["structure"]["pdb"])
    if not pdb_path.is_absolute():
        pdb_path = (REPO_ROOT / pdb_path).resolve()
    pdb_bytes = pdb_path.read_bytes()

    variant_ids = ["WT"] + [v["id"] for v in cfg.get("variants", [])]
    print(f"vertical slice variants={variant_ids} config={config_path}")

    # 1) Prepare
    manifests = {}
    for vid in variant_ids:
        print(f"prepare {vid}")
        manifests[vid] = prepare_frame.remote(config_yaml, pdb_bytes, vid)

    # 2) Protonate (map)
    print("protonate map")
    prot_results = list(
        protonate.starmap(  # type: ignore[attr-defined]
            [(config_yaml, vid, "0") for vid in variant_ids]
        )
    )
    prot_by = dict(zip(variant_ids, prot_results))

    flags = flag_protonation_mismatches(
        {vid: prot_by[vid].get("pka_states", []) for vid in variant_ids},
        cfg,
        cfg["variants"][0]["chain"],
        int(cfg["variants"][0]["resid"]),
    )

    # 3) APBS map
    print("apbs map")
    apbs_results = list(
        apbs_triple.starmap(  # type: ignore[attr-defined]
            [(config_yaml, vid, "0") for vid in variant_ids]
        )
    )
    apbs_by = dict(zip(variant_ids, apbs_results))

    # 4) Cluster map
    print("cluster map")
    cluster_results = list(
        cluster_cut.starmap(  # type: ignore[attr-defined]
            [(config_yaml, vid, "0") for vid in variant_ids]
        )
    )
    cluster_by = dict(zip(variant_ids, cluster_results))

    # 5) AIMNet2 (sequential; small GPU jobs)
    aimnet_by = {}
    for vid in variant_ids:
        print(f"aimnet {vid}")
        aimnet_by[vid] = aimnet_energy.remote(config_yaml, vid, "0")

    stage = {
        "variants": {
            vid: {
                "manifest": manifests[vid],
                "protonate": prot_by[vid],
                "apbs": apbs_by[vid],
                "cluster": cluster_by[vid],
                "aimnet": aimnet_by[vid],
            }
            for vid in variant_ids
        },
        "protonation_flags": flags,
    }
    print("aggregate")
    agg = aggregate_results.remote(config_yaml, stage)
    print("summary:", agg.get("summary_records"))

    local_out = HERE / "results"
    local_out.mkdir(parents=True, exist_ok=True)
    files = download_results.remote(config_yaml)
    for rel, data in files.items():
        dest = local_out / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        print(f"wrote {dest}")

    print("done vertical slice")
