"""Long-format table and step-1 two-row summary."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from apbs_run import ddg_elec, subtract_dx


def rows_from_stage_results(
    variant: str,
    frame: str,
    apbs: dict[str, Any],
    aimnet: dict[str, Any] | None,
    manifest: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    for mol, e in apbs["elec_kJ_mol"].items():
        rows.append(
            {
                "variant": variant,
                "frame": frame,
                "quantity": f"ELEC_{mol}",
                "value": e,
                "unit": "kJ/mol",
            }
        )
    rows.append(
        {
            "variant": variant,
            "frame": frame,
            "quantity": "ddG_binding_elec",
            "value": apbs["ddg_binding_kJ_mol"],
            "unit": "kJ/mol",
        }
    )
    if aimnet and aimnet.get("status") == "ok":
        rows.append(
            {
                "variant": variant,
                "frame": frame,
                "quantity": "E_aimnet",
                "value": aimnet["energy_eV"],
                "unit": "eV",
            }
        )
        if aimnet.get("interaction_eV") is not None:
            rows.append(
                {
                    "variant": variant,
                    "frame": frame,
                    "quantity": "dE_int",
                    "value": aimnet["interaction_eV"],
                    "unit": "eV",
                }
            )
    elif aimnet and aimnet.get("status") == "refused":
        rows.append(
            {
                "variant": variant,
                "frame": frame,
                "quantity": "aimnet_refused",
                "value": 1.0,
                "unit": "flag",
            }
        )
    rows.append(
        {
            "variant": variant,
            "frame": frame,
            "quantity": "rmsd",
            "value": manifest.get("rmsd"),
            "unit": "Angstrom",
        }
    )
    return rows


def write_long_table(rows: list[dict[str, Any]], path: Path) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


def two_row_summary(
    long_df: pd.DataFrame,
    wt_id: str = "WT",
    out_path: Path | None = None,
) -> pd.DataFrame:
    """One row per variant with ΔΔG_elec and ΔE_int vs WT."""
    bind = long_df[long_df["quantity"] == "ddG_binding_elec"].set_index("variant")
    if wt_id not in bind.index:
        raise KeyError(f"WT row {wt_id!r} missing from long table")
    wt_bind = float(bind.loc[wt_id, "value"])

    de = long_df[long_df["quantity"] == "dE_int"].set_index("variant")
    wt_de = float(de.loc[wt_id, "value"]) if wt_id in de.index else None

    summary_rows = []
    for variant in bind.index.unique():
        v_bind = float(bind.loc[variant, "value"])
        row = {
            "variant": variant,
            "ddG_binding_elec_kJ_mol": v_bind,
            "dddG_elec_kJ_mol": ddg_elec(v_bind, wt_bind),
        }
        if wt_de is not None and variant in de.index:
            v_de = float(de.loc[variant, "value"])
            row["dE_int_eV"] = v_de
            row["ddE_int_eV"] = v_de - wt_de
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_path, index=False)
    return summary


def write_versions(path: Path, versions: dict[str, str]) -> None:
    path.write_text(json.dumps(versions, indent=2), encoding="utf-8")


def make_difference_dx(
    wt_complex_dx: Path,
    variant_complex_dx: Path,
    out_dx: Path,
) -> Path:
    return subtract_dx(wt_complex_dx, variant_complex_dx, out_dx)
