"""Long-format table and classical / vertical-slice summaries."""

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


def classical_summary(
    long_df: pd.DataFrame,
    wt_id: str = "WT",
    out_path: Path | None = None,
) -> pd.DataFrame:
    """One row per variant with WT-referenced ΔΔG_elec (and ΔE_int if present).

    Ranked by |ΔΔG_elec| descending so the classical-scale table is easy to scan.
    Absolute ELEC energies stay only in the long table and are not ranked.
    """
    bind = long_df[long_df["quantity"] == "ddG_binding_elec"].copy()
    if bind.empty:
        raise KeyError("ddG_binding_elec missing from long table")
    # One value per variant (step 1–2: single frame).
    bind = bind.groupby("variant", as_index=True)["value"].mean()
    if wt_id not in bind.index:
        raise KeyError(f"WT row {wt_id!r} missing from long table")
    wt_bind = float(bind.loc[wt_id])

    de = long_df[long_df["quantity"] == "dE_int"]
    de_s = (
        de.groupby("variant")["value"].mean() if not de.empty else pd.Series(dtype=float)
    )
    wt_de = float(de_s.loc[wt_id]) if wt_id in de_s.index else None

    summary_rows = []
    for variant in bind.index.unique():
        v_bind = float(bind.loc[variant])
        row: dict[str, Any] = {
            "variant": variant,
            "ddG_binding_elec_kJ_mol": v_bind,
            "dddG_elec_kJ_mol": ddg_elec(v_bind, wt_bind),
        }
        if wt_de is not None and variant in de_s.index:
            v_de = float(de_s.loc[variant])
            row["dE_int_eV"] = v_de
            row["ddE_int_eV"] = v_de - wt_de
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary["abs_dddG_elec_kJ_mol"] = summary["dddG_elec_kJ_mol"].abs()
    summary = summary.sort_values(
        by=["abs_dddG_elec_kJ_mol", "variant"],
        ascending=[False, True],
    ).reset_index(drop=True)
    summary["rank"] = range(1, len(summary) + 1)
    summary = summary.drop(columns=["abs_dddG_elec_kJ_mol"])

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_path, index=False)
    return summary


# Back-compat alias used by step-1 vertical slice callers / tests.
two_row_summary = classical_summary


def write_versions(path: Path, versions: dict[str, Any]) -> None:
    path.write_text(json.dumps(versions, indent=2), encoding="utf-8")


def make_difference_dx(
    wt_complex_dx: Path,
    variant_complex_dx: Path,
    out_dx: Path,
) -> Path:
    return subtract_dx(wt_complex_dx, variant_complex_dx, out_dx)
