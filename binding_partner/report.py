"""CSV / JSON output for the classifier."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from .analyze import FileResult
from .vote import Summary

FILE_FIELDS = ["file", "entry_id", "status", "detail", "model", "n_copies"]
EVIDENCE_FIELDS = [
    "file", "entry_id", "model", "site_chain", "site_resname", "site_resnum", "ref_identity",
    "relative_sasa", "label", "partner_id", "component", "entity_id", "n_pairs", "min_dist",
    "in_assembly",
]
CONSENSUS_FIELDS = [
    "label", "component", "n_entries", "n_eligible", "fraction", "status", "confidence",
    "partners", "min_dist", "max_pairs", "entries",
]


def _round(value: float | None, digits: int = 3) -> float | str:
    return "" if value is None else round(value, digits)


def file_rows(results: list[FileResult]) -> list[dict]:
    return [
        {
            "file": r.path, "entry_id": r.entry_id, "status": r.status, "detail": r.detail,
            "model": r.model_kind, "n_copies": len(r.sites),
        }
        for r in results
    ]


def evidence_rows(results: list[FileResult]) -> list[dict]:
    rows: list[dict] = []
    for r in results:
        base = {"file": r.path, "entry_id": r.entry_id, "model": r.model_kind}
        for s in r.sites:
            site = {
                **base, "site_chain": s.chain, "site_resname": s.resname,
                "site_resnum": s.resnum, "ref_identity": _round(s.identity),
                "relative_sasa": _round(s.rsa), "in_assembly": True,
            }
            if not s.partners:
                rows.append({**site, "label": s.baseline})
            for p in s.partners:
                rows.append({
                    **site, "label": p.label, "partner_id": p.partner_id,
                    "component": p.component, "entity_id": p.entity_id,
                    "n_pairs": p.n_pairs, "min_dist": _round(p.min_dist, 2),
                })
        for p in r.crystal_only:
            rows.append({
                **base, "label": p.label, "partner_id": p.partner_id, "component": p.component,
                "entity_id": p.entity_id, "n_pairs": p.n_pairs,
                "min_dist": _round(p.min_dist, 2), "in_assembly": False,
            })
    return rows


def consensus_rows(summary: Summary) -> list[dict]:
    return [
        {**asdict(r), "min_dist": _round(r.min_dist, 2),
         "max_pairs": "" if r.max_pairs is None else r.max_pairs}
        for r in summary.rows
    ]


def _write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(out_dir: str | Path, results: list[FileResult], summary: Summary) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    _write_csv(out / "files.csv", FILE_FIELDS, file_rows(results))
    _write_csv(out / "evidence.csv", EVIDENCE_FIELDS, evidence_rows(results))
    _write_csv(out / "consensus.csv", CONSENSUS_FIELDS, consensus_rows(summary))
    payload = {k: v for k, v in asdict(summary).items() if k != "rows"}
    (out / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
