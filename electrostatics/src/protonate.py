"""PDB2PQR / PROPKA wrapper and pKa parse helpers."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Any


NET_CHARGE_RE = re.compile(
    r"(?:Total\s+charge|Net\s+charge)\s*[:=]?\s*([+-]?\d+\.?\d*)",
    re.IGNORECASE,
)
PKA_LINE_RE = re.compile(
    r"^\s*([A-Z]{3})\s+([A-Z0-9]?)\s+(\d+)\s+([+-]?\d+\.?\d*)",
)


def run_pdb2pqr(
    pdb_in: Path,
    out_dir: Path,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    """Run pdb2pqr30 once; write .pqr, .pdb, .log under out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "structure"
    pqr = out_dir / f"{stem}.pqr"
    pdb_out = out_dir / f"{stem}_protonated.pdb"
    log_path = out_dir / f"{stem}.log"
    ff = cfg.get("apbs", {}).get("forcefield", "PARSE")
    ph = float(cfg.get("pH", 7.0))

    cmd = [
        "pdb2pqr30",
        f"--ff={ff}",
        f"--with-ph={ph}",
        "--drop-water",
        "--keep-chain",
        f"--pdb-output={pdb_out}",
        "--log-level=INFO",
        str(pdb_in),
        str(pqr),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log_path.write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or ""),
        encoding="utf-8",
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"pdb2pqr30 failed ({proc.returncode}): {proc.stderr[-2000:]}"
        )
    if not pqr.exists():
        raise FileNotFoundError(f"Expected PQR missing: {pqr}")

    net_charge = parse_net_charge(log_path.read_text(encoding="utf-8"))
    if net_charge is None:
        net_charge = sum_pqr_charges(pqr)

    pka_path = out_dir / f"{stem}.pka"
    # Newer pdb2pqr may write propka output next to the pqr or only in the log.
    pka_states = parse_pka_table(log_path.read_text(encoding="utf-8"))
    if pka_path.exists():
        pka_states = parse_pka_table(pka_path.read_text(encoding="utf-8")) or pka_states

    return {
        "pqr": str(pqr),
        "pdb": str(pdb_out) if pdb_out.exists() else str(pdb_in),
        "log": str(log_path),
        "net_charge": float(net_charge),
        "pka_states": pka_states,
    }


def parse_net_charge(text: str) -> float | None:
    matches = NET_CHARGE_RE.findall(text)
    if not matches:
        return None
    return float(matches[-1])


def sum_pqr_charges(pqr_path: Path) -> float:
    total = 0.0
    for line in pqr_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        # PQR: charge is typically field after xyz / occupancy-like columns.
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            total += float(parts[-2])
        except ValueError:
            try:
                total += float(parts[8])
            except ValueError:
                continue
    return total


def parse_pka_table(text: str) -> list[dict[str, Any]]:
    """Best-effort PROPKA / pdb2pqr pKa rows."""
    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        m = PKA_LINE_RE.match(line)
        if not m:
            continue
        resname, chain, resid, pka = m.groups()
        rows.append(
            {
                "resname": resname,
                "chain": chain.strip() or "?",
                "resid": int(resid),
                "pka": float(pka),
            }
        )
    return rows


def flag_protonation_mismatches(
    states_by_variant: dict[str, list[dict[str, Any]]],
    cfg: dict[str, Any],
    site_chain: str,
    site_resid: int,
) -> list[dict[str, Any]]:
    """Flag titratable groups within radius whose state differs across variants.

    Step-1 heuristic: compare reported pKa relative to working pH (protonated if pKa > pH).
    """
    ph = float(cfg.get("pH", 7.0))
    radius = float(cfg.get("protonation", {}).get("flag_radius_A", 12.0))
    titratable = set(cfg.get("protonation", {}).get("titratable", []))
    # Without coordinates here, flag any titratable with differing binary state.
    # Coordinate filtering is applied in cluster/aggregate when manifests merge.
    _ = (radius, site_chain, site_resid)
    keys: dict[tuple, dict[str, str]] = {}
    for variant, rows in states_by_variant.items():
        for row in rows:
            if row["resname"] not in titratable:
                continue
            key = (row["chain"], row["resid"], row["resname"])
            state = "protonated" if row["pka"] > ph else "deprotonated"
            keys.setdefault(key, {})[variant] = state
    flags = []
    for key, by_var in keys.items():
        if len(set(by_var.values())) > 1:
            flags.append(
                {
                    "chain": key[0],
                    "resid": key[1],
                    "resname": key[2],
                    "states": by_var,
                    "note": "keep_both_states",
                }
            )
    return flags
