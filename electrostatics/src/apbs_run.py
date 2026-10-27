"""APBS input generation, grid identity checks, ELEC parse, dx differencing."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


ELEC_ENERGY_RE = re.compile(
    r"(?:Global net ELEC energy|Total electrostatic energy|Electostatic energy|"
    r"Electrostatic energy)\s*=\s*([+-]?\d+\.?\d*(?:[Ee][+-]?\d+)?)"
    r"(?:\s*kJ/?mol)?",
    re.IGNORECASE,
)

GRID_KEYS = ("dime", "glen", "gcent")


def grid_block_from_config(cfg: dict[str, Any]) -> dict[str, Any]:
    g = cfg["grid"]
    return {
        "dime": [int(x) for x in g["dime"]],
        "glen": [float(x) for x in g["glen"]],
        "gcent": [float(x) for x in g["gcent"]],
    }


def format_grid_lines(grid: dict[str, Any]) -> list[str]:
    dime = " ".join(str(x) for x in grid["dime"])
    glen = " ".join(f"{x:.6f}" for x in grid["glen"])
    gcent = " ".join(f"{x:.6f}" for x in grid["gcent"])
    return [
        f"    dime {dime}",
        f"    glen {glen}",
        f"    gcent {gcent}",
    ]


def extract_grid_block(apbs_in_text: str) -> dict[str, tuple]:
    """Parse dime / glen / gcent from a written APBS input."""
    out: dict[str, tuple] = {}
    for key in GRID_KEYS:
        m = re.search(rf"^\s*{key}\s+(.+)$", apbs_in_text, re.MULTILINE)
        if not m:
            raise ValueError(f"Missing {key} in APBS input")
        parts = m.group(1).split()
        if key == "dime":
            out[key] = tuple(int(x) for x in parts[:3])
        else:
            out[key] = tuple(float(x) for x in parts[:3])
    return out


def assert_identical_grids(inputs: dict[str, str]) -> dict[str, tuple]:
    """Abort if complex / protein / partner grid blocks differ."""
    parsed = {name: extract_grid_block(text) for name, text in inputs.items()}
    names = list(parsed)
    ref = parsed[names[0]]
    for name in names[1:]:
        if parsed[name] != ref:
            raise AssertionError(
                f"APBS grid mismatch: {names[0]}={ref} vs {name}={parsed[name]}"
            )
    return ref


def _pqr_chain_id(line: str) -> str:
    """Resolve chain ID from fixed-column or whitespace PQR/PDB layouts."""
    if len(line) > 21 and line[21].strip():
        return line[21].strip()
    parts = line.split()
    # pdb2pqr often: ATOM n name res chain resid x y z charge radius
    if len(parts) >= 5 and len(parts[4]) <= 2 and not parts[4].lstrip("-").isdigit():
        return parts[4]
    # fallback: residue field may embed chain as resseq insertion
    if len(parts) >= 5 and parts[4][:1].isalpha() and parts[4][1:].isdigit():
        return parts[4][:1]
    return ""


def split_pqr_by_chains(
    pqr_path: Path,
    protein_chains: set[str],
    partner_chains: set[str],
    out_dir: Path,
) -> dict[str, Path]:
    lines = pqr_path.read_text(encoding="utf-8").splitlines()
    buckets: dict[str, list[str]] = {"complex": [], "protein": [], "partner": []}
    seen_chains: set[str] = set()
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = _pqr_chain_id(line)
        seen_chains.add(chain)
        if chain in protein_chains or chain in partner_chains:
            buckets["complex"].append(line)
        if chain in protein_chains:
            buckets["protein"].append(line)
        if chain in partner_chains:
            buckets["partner"].append(line)
    paths = {}
    for name, atom_lines in buckets.items():
        path = out_dir / f"{name}.pqr"
        path.write_text("\n".join(atom_lines) + "\n", encoding="utf-8")
        paths[name] = path
        if not atom_lines:
            sample = "\n".join(lines[:5])
            raise RuntimeError(
                f"Split PQR {name} is empty from {pqr_path}. "
                f"seen_chains={sorted(seen_chains)} "
                f"protein={sorted(protein_chains)} partner={sorted(partner_chains)}. "
                f"sample:\n{sample}"
            )
    return paths


def write_apbs_input(
    name: str,
    pqr_path: Path,
    out_in: Path,
    cfg: dict[str, Any],
) -> str:
    grid = grid_block_from_config(cfg)
    pdie = float(cfg.get("solute_dielectric", 2.0))
    sdie = float(cfg.get("solvent_dielectric", 78.5))
    temp = float(cfg.get("temperature_K", 298.15))
    ionic = float(cfg.get("ionic_strength_M", 0.15))
    grid_lines = "\n".join(format_grid_lines(grid))
    text = f"""read
    mol pqr {pqr_path.name}
end
elec name {name}
    mg-manual
{grid_lines}
    mol 1
    lpbe
    bcfl sdh
    ion charge 1 conc {ionic:.4f} radius 2.0
    ion charge -1 conc {ionic:.4f} radius 2.0
    pdie {pdie:.4f}
    sdie {sdie:.4f}
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp {temp:.2f}
    calcenergy total
    calcforce no
    write pot dx {name}
end
quit
"""
    out_in.write_text(text, encoding="utf-8")
    return text


def parse_elec_energy(log_text: str) -> float:
    matches = ELEC_ENERGY_RE.findall(log_text)
    if not matches:
        # Fallback: any line mentioning ELEC energy with a float.
        loose = re.findall(
            r"ELEC[^\n]*?([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)",
            log_text,
            flags=re.IGNORECASE,
        )
        if loose:
            return float(loose[-1])
        raise ValueError(
            "No Global net ELEC energy found in APBS log. Tail:\n"
            + "\n".join(log_text.strip().splitlines()[-40:])
        )
    return float(matches[-1])


def run_apbs_triple(
    pqr_path: Path,
    out_dir: Path,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    """Three APBS jobs with identical grids; return ELEC energies and dx paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    protein = set(cfg["structure"]["protein_chains"])
    partner = set(cfg["structure"]["partner_chains"])
    pqrs = split_pqr_by_chains(pqr_path, protein, partner, out_dir)

    written: dict[str, str] = {}
    for name, path in pqrs.items():
        in_path = out_dir / f"{name}.in"
        written[name] = write_apbs_input(name, path, in_path, cfg)

    grid = assert_identical_grids(written)

    energies: dict[str, float] = {}
    dx_paths: dict[str, Path] = {}
    for name in ("complex", "protein", "partner"):
        in_path = out_dir / f"{name}.in"
        log_path = out_dir / f"{name}.log"
        proc = subprocess.run(
            ["apbs", str(in_path.name)],
            cwd=str(out_dir),
            capture_output=True,
            text=True,
            check=False,
        )
        log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
        log_path.write_text(log_text, encoding="utf-8")
        if proc.returncode != 0:
            raise RuntimeError(f"APBS {name} failed: {proc.stderr[-2000:]}")
        energies[name] = parse_elec_energy(log_text)
        dx = out_dir / f"{name}.dx"
        # APBS may write name.dx in cwd.
        if not dx.exists():
            candidates = list(out_dir.glob(f"*{name}*.dx"))
            if candidates:
                shutil.move(str(candidates[0]), str(dx))
        if not dx.exists():
            raise FileNotFoundError(f"APBS dx missing for {name}")
        dx_paths[name] = dx

    binding = energies["complex"] - energies["protein"] - energies["partner"]
    return {
        "grid": {k: list(v) for k, v in grid.items()},
        "elec_kJ_mol": energies,
        "ddg_binding_kJ_mol": binding,
        "dx": {k: str(v) for k, v in dx_paths.items()},
    }


def ddg_elec(variant_binding: float, wt_binding: float) -> float:
    """ΔΔG_elec = binding(variant) − binding(WT)."""
    return variant_binding - wt_binding


def subtract_dx(wt_dx: Path, var_dx: Path, out_dx: Path) -> Path:
    """Write WT − variant potential difference."""
    from gridData import Grid

    g_wt = Grid(str(wt_dx))
    g_var = Grid(str(var_dx))
    if g_wt.grid.shape != g_var.grid.shape:
        raise ValueError(
            f"dx shape mismatch: {g_wt.grid.shape} vs {g_var.grid.shape}"
        )
    diff = g_wt.grid - g_var.grid
    g_out = Grid()
    g_out.grid = diff
    g_out.origin = g_wt.origin
    g_out.delta = g_wt.delta
    out_dx.parent.mkdir(parents=True, exist_ok=True)
    g_out.export(str(out_dx), file_format="dx")
    return out_dx
