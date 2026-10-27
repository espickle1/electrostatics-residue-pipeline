"""AIMNet2 single-point charges and interaction energy on a cluster."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from cluster import assert_charge_handoff, refuse_if_unsupported_elements


def _read_xyz_from_pdb(pdb_path: Path) -> tuple[list[str], list[list[float]]]:
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        if len(line) >= 78 and line[76:78].strip():
            raw = line[76:78].strip()
            elem = raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
        else:
            letters = "".join(c for c in name if c.isalpha())
            if letters.upper().startswith("CL"):
                elem = "Cl"
            elif letters.upper().startswith("BR"):
                elem = "Br"
            else:
                elem = letters[0].upper() if letters else "X"
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        symbols.append(elem)
        coords.append([x, y, z])
    return symbols, coords


def run_aimnet(
    cluster_pdb: Path,
    out_dir: Path,
    cfg: dict[str, Any],
    net_charge_int: int,
    protein_chains: set[str] | None = None,
    partner_chains: set[str] | None = None,
) -> dict[str, Any]:
    """Run AIMNet2; fail if charge handoff or unsupported elements."""
    out_dir.mkdir(parents=True, exist_ok=True)
    symbols, coords = _read_xyz_from_pdb(cluster_pdb)
    refuse = refuse_if_unsupported_elements(set(symbols))
    if refuse["quantum_refused"]:
        result = {
            "status": "refused",
            "metal_check": refuse,
            "energy_eV": None,
            "charges": None,
        }
        (out_dir / "aimnet_result.json").write_text(
            json.dumps(result, indent=2), encoding="utf-8"
        )
        return result

    import numpy as np
    from ase import Atoms
    from aimnet.calculators import AIMNet2ASE

    model = cfg.get("aimnet", {}).get("model", "aimnet2")
    atoms = Atoms(symbols=symbols, positions=coords)
    atoms.calc = AIMNet2ASE(model, charge=int(net_charge_int))
    energy = float(atoms.get_potential_energy())
    charges = atoms.calc.results.get("charges")
    if charges is None:
        raise RuntimeError("AIMNet2 did not return atomic charges")
    charges_list = [float(c) for c in np.asarray(charges)]
    qsum = float(sum(charges_list))
    assert_charge_handoff(float(net_charge_int), qsum, tol=0.25)

    # Interaction proxy: E(complex cluster) only for step 1; split energies optional.
    interaction = None
    if protein_chains and partner_chains:
        interaction = _interaction_energy(
            cluster_pdb, model, net_charge_int, protein_chains, partner_chains
        )

    result = {
        "status": "ok",
        "model": model,
        "net_charge": int(net_charge_int),
        "charge_sum": qsum,
        "energy_eV": energy,
        "interaction_eV": interaction,
        "n_atoms": len(symbols),
        "metal_check": refuse,
        "charges": charges_list,
    }
    (out_dir / "aimnet_result.json").write_text(
        json.dumps(
            {k: v for k, v in result.items() if k != "charges"}
            | {"charges_head": charges_list[:20]},
            indent=2,
        ),
        encoding="utf-8",
    )
    return result


def _interaction_energy(
    cluster_pdb: Path,
    model: str,
    net_charge_int: int,
    protein_chains: set[str],
    partner_chains: set[str],
) -> float:
    """ΔE_int ≈ E(all) − E(protein atoms) − E(partner atoms) with charge split by PQR-like neutrality heuristic."""
    import numpy as np
    from ase import Atoms
    from aimnet.calculators import AIMNet2ASE

    # Build index masks from PDB chain IDs.
    symbols_all: list[str] = []
    coords_all: list[list[float]] = []
    chains: list[str] = []
    for line in cluster_pdb.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21].strip()
        name = line[12:16].strip()
        if len(line) >= 78 and line[76:78].strip():
            raw = line[76:78].strip()
            elem = raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
        else:
            letters = "".join(c for c in name if c.isalpha())
            elem = letters[0].upper() if letters else "C"
        symbols_all.append(elem)
        coords_all.append(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
        chains.append(chain)

    def energy_subset(mask: list[bool], charge: int) -> float:
        sym = [s for s, m in zip(symbols_all, mask) if m]
        xyz = [c for c, m in zip(coords_all, mask) if m]
        if not sym:
            return 0.0
        atoms = Atoms(symbols=sym, positions=xyz)
        atoms.calc = AIMNet2ASE(model, charge=charge)
        return float(atoms.get_potential_energy())

    mask_all = [True] * len(symbols_all)
    mask_prot = [c in protein_chains for c in chains]
    mask_part = [c in partner_chains for c in chains]

    # Neutral fragments when unknown; complex uses handoff charge.
    e_all = energy_subset(mask_all, net_charge_int)
    e_prot = energy_subset(mask_prot, 0)
    e_part = energy_subset(mask_part, 0)
    _ = np  # kept for potential charge arrays later
    return e_all - e_prot - e_part
