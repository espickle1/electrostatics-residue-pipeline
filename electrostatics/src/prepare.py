"""PDBFixer clean / mutate / short restrained minimise."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def _restore_chain_ids(pdb_path: Path, expected_ids: list[str]) -> None:
    """Rewrite chain IDs in a PDB to match config order of appearance."""
    lines = pdb_path.read_text(encoding="utf-8").splitlines()
    seen: list[str] = []
    mapping: dict[str, str] = {}
    out: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            old = line[21] if len(line) > 21 else " "
            if old not in mapping:
                if len(seen) >= len(expected_ids):
                    raise RuntimeError(
                        f"More chains in {pdb_path} than expected {expected_ids}"
                    )
                mapping[old] = expected_ids[len(seen)]
                seen.append(old)
            new = mapping[old]
            line = line[:21] + new + line[22:]
        out.append(line)
    if set(mapping.values()) != set(expected_ids):
        raise RuntimeError(
            f"Chain restore incomplete: {mapping} vs expected {expected_ids}"
        )
    pdb_path.write_text("\n".join(out) + "\n", encoding="utf-8")


def _keep_chains(pdb_text: str, chains: set[str]) -> str:
    out: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if line[21] not in chains:
                continue
            if line[17:20] == "HOH":
                continue
            out.append(line)
        elif line.startswith("END"):
            out.append("END")
        elif line.startswith("TER"):
            if out and out[-1].startswith(("ATOM", "HETATM")):
                out.append(line)
    if not out or not out[-1].startswith("END"):
        out.append("END")
    return "\n".join(out) + "\n"


def prepare_variant(
    pdb_path: Path,
    out_pdb: Path,
    cfg: dict[str, Any],
    variant: dict[str, Any] | None,
) -> dict[str, Any]:
    """Clean PDB, optionally mutate one site, restrained-minimise, write PDB."""
    from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, unit
    from openmm.app import ForceField, Modeller, PDBFile, Simulation
    from pdbfixer import PDBFixer

    protein = set(cfg["structure"]["protein_chains"])
    partner = set(cfg["structure"]["partner_chains"])
    keep = protein | partner

    raw = Path(pdb_path).read_text(encoding="utf-8", errors="replace")
    filtered = _keep_chains(raw, keep)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    tmp_in = out_pdb.parent / f"_tmp_in_{out_pdb.stem}.pdb"
    tmp_in.write_text(filtered, encoding="utf-8")

    fixer = PDBFixer(filename=str(tmp_in))
    fixer.findMissingResidues()
    # Do not build entire missing terminal loops for the demo slice.
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)

    if variant is not None:
        # PDBFixer expects three-letter codes: ORIG-resid-NEW
        from_aa = variant["from_resname"].upper()
        to_aa = variant["to_resname"].upper()
        if from_aa not in THREE_TO_ONE or to_aa not in THREE_TO_ONE:
            raise ValueError(f"Unsupported mutation residues {from_aa}->{to_aa}")
        resid = int(variant["resid"])
        chain_id = variant["chain"]
        fixer.applyMutations([f"{from_aa}-{resid}-{to_aa}"], chain_id)

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(float(cfg.get("pH", 7.0)))

    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    system = forcefield.createSystem(modeller.topology)

    k = float(cfg.get("prepare", {}).get("restraint_k_kcal", 10.0))
    restraint = CustomExternalForce(f"{k}*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for atom in modeller.topology.atoms():
        if atom.element.symbol == "H":
            continue
        pos = modeller.positions[atom.index]
        restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])
    system.addForce(restraint)

    integrator = LangevinMiddleIntegrator(
        float(cfg.get("temperature_K", 298.15)) * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    try:
        platform = Platform.getPlatformByName("CPU")
    except Exception:
        platform = Platform.getPlatformByName("Reference")

    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    max_iter = int(cfg.get("prepare", {}).get("minimise_max_iterations", 200))
    simulation.minimizeEnergy(maxIterations=max_iter)
    positions = simulation.context.getState(getPositions=True).getPositions()

    with open(out_pdb, "w", encoding="utf-8") as fh:
        PDBFile.writeFile(simulation.topology, positions, fh)
    tmp_in.unlink(missing_ok=True)

    # PDBFixer/OpenMM may renumber chain IDs (e.g. D -> B). Restore config IDs
    # in order of appearance so APBS protein/partner splits stay correct.
    expected_ids = list(cfg["structure"]["protein_chains"]) + list(
        cfg["structure"]["partner_chains"]
    )
    _restore_chain_ids(out_pdb, expected_ids)

    manifest = {
        "variant": "WT" if variant is None else variant["id"],
        "rmsd": None,
        "rmsd_note": "step1_no_md_alignment",
        "source_pdb": str(pdb_path),
        "output_pdb": str(out_pdb),
        "mutation": None
        if variant is None
        else {
            "chain": variant["chain"],
            "resid": variant["resid"],
            "from": variant["from_resname"],
            "to": variant["to_resname"],
        },
    }
    out_pdb.with_suffix(".json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
