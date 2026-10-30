"""Hand-built PDB text for tests: straight peptides plus hetero groups at chosen coordinates.

A residue starting at x0 has N (x0), CA (x0+1.46), C (x0+2.98), O (x0+2.98, y0+1.2) and
CB (x0+1.46, y0-1.5); the next residue starts 4.31 A further along x. With a peptide of
chain A at y0=0, the second residue has CA (5.77, 0, 0) and CB (5.77, -1.5, 0).
"""

from __future__ import annotations

from pathlib import Path

STEP = 4.31

A_SEQ = ["ALA", "GLU", "LYS", "LEU", "SER", "GLY", "VAL", "TRP"]
B_SEQ = ["TRP", "PHE", "TYR", "HIS", "ARG", "ASP", "CYS", "MET"]  # no overlap with A_SEQ


class PdbBuilder:
    def __init__(self, header: list[str] | None = None):
        self.lines: list[str] = list(header or [])
        self.serial = 0

    def atom(self, name, resname, chain, resnum, xyz, element, *, het=False, occ=1.0, alt=" "):
        self.serial += 1
        label = f" {name:<3}" if len(element) == 1 else f"{name:<4}"
        record = "HETATM" if het else "ATOM"
        self.lines.append(
            f"{record:<6}{self.serial:>5} {label}{alt}{resname:>3} {chain}{resnum:>4}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{occ:6.2f}{10.0:6.2f}          {element:>2}"
        )
        return self

    def peptide(self, chain, residues, x0=0.0, y0=0.0, z0=0.0, start=1, het=(), no_cb=()):
        x = x0
        for i, name in enumerate(residues):
            num = start + i
            kw = {"het": name in het}
            self.atom("N", name, chain, num, (x, y0, z0), "N", **kw)
            self.atom("CA", name, chain, num, (x + 1.46, y0, z0), "C", **kw)
            self.atom("C", name, chain, num, (x + 2.98, y0, z0), "C", **kw)
            self.atom("O", name, chain, num, (x + 2.98, y0 + 1.2, z0), "O", **kw)
            if name != "GLY" and num not in no_cb:
                self.atom("CB", name, chain, num, (x + 1.46, y0 - 1.5, z0), "C", **kw)
            x += STEP
        return self

    def group(self, chain, resname, resnum, atoms):
        """A hetero group; atoms is a list of (name, element, (x, y, z))."""
        for name, element, xyz in atoms:
            self.atom(name, resname, chain, resnum, xyz, element, het=True)
        return self

    def text(self) -> str:
        return "\n".join(self.lines) + "\nEND\n"

    def write(self, path: Path) -> Path:
        path.write_text(self.text(), encoding="utf-8")
        return path


def cryst1(spacegroup: str = "P 2", cell: float = 50.0, z: int = 2) -> str:
    return (
        f"CRYST1{cell:9.3f}{cell:9.3f}{cell:9.3f}{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} "
        f"{spacegroup:<11}{z:>4}"
    )


def biomt(chains: str) -> list[str]:
    """REMARK 350 header: biological assembly 1 is the listed chains, identity operator."""
    return [
        "REMARK 350 BIOMOLECULE: 1",
        f"REMARK 350 APPLY THE FOLLOWING TO CHAINS: {chains}",
        "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000",
        "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000",
        "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000",
    ]


# Three heavy atoms 3-4 A below the CB of residue 2 of chain A: a qualifying small molecule.
LIGAND_ATOMS = [
    ("C1", "C", (5.0, -4.5, 0.0)),
    ("C2", "C", (6.5, -4.5, 0.0)),
    ("O1", "O", (5.8, -5.5, 0.0)),
]
