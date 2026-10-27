"""One-shot helper: WT 8 A neighborhood of the substituted site A:GLU73."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

pdb_lines = Path(__file__).with_name("1BRS.pdb").read_text().splitlines()
atoms: list[tuple] = []
for line in pdb_lines:
    if not line.startswith("ATOM"):
        continue
    chain = line[21]
    if chain not in ("A", "D"):
        continue
    resname = line[17:20].strip()
    resseq = int(line[22:26])
    name = line[12:16].strip()
    alt = line[16]
    if alt not in (" ", "A"):
        continue
    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
    atoms.append((chain, resseq, resname, name, x, y, z))

by_res: dict[tuple, list] = defaultdict(list)
for a in atoms:
    if a[3].startswith("H"):
        continue
    by_res[(a[0], a[1], a[2])].append(a)

site = by_res[("A", 73, "GLU")]
cutoff = 8.0
cluster: set[tuple] = set()
for k, v in by_res.items():
    hit = False
    for a in v:
        for s in site:
            dx = a[4] - s[4]
            dy = a[5] - s[5]
            dz = a[6] - s[6]
            if dx * dx + dy * dy + dz * dz <= cutoff * cutoff:
                hit = True
                break
        if hit:
            break
    if hit:
        cluster.add(k)

n_heavy = sum(len(by_res[k]) for k in cluster)
ca = next(a for a in site if a[3] == "CA")
items = sorted(cluster, key=lambda k: (k[0], k[1]))
print("residues", len(items), "heavy_atoms", n_heavy)
print("gcent", round(ca[4], 3), round(ca[5], 3), round(ca[6], 3))
for k in items:
    print(f"  - {k[0]}:{k[2]}{k[1]}")
