# electrostatics-residue-pipeline

## Binding-partner classifier

`binding_partner/` decides, from solved structures only and by distance, whether one residue
contacts another protein, the same protein (another subunit), a ligand, an ion, or nothing. The
design and the agreed rules are in [`docs/binding-partner-pipeline-talk.md`](docs/binding-partner-pipeline-talk.md).

```bash
pip install -e '.[dev]'

# Residue GLU 73 of chain A in the reference file, classified across a list of structures.
python -m binding_partner --ref ref.cif --chain A --resname ARG --resnum 83 \
    --biounit-glob '*_3mer*' --out out ref.cif other1.cif other2.pdb

pytest                      # synthetic structures, no data files needed
BRS_PDB=path/to/1BRS.pdb pytest tests/test_brs_demo.py
```

Outputs in `--out`: `files.csv` (per-file status), `evidence.csv` (one row per file, placed copy
and partner), `consensus.csv` (one row per label), `summary.json`. Thresholds and component lists
are in [`binding_partner/config.example.yaml`](binding_partner/config.example.yaml); pass a copy
with `--config`.
