# electrostatics-residue-pipeline

Config-driven residue-substitution electrostatics pipeline. Solvers stay upstream
(PDB2PQR, APBS, PDBFixer, OpenMM, AIMNet2). This repo owns configuration, file
layout, charge/grid checks, and the merged table. GPU steps run on Modal.

## Demo system

Public PDB **1BRS** (barnase–barstar). Biology lives only in
[`electrostatics/config.example.yaml`](electrostatics/config.example.yaml).

## Step 2 (classical scale) — current default

Same prepare → PDB2PQR → APBS path over the configured variant list. One
minimized frame per variant. No MD and no AIMNet2. Emits a ranked
`summary_classical.csv` (WT-referenced `ΔΔG_elec`).

```bash
# local fixture checks
pip install -r electrostatics/requirements-dev.txt
python -m pytest electrostatics/tests

# classical scale on Modal (default mode)
python -m modal run electrostatics/app.py

# step-1 vertical slice (adds cluster + AIMNet2)
python -m modal run electrostatics/app.py --mode vertical_slice
```

Copy `electrostatics/config.example.yaml` to `electrostatics/config.yaml` when
pointing at your own complex. Heavy products land on Modal Volume
`electrostatics-results`; downloaded summaries go to `electrostatics/results/`.

## Build order

1. Vertical slice (WT + one variant + AIMNet2) — done
2. Classical scale (full variant list, classical only) — this milestone
3. Noise floor (OpenMM sampling) — not started
4. Quantum at the same frames — not started
5. Readout / ranks across frames — not started
