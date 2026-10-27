# electrostatics-residue-pipeline

Config-driven residue-substitution electrostatics pipeline. Solvers stay upstream
(PDB2PQR, APBS, PDBFixer, OpenMM, AIMNet2). This repo owns configuration, file
layout, charge/grid checks, and the merged table. GPU steps run on Modal.

## Step 1 (vertical slice)

Wild type plus one demo variant on public PDB **1BRS** (barnase–barstar), one
minimized frame each, no MD. Emits a two-row summary and a WT−variant `.dx`.

```bash
# local fixture checks
pip install -r electrostatics/requirements-dev.txt
pytest electrostatics/tests

# end-to-end on Modal (requires `modal` auth)
modal run electrostatics/app.py
```

Copy [`electrostatics/config.example.yaml`](electrostatics/config.example.yaml)
to `electrostatics/config.yaml` and point `structure.pdb` at your complex when
ready. Heavy products land on Modal Volume `electrostatics-results`; downloaded
summaries go to `electrostatics/results/`.
