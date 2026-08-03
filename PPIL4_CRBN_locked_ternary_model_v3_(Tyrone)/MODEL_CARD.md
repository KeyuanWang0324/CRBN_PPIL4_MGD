# PPIL4–CRBN Locked Ternary Docking Model v3 (Tyrone)

## Identity

- Owner: Tyrone
- Protocol version: `ppil4-lock-v3-45ee93aea3ab`
- Purpose: rank 30 specified molecular-glue candidates using a reproducible
  PPIL4–CRBN–ligand HADDOCK3 ternary-docking workflow.
- Positive control: FPFT-2216.
- Negative controls: the three structures in
  `controls/negative_controls.csv`, including FPFT-2216 without glutarimide.

## Packaged contents

The folder contains the locked receptor and 9DWV reference structures, fixed
restraints, configuration files, the exact 30-candidate queue, and all four
prepared ligand files for every candidate. `MANIFEST.sha256` protects the
packaged artifacts. Large generated docking runs are deliberately excluded
from GitHub.

## Scoring and filters

All successful minimized models are sorted by HADDOCK score. `mean_best_10` is
the arithmetic mean of the ten lowest scores. The intended candidate filter is
DockQ >= 0.45 and mean_best_10 < -475.3. Detailed ranking and tie rules are in
`RANKING_RULES.md`.

## Validation status

Status: **validation pending; candidate gate closed**.

The completed six-core FPFT-2216 trial produced 172 successful minimized
models, mean_best_10 = -115.719300, best HADDOCK = -125.639000, and DockQ =
0.042000 for the best-energy model. The maximum DockQ observed across that run
was 0.409, so no model met the intended DockQ >= 0.45 threshold. This result is
recorded for provenance but does not constitute the independent one-core
golden-reference reproduction required by `PROTOCOL_LOCK.md`.

Do not run or rank the candidate queue until an independent FPFT-2216 run is
approved in `REFERENCE_FPFT.json` with the same protocol version.

## Entry point

From the repository root:

```bash
python3 'PPIL4_CRBN_locked_ternary_model_v3_(Tyrone)/11_run_locked_ppil4_ternary_campaign_(Tyrone).py' --prepare
```

After the validation gate is formally opened, add `--run` to execute the
prepared queue. Generated runs default to `locked_campaign_runs/`, which is
excluded from Git.
