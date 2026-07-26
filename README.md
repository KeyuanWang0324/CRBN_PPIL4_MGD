# CRBN_MGD

> **CRBN** (Cereblon) **M**olecular **G**lue **D**egrader project.
>
> _This README is a starting scaffold — update the description and sections below to match the project._

## Overview

Short description of what this repository is for (e.g. modeling, screening, or
analysis of cereblon-based molecular glue degraders). _TODO: fill in._

## Repository structure

```
CRBN_MGD/
├── README.md        # this file
└── ...              # TODO: add code / data / notebooks
```

## Guided PPIL4 HADDOCK screen (pipeline step 10)

Files ending in `_(Tyrone)` are the 9DWV-guided exploratory CRBN--candidate--
PPIL4-RRM screen.  The four scripts cover preparation, pre-HADDOCK triage,
sequential execution, and result aggregation.  The three corresponding CSV
files record the 30-candidate shortlist, all campaign statuses, and all final
HADDOCK model scores. `10_guided_haddock_completed_candidates_ranked_(Tyrone).csv`
is ranked by the mean score of each candidate's best ten final models.

These calculations are exploratory structural hypotheses. They are not
experimental confirmation of PPIL4 recruitment, binding affinity, or efficacy.

## Getting started

```bash
git clone https://github.com/KeyuanWang0324/CRBN_MGD.git
cd CRBN_MGD
# TODO: environment setup (e.g. conda/uv, dependencies)
```

## Usage

_TODO: how to run the main analysis / pipeline._

## Data

_TODO: where the data lives, formats, and any access notes._

## Contact

Maintainer: [@KeyuanWang0324](https://github.com/KeyuanWang0324)
