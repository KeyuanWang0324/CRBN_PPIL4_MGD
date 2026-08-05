# CRBN_MGD — CRBN–PPIL4 Molecular Glue Degrader discovery

A computational pipeline that designs and ranks candidate **CRBN molecular glue
degraders (MGDs) of PPIL4**, as a structure-based hypothesis for treating
RAS-driven cancer through the spliceosome.

> **What this is / isn't.** This is a *structural-plausibility enrichment*
> pipeline — it prioritizes chemotypes that are geometrically consistent with
> forming a CRBN–ligand–PPIL4 ternary complex. It is **not** a validated
> predictor of degradation or efficacy, and it contains **no experimental
> confirmation** of PPIL4 recruitment, binding, or anti-cancer activity. Read
> the rankings as triage, not truth.

## Scientific rationale (and its honest boundary)

1. **RAS-driven cancers are selectively dependent on efficient pre-mRNA splicing**,
   driven by oncogenic RAS's elevated transcriptional load (Nat. Commun. 2026;
   the "spliceosome addiction" logic of Nature 2015 for MYC).
2. **PPIL4 is a functional component of the spliceosome catalytic core** — shown
   here independently by DepMap co-dependency (it co-varies with the Prp19/NTC,
   DHX38, CDC40 catalytic-activation module) and by the literature.
3. **PPIL4 is a druggable CRBN glue target** — the CRBN–FPFT-2216–PPIL4 ternary
   complex is solved (PDB **9DWV**).
4. **Hypothesis (what this project tests):** a *selective* PPIL4 degrader could
   impair splicing in a way RAS-driven cancers are especially vulnerable to.
   FPFT-2216 is promiscuous (also degrades IKZF1/3, CK1α, PDE6D, DTWD1), so a
   selective glue is needed to test the PPIL4-specific hypothesis.

**Boundary (do not overstate):** PPIL4 is broadly *common-essential* and is
**not** shown to be KRAS-selective in DepMap knockout data (nor are the validated
factors SF3B1/RBM39 — full-KO CRISPR can't resolve this class of selectivity).
Any tumor selectivity would arise from partial-degradation dosing and RAS's
transcriptional-load dependence, and **remains to be tested** with isogenic
RAS-on/off or sub-lethal splicing assays.

## Pipeline

| Stage | Script | Role | Interpreter |
|---|---|---|---|
| 00 | `00_validate_*` | ground-truth interface / SMILES validation (5 Å heavy-atom contacts vs 9DWV) | system |
| 01 | `01_generate_thalidomide_analogs` | analog library (glutarimide degron fixed; benzo-ring substituents vary) | system (rdkit) |
| 02 | `02_CRBN_binder_test` | CRBN-glue chemotype RF classifier (pooled ChEMBL targets) | system |
| 03 | `03_crbn_binder_scaffold_model` | CRBN-binder RF (scaffold split) + builds the active-candidate list | system |
| 04 | `04_vina_dock_candidates` | Vina docking into CRBN pocket (+ PPIL4) → poses, contacts, screening scores | system (vina/meeko) |
| 05 | `05_haddock3_ternary_novel_candidate` | lite ligand-free HADDOCK3 ternary — fast pre-filter | `.venv-haddock3` |
| 06 | `06_haddock3_ternary_complete` | complete ligand-free HADDOCK3 — coarse pre-filter / candidate ledger | `.venv-haddock3` |
| **08** | `08_haddock3_ternary_with_ligand` | **ligand-inclusive 3-body HADDOCK3 — the final ranking stage (LOCKED, `ppil4_lock_v1`)** | `.venv-haddock3` |
| 07 | `07_extract_top_structures` | extract + chain-color the top ternary models (PyMOL `.pse`) | system + PyMOL |
| 09 | `09_score_vs_fpft2216_reference` | positive control: FPFT-2216 through the pipeline, scored vs 9DWV | `.venv-haddock3` |
| 10 | `10_*_(Tyrone)` | independent 9DWV-*guided* HADDOCK screen — a parallel cross-check pipeline | — |

**Why 08 is the ranking stage, not 05/06:** 05/06 never put the ligand in the CNS
topology, so candidates sharing the glutarimide degron collapse to identical
scores (structurally insensitive to chemistry). 08 docks the ligand as a real
third body, so its chemistry differentiates candidates.

## The locked docking protocol — `PROTOCOL_LOCK.md`

08 and Tyrone's screen (10) originally used different protocols, so their scores
weren't comparable. [`PROTOCOL_LOCK.md`](PROTOCOL_LOCK.md) pins **one** protocol
(`ppil4_lock_v1`): fixed restraints (`ambig_FIXED.tbl`, identical for every
molecule — removes the per-candidate Vina noise), a pinned `iniseed`, fixed
inputs/versions, and a single output schema. An **FPFT-2216 golden reference**
(`REFERENCE_FPFT.json`, §5) must reproduce before candidate runs are trusted.

> **Status:** the lock is wired into 08 and committed. Scores produced under the
> *old* per-candidate-restraint protocol (the current `08_final_*` finalists and
> `final_buy_list`) are **provisional** until re-run under the lock.

## Key outputs

- `final_buy_list_(Ryan).csv` — the deliverable: BUY / ADD / CONSIDER candidates + positive & negative controls, each with a one-line rationale.
- `08_final_ternary_with_ligand_results_(Ryan).csv` — 08's ranked finalists.
- `ternary_structures/` — extracted, chain-colored ternary models (`.pse`).
- `PROTOCOL_LOCK.md`, `ambig_FIXED.tbl`, `crbn_active_FIXED.txt`, `ppil4_active_FIXED.txt` — the locked protocol bundle.

## Environment

Two interpreters (kept separate because HADDOCK3's venv lacks the cheminformatics stack):

```bash
# system Python 3.12 — rdkit, scikit-learn, vina, meeko, requests, biopython, PyMOL
#   used by 00–04, 07, 09
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 "01_generate_thalidomide_analogs_(Ryan).py"

# .venv-haddock3 — the haddock3 / pdb-tools CLIs
#   used by 05, 06, 08 (scripts auto-relaunch into it if not active)
source .venv-haddock3/bin/activate
python3 "08_haddock3_ternary_with_ligand_(Ryan).py"
```

Most scripts detect the wrong interpreter and re-exec into the correct one automatically.

## Limitations (read before trusting any number)

- **No experimental validation.** PPIL4 has no public CRBN-glue activity data; every score is a structural proxy.
- **HADDOCK score is size-biased and interface-dominated** — known glues rank mid-pack; the positive/negative controls in the buy list are the intended calibration.
- **The positive control is partly circular** (PPIL4 restraints derive from 9DWV); a neutral-restraint control is planned to test whether the protocol finds the interface unaided.
- Candidate rankings are **enrichment/triage**, not potency predictions.

## Contact

- Maintainer: [@KeyuanWang0324](https://github.com/KeyuanWang0324) (Ryan)
- Guided PPIL4 screen (`10_*`): Tyrone
