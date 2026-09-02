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

1. **RAS-driven cancers may be selectively dependent on efficient pre-mRNA splicing**,
   driven by oncogenic RAS's elevated transcriptional load (Nat. Commun. 2026). Note the
   asymmetry: "spliceosome addiction" is *established* for MYC (Nature 2015); for RAS this
   is emerging evidence and an argument by analogy.
2. **PPIL4 is a spliceosome-associated peptidyl-prolyl isomerase** (cyclophilin family,
   with an RRM domain) acting in catalytic *activation* of the spliceosome rather than
   forming part of the catalytic core. Our own DepMap co-dependency analysis places it
   alongside the Prp19/NTC complex and the second-step factors DHX38 and CDC40 — our
   analysis, not an established result, and a weak signal for a broadly essential gene.
3. **PPIL4 is a druggable CRBN glue target** — the CRBN–FPFT-2216–PPIL4 ternary
   complex is solved (PDB **9DWV**). PPIL4 was surfaced as a CRBN neosubstrate by
   Petzold, Gainza *et al.*, who mined the CRBN target space computationally and then
   confirmed compound-dependent recruitment of several previously uncharacterised
   neosubstrates, PPIL4 among them (*Science*, 2025).
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
| 01 | `01_generate_thalidomide_analogs` | analog library (glutarimide CRBN-binding warhead fixed; benzo-ring substituents vary) | system (rdkit) |
| 02 | `02_CRBN_binder_test` | CRBN-glue chemotype RF classifier (pooled ChEMBL targets) | system |
| 03 | `03_crbn_binder_scaffold_model` | CRBN-binder RF (scaffold split) + builds the active-candidate list | system |
| 04 | `04_vina_dock_candidates` | Vina docking into CRBN pocket (+ PPIL4) → poses, contacts, screening scores | system (vina/meeko) |
| 05 | `05_haddock3_ternary_novel_candidate` | lite ligand-free HADDOCK3 ternary — fast pre-filter | `.venv-haddock3` |
| 06 | `06_haddock3_ternary_complete` | complete ligand-free HADDOCK3 — coarse pre-filter / candidate ledger | `.venv-haddock3` |
| **08** | `08_haddock3_ternary_with_ligand` | **ligand-inclusive 3-body HADDOCK3 — the final ranking stage (LOCKED, `ppil4_lock_v1`)** | `.venv-haddock3` |
| 07 | `07_extract_top_structures` | extract + chain-color the top ternary models (PyMOL `.pse`) | system + PyMOL |
| 09 | `09_score_vs_fpft2216_reference` | positive control: FPFT-2216 through the pipeline, scored vs 9DWV | `.venv-haddock3` |
| 10 | `10_*_(Tyrone)` | independent 9DWV-*guided* HADDOCK screen — a parallel cross-check pipeline | — |
| 11 | `11_assemble_buy_list` | assembles the buy list (top-N over everything docked under the lock) | system (rdkit) |
| 12 | `12_annotate_purchase_list` | order sheet + SDF + structure sheet; PubChem registry lookup | system (rdkit) |
| 13 | `13_render_ternary_views` | static PyMOL renders, one shared camera | PyMOL |
| 14 | `14_build_candidate_panels` | per-candidate panels injected into the project page | system (rdkit) |
| 15 | `15_extract_viewer_structures` | trimmed PDBs with PyMOL `dss` secondary structure | PyMOL |
| 16 | `16_export_pymol_gltf` | PyMOL cartoon geometry as binary glTF for the web viewer | PyMOL |

**Why 08 is the ranking stage, not 05/06:** 05/06 never put the ligand in the CNS
topology, so candidates sharing the glutarimide warhead collapse to identical
scores (structurally insensitive to chemistry). 08 docks the ligand as a real
third body, so its chemistry differentiates candidates.

## The locked docking protocol — `PROTOCOL_LOCK.md`

08 and Tyrone's screen (10) originally used different protocols, so their scores
weren't comparable. [`PROTOCOL_LOCK.md`](PROTOCOL_LOCK.md) pins **one** protocol
(`ppil4_lock_v1`): fixed restraints (`ambig_FIXED.tbl`, identical for every
molecule — removes the per-candidate Vina noise), a pinned `iniseed`, fixed
inputs/versions, and a single output schema. An **FPFT-2216 golden reference**
(`REFERENCE_FPFT.json`, §5) must reproduce before candidate runs are trusted.

> **Status: complete.** All 27 molecules — 18 designed candidates, 3 cross-check,
> 6 controls — have been docked under `ppil4_lock_v1`, and every row of
> `08_final_*`, `08_crosscheck_*`, `08_controls_*` and the buy list carries that
> `protocol_version`. The FPFT-2216 golden reference reproduces byte-identically.
> The only pre-lock artefact left is `final_buy_list_lock_v1_(Ryan).csv`, the
> hand-assembled 10-candidate list that `final_buy_list_lock_v1_top20_(Ryan).csv`
> supersedes.

## Cross-checking the two screens on the same molecules

Ryan's funnel and Tyrone's guided screen (10) rank different molecules, so agreement was
never directly testable. Three of Tyrone's picks — Candidates 118, 45 and 6 — are now docked
through Ryan's wrapper under the lock as `tyrone_118` / `tyrone_45` / `tyrone_6`
(see [`PROTOCOL_LOCK.md` §8a](PROTOCOL_LOCK.md)). All three happen to be in 01's library
already under a `cand_*` name; none had ever been docked in 08.

**They do not agree.** All three score below the ~-112 noise floor (`mean_best_10` -98.97 to
-108.56), i.e. worse than every one of Ryan's 18 and worse than a degron-dead negative. That
is a statement about the *score*, not a verdict on the molecules: all three are bare
phthalimides with no PPIL4-facing arm, and the HADDOCK score is size- and
interface-dominated (see Limitations). The disagreement is real rather than a setup bug —
the FPFT-2216 golden reference reproduced byte-identically in the same session.

> **Tyrone's `completed_results_with_controls` CSV is not merged into the buy list, and
> must not be.** It is `ppil4-lock-v3`, not `ppil4_lock_v1`: the §5 acceptance test fails
> between the two (same FPFT-2216, DockQ 0.376 vs 0.042, different cluster), and on the
> three molecules run under *both* protocols v3 scores the same molecule 14.5 units more
> negative on average (sd 6.3 — not a constant offset). Molecules from that screen enter the
> buy list only by being re-docked under v1. See [`PROTOCOL_LOCK.md` §8b](PROTOCOL_LOCK.md).

## Key outputs

- `Final Buy List for Purchase.csv` — **the order sheet**: 编号 / 分子式 / 结构式 (SMILES) / 分子量 / CAS号 / 纯度级别 / 质量要求 for the 20 candidates and 5 controls.
- `Final Buy List for Purchase.sdf` — **send this one**: the same 25 as real 2D structures with the order fields attached. The list contains five pairs of constitutional isomers, so the formula alone is not a specification.
- `Final Buy List Structures.svg` — printable structure sheet.
- `final_buy_list_lock_v1_top20_(Ryan).csv` — **the deliverable**: the 20 best-scoring molecules docked under the lock (Ryan's 18 + 2 cross-check) + the 5 controls, each with a one-line rationale. Generated by 11, so adding a molecule is a re-run, not a hand edit.

> **Two naming schemes, on purpose.** The buy list names molecules by whose screen they came
> from — `RW_<n>` for Ryan's, `TZ_<n>` for Tyrone's (controls keep their descriptive names).
> Everything upstream — run directories, 04/05/08/09's CSVs — still uses the pipeline names
> (`cand_708`, `tyrone_6`). The buy list carries a **`run_name`** column holding the pipeline
> name, so every row traces back to the run that produced it. `RW_708` = `cand_708`;
> `TZ_6` = `tyrone_6`.
- `final_buy_list_lock_v1_(Ryan).csv` — the earlier hand-assembled 10-candidate list it supersedes.
- `final_buy_list_crosscheck_(Ryan).csv` — the same schema for the cross-check molecules alone (11).
- `08_final_ternary_with_ligand_results_(Ryan).csv` — 08's ranked finalists.
- `08_crosscheck_results_(Ryan).csv` — 08's locked-schema results for the cross-check set.
- `ternary_structures/` — extracted, chain-colored ternary models (`.pse`).
- `PROTOCOL_LOCK.md`, `ambig_FIXED.tbl`, `crbn_active_FIXED.txt`, `ppil4_active_FIXED.txt` — the locked protocol bundle.

## Environment

**Three** interpreters, not two — PyMOL is its own:

```bash
# 1. system Python 3.12 — rdkit, scikit-learn, vina, meeko, requests, biopython
#    used by 00–04, 11, 12, 14
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 "01_generate_thalidomide_analogs_(Ryan).py"

# 2. .venv-haddock3 — the haddock3 / pdb-tools CLIs and haddock.libs
#    used by 05, 06, 08, 09  (scripts auto-relaunch into it if not active)
source .venv-haddock3/bin/activate
python3 "08_haddock3_ternary_with_ligand_(Ryan).py"

# 3. PyMOL's own bundled Python — used by 07, 13, 15, 16
/Applications/PyMOL.app/Contents/MacOS/PyMOL -cq "16_export_pymol_gltf_(Ryan).py"
```

Most scripts detect the wrong interpreter and re-exec into the correct one automatically.

> **PyMOL is not importable from the system Python here, and `pip install
> pymol-open-source` does not fix it** — the macOS wheel ships rpaths pointing at
> the build machine's environment, so it fails on a missing `libpng16`. The
> installed **PyMOL.app** works headlessly and is what 07/13/15/16 use. Note also
> that PyMOL's `run` rebinds `__file__` to its own package directory, which is why
> those scripts resolve the project directory from `cwd`/`$ISEF_DIR` instead.

## What the ranking says

Four conclusions the numbers support, each with its limit. Correlations are Spearman.

1. **The ranking axis is substantially a size axis** — ρ = −0.58 against heavy-atom count across
   the 26 plotted molecules, still −0.39 inside the candidates' 30–34 atom band. Not purely size:
   TM-04-064-02 is the largest molecule here (43 atoms) and scores near the bottom.
2. **Z6466608628 is the comparison worth making** — at 30 heavy atoms it is the only size-matched
   reference, so beating it is not automatically a size artefact the way beating FPFT-2216 (20)
   would be. 7 of 20 clear it.
3. **The three rules are near-independent**, so their intersection is small by construction: score
   vs structure axis is ρ = −0.08; individually 6, 9 and 7 of 18 pass A, B and C, but only 4 pass
   all three. That scarcity is a property of asking three unrelated questions, not evidence the
   four are exceptional.
4. **Score and funnel quality disagree at the top** — the two best-scoring molecules both fail the
   structure rule, and cluster convergence runs slightly *against* score (ρ = +0.20).

Suggestive but **not** established: all four molecules passing everything share a mono-oxo
isoindolinone core and none of the six di-oxo phthalimides passes — but that is 4/12 vs 0/6,
Fisher exact p = 0.25. The di-oxo set actually scores *better* on average (−126.1 vs −124.5); what
separates them is the structure axis, not the score.

## Limitations (read before trusting any number)

- **No experimental validation of these designs.** None of the ranked candidates has been made or
  tested; every score is a structural proxy. Note this is *not* the same as saying the target is
  unvalidated — PPIL4 is a demonstrated CRBN-glue target (FPFT-2216, Z6466608628; ternary complex
  solved as PDB 9DWV). A PPIL4-**selective** glue has been reported too — Z6466608628 (a control here), which recruits
  PPIL4 at EC50 0.34 µM but degrades it only modestly and is presented as a lead for optimisation.
  The gap is a *potent* selective degrader, not the absence of any.
- **HADDOCK score is size-biased and interface-dominated.** Measured on this set: score tracks
  heavy-atom count at Spearman ρ = −0.58 across the 26 plotted molecules, and still −0.39 inside
  the candidates' narrow 30–34 atom band. The control set does **not** remove that bias — it only
  makes scores comparable across runs. One piece of luck: the bar itself, Z6466608628 at 30 heavy
  atoms, *is* size-matched to the candidates, so clearing it is not purely a size effect. The noise
  floor is not (degron-dead negative 23 atoms, FPFT-2216 20), so part of why those sit low is simply
  that they are small.
- **The positive control is partly circular** (PPIL4 restraints derive from 9DWV); a neutral-restraint control is planned to test whether the protocol finds the interface unaided.
- Candidate rankings are **enrichment/triage**, not potency predictions.

## Contact

- Keyuan (Ryan) Wang — design & docking pipeline — [@KeyuanWang0324](https://github.com/KeyuanWang0324)
- Zhentai (Tyrone) Zong — independent guided PPIL4 screen (`10_*`)
