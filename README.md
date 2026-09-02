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

## Project website

`CRBN_Project_site/` is the source of the public write-up, deployed to a **separate repo**
([KeyuanWang0324/CRBN_Project](https://github.com/KeyuanWang0324/CRBN_Project)) and served at
<https://keyuanwang0324.github.io/CRBN_Project/>. Deployment is a file copy into that repo — the
two are not linked, so they can drift if only one is updated.

| File | What it is |
|---|---|
| `index.html` | the page. ~89% of it is generated: the 27 candidate/control records sit between marker comments and are rewritten by 14. Only the prose, CSS and score chart outside those markers are hand-edited. |
| `glb/` | 27 binary-glTF meshes exported from PyMOL by 16 (~27 MB) |
| `structures.json` | index of mesh URLs, each with that structure's ligand centre |
| `viewer.js` | three.js viewer; builds on record open, disposes on close (WebGL contexts are capped ~16, the page has 27) |
| `glossary.js` | 44 hover definitions, each with an analogy; annotates text nodes at load |
| `lib/` | vendored three.js + GLTFLoader + OrbitControls — the page makes **no** external request |
| `img/` | static PyMOL renders, used as the no-JavaScript fallback |

Rebuild order after any pipeline change: **11 → 12 → 13 → 16 → 14**. 14 must run last: it injects
the records and stamps a content hash onto `viewer.js`, `glossary.js`, `structures.json` and `lib/`
so browsers re-fetch exactly what changed. GitHub Pages takes roughly 30–60 s to serve a push.

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

## Ordering, and how the vendor sheet was checked

Only **3 of the 25** compounds on the order sheet are registered substances (TZ_6 244057-33-8,
Z6466608628 2869183-63-9, negctrl_no_ppil4_arm 26581-81-7), confirmed by exact-structure PubChem
lookup. The other 22 are novel and can only be made to order, so **the structure is the
specification** — which is why `Final Buy List for Purchase.sdf` exists and is the file to send.

The sheet contains **five pairs of constitutional isomers** — `RW_491`/`RW_2392`,
`RW_2123`/`RW_1668`, `RW_381`/`RW_1092`, `RW_1078`/`RW_80`, and `RW_72`/`Z6466608628` — so the
molecular formula is *not* a specification on its own. All 25 InChIKeys are unique; quote those.

The vendor returned `副本生物物料采购(3).xlsx` with a drawn structure for every row. All 25 were
checked against our own depictions, scaffold-aligned so substituent positions compare directly:
**25/25 correct, including all five isomer pairs.** Their drawing tool suppresses hydrogens, so
`RW_1959`'s hydroxyl and `RW_2392`'s carboxylic acid render as a bare `O` — verified against
`RW_491`, where a real methyl is drawn.

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

## References

Harvard style, matching the project page. Sources whose full author list could not be verified are
cited by organisation in the online-resource form rather than guessed.

**Target, prior art and rationale**

- Baek, K., Metivier, R.J., Roy Burman, S.S., Bushman, J.W., Yoon, H., Lumpkin, R.J., Abeja, D.M.,
  Lakshminarayan, M., Yue, H., Ojeda, S., Verano, A.L., Gray, N.S., Donovan, K.A. and Fischer, E.S.
  (2025) 'Unveiling the hidden interactome of CRBN molecular glues', *Nature Communications*, 16,
  6831. doi:[10.1038/s41467-025-62099-w](https://doi.org/10.1038/s41467-025-62099-w).
  — Z6466608628 and the CRBN–FPFT-2216–PPIL4 structure. **The same paper covers both**; they were
  previously cited separately.
- Hsu, T.Y. et al. (2015) 'The spliceosome is a therapeutic vulnerability in MYC-driven cancer',
  *Nature*, 525(7569), pp. 384–388. doi:[10.1038/nature14985](https://doi.org/10.1038/nature14985).
- Petzold, G., Gainza, P. et al. (2025) 'Mining the CRBN target space redefines rules for molecular
  glue–induced neosubstrate recognition', *Science*, 389(6755).
  doi:[10.1126/science.adt6736](https://doi.org/10.1126/science.adt6736). — where PPIL4 came from.
- RCSB Protein Data Bank (2025) *Entry 9DWV*. Available at: <https://www.rcsb.org/structure/9DWV>
  (Accessed: 2 September 2026).

**Methods and software**

- Basu, S. and Wallner, B. (2016) 'DockQ: a quality measure for protein–protein docking models',
  *PLOS ONE*, 11(8), e0161879. doi:[10.1371/journal.pone.0161879](https://doi.org/10.1371/journal.pone.0161879).
- Dominguez, C., Boelens, R. and Bonvin, A.M.J.J. (2003) 'HADDOCK: a protein–protein docking
  approach based on biochemical or biophysical information', *Journal of the American Chemical
  Society*, 125(7), pp. 1731–1737. doi:[10.1021/ja026939x](https://doi.org/10.1021/ja026939x).
- Eberhardt, J., Santos-Martins, D., Tillack, A.F. and Forli, S. (2021) 'AutoDock Vina 1.2.0: new
  docking methods, expanded force field, and Python bindings', *Journal of Chemical Information and
  Modeling*, 61(8), pp. 3891–3898. doi:[10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203).
- Jumper, J. et al. (2021) 'Highly accurate protein structure prediction with AlphaFold', *Nature*,
  596(7873), pp. 583–589. doi:[10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2).
- Pedregosa, F. et al. (2011) 'Scikit-learn: machine learning in Python', *Journal of Machine
  Learning Research*, 12, pp. 2825–2830.
- RDKit (2026) *RDKit: open-source cheminformatics*. Available at: <https://www.rdkit.org/>
  (Accessed: 2 September 2026).
- Schrödinger, LLC (2026) *The PyMOL molecular graphics system*, version 3.1.6. Available at:
  <https://pymol.org/> (Accessed: 2 September 2026).
- Schüttelkopf, A.W. and van Aalten, D.M.F. (2004) 'PRODRG: a tool for high-throughput
  crystallography of protein–ligand complexes', *Acta Crystallographica Section D*, 60(8),
  pp. 1355–1363. doi:[10.1107/S0907444904011679](https://doi.org/10.1107/S0907444904011679).
- three.js (2026) *three.js JavaScript 3D library*, r140. Available at: <https://threejs.org/>
  (Accessed: 2 September 2026).
- Trott, O. and Olson, A.J. (2010) 'AutoDock Vina: improving the speed and accuracy of docking with
  a new scoring function, efficient optimization, and multithreading', *Journal of Computational
  Chemistry*, 31(2), pp. 455–461. doi:[10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334).
- van Zundert, G.C.P. et al. (2016) 'The HADDOCK2.2 web server: user-friendly integrative modeling
  of biomolecular complexes', *Journal of Molecular Biology*, 428(4), pp. 720–725.
  doi:[10.1016/j.jmb.2015.09.014](https://doi.org/10.1016/j.jmb.2015.09.014).

**Databases and compound sources**

- Broad Institute (2026) *DepMap portal*. Available at: <https://depmap.org/> (Accessed: 2 September 2026).
- European Bioinformatics Institute (2026) *ChEMBL database*. Available at:
  <https://www.ebi.ac.uk/chembl/> (Accessed: 2 September 2026).
- MedChemExpress (2026) *Z6466608628 (HY-175599)*. Available at: <https://www.medchemexpress.com/>
  (Accessed: 2 September 2026).
- National Center for Biotechnology Information (2026) *PubChem*. Available at:
  <https://pubchem.ncbi.nlm.nih.gov/> (Accessed: 2 September 2026).

> **Not listed:** the *Nature Communications* 2026 paper cited above for the RAS–splicing
> dependency. It could not be verified; the claim is hedged as emerging evidence rather than
> asserted, and an unverified source is worse than none. Supply the citation and it goes in.

## Contact

- Keyuan (Ryan) Wang — design & docking pipeline — [@KeyuanWang0324](https://github.com/KeyuanWang0324)
- Zhentai (Tyrone) Zong — independent guided PPIL4 screen (`10_*`)
