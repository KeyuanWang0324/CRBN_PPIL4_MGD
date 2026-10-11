# Concrete fit of 2026 molecular-glue and deep-learning QSP models to the PPIL4–CRBN project

## What exists now, and what can actually be estimated?

### Takeaway
The active repository already contains a complete, locked structural-triage dataset and a useful but non-identifying Z6466608628 pilot. It does not yet contain the independently replicated concentration-by-time PPIL4-abundance data, absolute protein abundances, ternary binding measurements, or intracellular exposure data required to identify a full molecular-glue kinetic model. A reduced empirical dose-response model is therefore the only defensible fit now; the 2026 molecular-glue ODE model becomes useful after the planned HiBiT/Western experiments, whereas the 2026 deep-learning–QSP model is a PROTAC precedent rather than a directly transferable PPIL4 predictor.

### Cited Findings
- The repository states that all 27 molecules (18 designed candidates, three cross-check compounds, and six controls) completed the ligand-inclusive `ppil4_lock_v1` stage, with FPFT-2216 reproducibility required before accepting runs; it explicitly defines the output as structural-plausibility triage rather than degradation or efficacy prediction. — [Active repository README](../../README.md)
- The stage-08 control table contains FPFT-2216, Z6466608628, TM-04-064-02, and three negative controls. Mean best-10 HADDOCK scores are −111.382, −125.610, −109.961, and −105.097 for FPFT-2216, Z6466608628, TM-04-064-02, and `negctrl_no_ppil4_arm`, respectively; these are ranking proxies, not measured affinities. — [Locked control results](../../08_controls_results_\(Ryan\).csv); [protocol boundary](../../README.md)
- The final locked buy-list table contains compound identity, chemistry, HADDOCK score and dispersion, BSA, cluster fraction, AIR and van der Waals terms, Vina scores, degron-contact overlap, RF probability, and rule decisions. Four candidates pass all three project rules; the current three-compound recommendation is RW_1987, RW_465, and RW_1992. — [Locked top-20 table](../../final_buy_list_lock_v1_top20_\(Ryan\).csv); [repository interpretation](../../README.md)
- The available cellular pilot is one plate with four technical wells per condition and four sequential rereads, not four biological replicates. It tested Z64, TM, and the negative control at only 5 and 10 µM in 293FT-GFP and HeLa-GFP, and treatment was confounded with plate row. — [Fluorescence analysis](../../reports/2026-10-07_GFP_Fluorescence_Assay_Analysis.md); [statistical audit](../../reports/2026-10-08_Z64_Fluorescence_Western_Blot_Statistics_and_Correlation.md)
- The same pilot's Western blot has one lane per condition. Approximate values suggest a dose-associated endogenous-PPIL4 decline in SY5Y, but several bands are saturated, loading varies, and anti-GFP and anti-PPIL4 measurements disagree in 293FT at 10 µM; no Western p-value or reliable DC50/Dmax is supported. — [Western interpretation](../../experimental_results/2026-10-06_Z64_positive_control_western_blot/2026-10-07_Z6466608628_Positive_Control_Western_Blot_Interpretation.md); [densitometry audit](../../reports/2026-10-08_Z64_Fluorescence_Western_Blot_Statistics_and_Correlation.md)
- Baek et al. measured a biochemical PPIL4-recruitment EC50 of 0.34 µM for Z6466608628, observed selective PPIL4 recruitment versus the broader FPFT-2216 profile, and reported only modest/incomplete cellular PPIL4 loss. Recruitment EC50 is not cellular degradation DC50. — [Baek et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12289921/)
- The 2026 biophysical study models free glue, effector/E3, POI, binary complexes, and ternary complex with mass-action ODEs, plus baseline POI synthesis/degradation and a lumped ternary-complex-to-degradation rate. It classifies glue–effector affinity and cooperativity as glue-intrinsic, while protein abundance, baseline POI turnover, baseline effector–POI binding, and effective degradation from ternary complex are system/extrinsic quantities. — [Chae et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13126642/)
- That study fitted ensembles rather than relying on a single best parameter vector, because published relative-abundance assays do not precisely identify all absolute kinetic parameters and some parameters are correlated. Its public implementation is MATLAB and includes fitting to BRET, affinity-selection/iMRM, and HiBiT datasets. — [Chae et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13126642/); [official model repository](https://github.com/jcd111/MG_Manuscript_Model)
- The 2026 deep-learning–QSP study is PROTAC-specific: it uses decomposed full-PROTAC, warhead, linker, and E3-ligand descriptors; predicted POI and E3 binary affinities; POI turnover; and learned corrections to Hook-model parameters. Its dataset has 444 measurements for 272 PROTACs, 53 POIs, and two E3 ligases and was randomly split by measurement. — [Goo et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888238/); [official code](https://github.com/swgoo/protac_deep_qsp)
- The deep-learning–QSP paper reports substantially better generalization for DC50 than Dmax and attributes poor Dmax generalization to missing cell-line, E3-expression, proteasome, duration, and detection-method context. Its official reproduction setup uses an NVIDIA GPU, Docker, and a CUDA PyTorch image. — [Goo et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888238/); [official code requirements](https://github.com/swgoo/protac_deep_qsp)

### Inferences

#### Parameter status for this project

| Quantity | Status now | What would make it estimable | Defensible use |
|---|---|---|---|
| Structural rank / pose robustness | Available for all locked compounds | No new experiment; retain protocol lock | Prioritization and stratified compound selection only |
| Molecular descriptors | Available or directly computable from SMILES | RDKit calculation with pinned version | Covariates; not efficacy labels |
| Apparent DC50 and Dmax | Not estimable from the pilot | 7–8 concentration degradation curves at fixed times, ≥3 independent experiments | Compound × cell line × time-specific empirical endpoints |
| Onset/recovery rates | Not estimable | Multi-timepoint HiBiT/Western, with optional washout | Apparent degradation and recovery kinetics |
| Baseline PPIL4 turnover, γ | Not available | Pulse-chase or carefully controlled synthesis-block experiment, preferably orthogonal to the degradation assay | Fix or tightly constrain rather than co-fit initially |
| Baseline PPIL4 and CRBN abundance, P_T and E_T | Not available in absolute units | Calibrated quantitative Western, targeted proteomics, or calibrated endogenous tags | Cell-line-specific fixed inputs; relative blots alone are insufficient |
| Glue–CRBN K_D | Not measured for candidates | Binary SPR/BLI/ITC or a validated competition assay | Fit/fix separately for each compound |
| PPIL4 binding to glue–CRBN and cooperativity | Not measured | Two-dimensional ternary recruitment matrix or orthogonal binary/ternary biophysics (for example TR-FRET plus SPR/BLI) | Needed to interpret potency mechanistically |
| Effective ternary-complex degradation rate, k_cat-like lumped term | Not separately identifiable from current data | Ternary occupancy/recruitment plus dense degradation time course; ubiquitination time course improves interpretation | Ensemble estimate, not a unique molecular constant |
| Intracellular free compound/exposure | Not measured | Intracellular concentration, permeability/stability, or at minimum washout persistence | Needed to distinguish slow clearance from durable catalytic action |
| Competitive neosubstrate parameters | Not identifiable | Quantitative proteomics plus expression/binding information for competing CRBN substrates | Defer unless selectivity data show competition matters |
| RAS-dependence parameters | Not present in either degrader model | Matched/isogenic RAS-on/off experiment measuring PPIL4 loss, splicing phenotype, and cell outcome | A downstream pharmacodynamic layer, not a docking-derived parameter |

The current docking metrics can be used as preregistered covariates or priors for ranking, but they cannot be converted into K_D, cooperativity, k_cat, DC50, or Dmax. With only a handful of tested PPIL4 compounds, fitting any learned mapping from docking features to degradation will be underdetermined; the scientifically useful test is whether the frozen docking rank enriches for active compounds, not whether a neural network can interpolate three or six molecules.

The Z64 pilot should be treated as assay-development evidence. Its two doses and one biological run cannot distinguish a nonmonotonic hook-like response from row effects, toxicity, fluorescence interference, protein turnover, or cell-line differences. The approximate SY5Y Western trend can motivate the concentration range but should not initialize a precise kinetic parameter.

### Gaps
- The repository does not document treatment duration for the 2026-10-06 Western blot, original unsaturated imager files, antibody identities, or independent replicate blots. — [Western file audit](../../experimental_results/2026-10-06_Z64_positive_control_western_blot/2026-10-07_Z6466608628_Positive_Control_Western_Blot_Interpretation.md)
- The repository does not document whether the GFP construct is a PPIL4 fusion in enough detail to validate bulk GFP as a PPIL4-abundance surrogate. — [Fluorescence analysis](../../reports/2026-10-07_GFP_Fluorescence_Assay_Analysis.md)
- No matched PPIL4 HiBiT dataset, absolute CRBN/PPIL4 abundance, candidate ternary recruitment assay, ubiquitination kinetics, intracellular exposure, or RAS-isogenic response dataset was located in the active repository.
- The 2026 biophysical framework demonstrates qualitative/mechanistic utility under parameter ensembles, but it does not make every fitted kinetic constant uniquely identifiable from relative HiBiT curves alone. — [Chae et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13126642/)

## What minimum experiment would make a useful model?

### Takeaway
Use a staged design. First establish reproducible, direct PPIL4 degradation and empirical DC50/Dmax. Then collect a compact concentration-by-time dataset for active compounds and fit a reduced turnover ODE with ensemble uncertainty. Only after that add ternary biophysics, selectivity, washout, and RAS-dependent phenotypes. Trying to estimate the full molecular-glue network in the first experiment would spend scarce wells on parameters the readout cannot identify.

### Cited Findings
- The molecular-glue model fits relative POI remaining from HiBiT/iMRM and shows that ternary affinity influences Dmax while affinity/cooperativity and glue–effector binding influence potency; nevertheless, system-specific protein levels and degradation capacity bound performance. — [Chae et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13126642/)
- A two-dimensional biochemical matrix can in principle recover molecular-glue thermodynamic parameters more directly than a one-dimensional cellular dose curve. — [Couturier et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12464983/)
- Baek et al. used FPFT-2216 as a broad recruiter reference and Z6466608628 as a more PPIL4-selective recruiter; the latter showed a 0.34 µM biochemical recruitment EC50 but incomplete cellular loss even at high concentration, supporting a range that spans below and above 0.34 µM without assuming it equals cellular DC50. — [Baek et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12289921/)
- The repository's present pilot shows that bulk GFP can be altered by cell number or optical interference and recommends direct endogenous PPIL4 measurement, randomized/block-balanced plates, independent biological repeats, viability/cell count, compound-only wells, and parental GFP-negative controls. — [Fluorescence analysis](../../reports/2026-10-07_GFP_Fluorescence_Assay_Analysis.md)
- The RAS rationale is a distinct downstream biological hypothesis: the 2026 study reports spliceosome induction as a vulnerability in RAS-driven senescence and cancer, but it does not establish PPIL4 degradation as RAS-selective. — [Wagner et al., 2026](https://www.nature.com/articles/s41467-026-71564-z)

### Inferences

#### Stage A — direct degradation screen and assay qualification

- Compounds: Z6466608628 (selective literature anchor), FPFT-2216 (broader positive reference), `negctrl_no_ppil4_arm` (matched structural negative), and the three current candidates RW_1987, RW_465, and RW_1992. TM-04-064-02 may remain a comparator, but should not displace the two literature-anchored controls.
- Readout: endogenous PPIL4-HiBiT if a validated endogenous knock-in is available; otherwise quantitative, unsaturated Western with total-protein normalization and prespecified PPIL4 band boundaries. Measure viability/cell count in matched wells. GFP should be secondary until validated against direct PPIL4.
- Concentrations: vehicle plus 0.01, 0.03, 0.1, 0.3, 1, 3, and 10 µM. This brackets the 0.34 µM biochemical Z64 recruitment EC50 and retains the literature-relevant 10 µM ceiling. If solubility or cytotoxicity fails QC, censor the affected concentrations instead of forcing a Hill fit.
- Times: 6 h and 24 h for the first screen. These do not identify a kinetic model but reveal early versus delayed effects and provide compound-specific ranges for Stage B.
- Replication: at least three independent cell-culture/treatment experiments on different days; technical duplicates or triplicates within each experiment. The biological experiment is the inferential unit.
- Layout/QC: randomized or blocked placement across rows and columns; distributed vehicle wells; constant DMSO; compound-only/no-cell interference wells where fluorescence/luminescence is used; prespecified saturation, bubble, and viability exclusions.
- Fit separately for each compound, cell line, and time using a four-parameter Hill/Emax curve. Report DC50 only when the curve spans both shoulders and Dmax only when a plateau is supported; otherwise report censored bounds (for example, DC50 >10 µM) and observed maximum loss.

#### Stage B — compact kinetic dataset for active compounds

- Advance Z64 plus no more than two active novel candidates and the no-PPIL4-arm negative.
- At each compound, choose vehicle plus three informative concentrations from Stage A: approximately DC20, DC50, and the highest non-cytotoxic plateau concentration.
- Measure PPIL4 at 0, 1, 2, 4, 8, and 24 h; add 48 h only if 24 h has not approached plateau or recovery. Three independent biological experiments remain the minimum.
- Fit a reduced model first: baseline PPIL4 synthesis and turnover plus a concentration-dependent effective induced-loss term. Estimate apparent onset and recovery with bootstrap or ensemble/profile intervals. Do not fit separate K_D1, K_D2, cooperativity, intracellular free drug, and k_cat from this dataset alone.
- Optional high-value washout: expose for 6 h at approximately DC50 and plateau concentration, wash, then measure immediately and 2, 6, 18, and 42 h after washout. This separates durable intracellular/ternary engagement from rapid rebound, but it still does not uniquely identify intracellular exposure without drug measurements.

#### Stage C — mechanism and selectivity

- For an active candidate, repeat a compact dose/time condition with CRBN loss or knockdown/rescue and with proteasome or cullin-neddylation inhibition, including inhibitor-only toxicity controls. This tests whether fitted induced loss is actually CRBN–CRL4/proteasome dependent.
- Confirm direct recruitment with a PPIL4–CRBN ternary assay. A 2D matrix varying both protein components and glue concentration is the appropriate experiment if the aim is to estimate affinities/cooperativity rather than merely rank compounds.
- Quantify PPIL4 and known FPFT-2216 off-target neosubstrates in the same samples, then use unbiased proteomics for only the strongest candidate. Selectivity is an outcome vector, not a single correction to HADDOCK score.

#### Stage D — RAS-dependence as a separate causal layer

- Use a matched/isogenic RAS-on versus RAS-off pair, or a carefully matched mutant/wild-type panel if an isogenic system is unavailable.
- At sublethal concentrations that produce matched PPIL4 degradation, measure PPIL4 abundance, one prespecified splicing phenotype, and viability/proliferation in the same biological experiments. A degradation-resistant PPIL4 rescue would be the strongest causal control if technically feasible.
- Fit a mixed model with experiment/day as a random effect and fixed effects for RAS state, compound, concentration, time, and the RAS-by-treatment interaction. Link outcome to measured PPIL4 remaining rather than nominal dose alone.
- Do not call a difference in viability “RAS dependence” if the RAS-on cells also have greater PPIL4 degradation or different CRBN expression; those variables must be measured and included as mediators/covariates.

### Gaps
- The exact cell line and tagging strategy for the planned HiBiT assay are not documented, so the design must be finalized around available validated cells.
- Available compound amount, solubility limits, plate capacity, and access to CRBN perturbation/proteomics are unknown; these determine whether all three novel candidates or only a stratified subset can be run.
- A precise Stage-B concentration set cannot be chosen until Stage-A viability and degradation curves are observed.

## Which model should be used, what should be deferred, and how would implementation proceed?

### Takeaway
For ISEF, the strongest route is an interpretable, staged “structure-to-degradation” study: preserve the frozen `ppil4_lock_v1` ranking, estimate empirical DC50/Dmax with proper replication, then fit a reduced molecular-glue turnover model and test whether structural rank predicts prospective activity. Adopt concepts and equations from the 2026 biophysical molecular-glue paper. Do not use the 2026 PROTAC deep-learning–QSP network to issue PPIL4 candidate predictions; its molecular representation, labels, split, and learned parameters do not match monovalent CRBN glues.

### Cited Findings
- The molecular-glue biophysical model is modality-matched, explicitly supports HiBiT-like relative abundance data, and has public simulation/fitting code, but its authors caution that absolute kinetic parameters are not precisely identified by much of the published relative-abundance data. — [Chae et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13126642/); [official code](https://github.com/jcd111/MG_Manuscript_Model)
- The deep-learning–QSP implementation expects PROTAC-specific component descriptors and two predicted binary affinities, learns on a random measurement split, and acknowledges that repeated compound–target observations share identical features. Those choices do not establish transfer to a monovalent molecular-glue chemotype or to unseen PPIL4 candidates. — [Goo et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888238/)
- Harder validation materially reduces apparent performance in degrader ML: a 2026 geometric PROTAC classifier fell from 80.81% accuracy under random split to 65.62% under cluster split and 64.08% under temporal split. — [SE(3)-PROTACs](https://pmc.ncbi.nlm.nih.gov/articles/PMC13158127/)
- The PPIL4 repository itself demonstrates why a local model must preserve endpoint boundaries: HADDOCK score is measurably size-biased and the authors explicitly state that the rankings are enrichment/triage rather than potency predictions. — [Active repository README](../../README.md)

### Inferences

#### Decision matrix

| Option | Modality fit | Data fit today | Interpretability | Main risk | Local-compute burden | Recommendation |
|---|---:|---:|---:|---|---:|---|
| Empirical Hill/Emax per time | High for measured endpoint | Medium after one proper experiment | High | Cannot separate mechanism | Low | Start here |
| Reduced PPIL4 turnover ODE | High | Low now; high after Stage B | High | Correlated effective parameters | Low–medium | Core ISEF mechanistic model |
| Full 2026 molecular-glue biophysical ODE | Very high | Low until ternary biophysics and abundance data exist | Very high | Practical non-identifiability | Medium; MATLAB reference can be ported/validated locally | Use as framework; fit in stages/ensembles |
| 2026 PROTAC deep-learning–QSP as published | Low–medium | Very low | Medium | PROTAC-to-glue domain shift; random-split optimism; weak Dmax transfer | High; official setup expects NVIDIA/CUDA | Do not use for PPIL4 predictions |
| Retrained CRBN molecular-glue ML | Potentially high | Very low locally | Variable | Too few standardized PPIL4 labels; leakage across close analogues | Medium–high | Defer; at most a separate public-data benchmark |
| Docking score mapped directly to DC50/Dmax | Low | None | Superficially high | Converts structural triage into unsupported potency | Low | Reject |

#### Phased implementation plan

1. **Freeze and register existing evidence.** Hash/export the current `ppil4_lock_v1` candidate/control table and prespecify the primary comparison: enrichment of direct PPIL4 degraders among the locked top/rule-passing candidates versus controls. No re-ranking after wet-lab labels are seen.
2. **Create one tidy experimental schema.** One row per biological experiment × cell line × compound × concentration × time × replicate, with raw and normalized PPIL4, viability/cell count, plate position, DMSO, assay method, batch, and QC flags.
3. **Fit empirical endpoints.** Estimate compound/time/cell-specific DC50 and Dmax with hierarchical or per-experiment fits; retain censoring and failed-QC states. Compare candidate activity with the frozen structural rank using rank correlation and top-k enrichment only after direct degradation is established.
4. **Fit a reduced ODE.** Fix or independently constrain baseline turnover; estimate an effective concentration-response function and induced degradation rate with bootstrap/profile/ensemble uncertainty. Validate synthetic-parameter recovery before fitting biological data.
5. **Add biophysical parameters only when measured.** Insert measured glue–CRBN affinity, ternary recruitment/cooperativity, and absolute E/P abundance sequentially; compare predictive checks and parameter profiles after each addition. Retain parameter ensembles when multiple solutions fit similarly.
6. **Attach selectivity and RAS outcomes as downstream modules.** Do not fold them into a universal score. Model selectivity as multiple protein outcomes and RAS dependence as an interaction conditioned on measured PPIL4 loss.
7. **Defer deep learning.** A later public-data project could reproduce the PROTAC QSP paper locally, then test scaffold/target/time splits and audit failure. It should remain method development, not evidence that the PPIL4 candidates will degrade PPIL4.

#### Candidate files/scripts to add later (no files added in this task)

- `modeling/configs/ppil4_model_v1.yaml` — units, concentration/time grids, fixed versus fitted parameters, seeds, and assay-specific error models.
- `modeling/01_build_experimental_long_table.py` — validate IDs, units, biological-replicate structure, QC fields, and censoring.
- `modeling/02_fit_hill_dc50_dmax.py` — per-time empirical fits plus confidence intervals and diagnostic plots.
- `modeling/03_fit_ppil4_turnover_ode.py` — reduced ODE using local SciPy integration/optimization or Bayesian sampling.
- `modeling/04_identifiability_and_profiles.py` — multi-start fits, profile likelihoods, parameter correlations, and synthetic recovery.
- `modeling/05_compare_structure_to_degradation.py` — preregistered rank/enrichment analysis; never train and evaluate on the same tiny candidate panel.
- `modeling/06_ras_dependence_mixed_model.py` — RAS-state interaction analysis linked to measured PPIL4 remaining.
- `modeling/tests/test_synthetic_recovery.py` — assert recovery of identifiable synthetic parameters and demonstrate non-identifiability of deliberately confounded ones.
- `experimental_results/model_input/ppil4_degradation_long.csv` — immutable tidy input derived from raw assay files, not hand-entered summary means.
- `reports/PPIL4_model_card.md` — endpoints, training/fitting data, exclusions, parameter identifiability, validation boundary, and prohibited claims.

All docking, simulation, fitting, resampling, and any future deep-learning reproduction should run on the user's local machine under pinned environments. The mechanistic ODE and bootstrap work is modest compared with HADDOCK; the published deep-QSP reproduction is the computationally expensive branch and should be deferred unless a local NVIDIA/CUDA setup is available.

### Gaps
- The molecular-glue paper's official reference implementation is MATLAB; the user's available MATLAB license/toolboxes are unknown. A Python port would need numerical validation against published/reference trajectories before biological use. — [Official model repository](https://github.com/jcd111/MG_Manuscript_Model)
- The published deep-QSP model has not been validated on molecular glues, PPIL4, or this thalidomide-analogue series. No defensible conversion from the repository's HADDOCK score to the model's required K_D inputs was found.
- With only three recommended candidates plus controls, this project can test prospective enrichment and mechanistic consistency, but cannot train a reliable compound-level deep network.
