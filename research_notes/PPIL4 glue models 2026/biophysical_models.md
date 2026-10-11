# 2026 molecular-glue biophysical model for CRBN–PPIL4

## What does “2026 molecular-glue biophysical model” most likely refer to?

### Takeaway
The phrase is not merely a generic description: it closely matches Chae et al., *Biophysical Modeling Elucidates Mechanistic Principles for Rational Molecular Glue Design*, a peer-reviewed 2026 JCIM paper. This is the most likely intended model because it explicitly includes a CRBN molecular-glue degrader module, cellular degradation, competition among neosubstrates, and open author code.

### Cited Findings
- Chae et al. published the paper online on 13 April 2026 and in JCIM 66(8), 4888–4902; DOI `10.1021/acs.jcim.5c03125`. The model is a unified mass-action framework for glue stabilizers and degraders, validated against multiple published biochemical and cellular data sets. — [Primary paper](https://doi.org/10.1021/acs.jcim.5c03125); [PubMed record](https://pubmed.ncbi.nlm.nih.gov/41973479/)
- The degrader module explicitly represents CRBN as the effector, a glue, a protein of interest, binary complexes, ternary-complex formation, basal target synthesis/degradation, and glue-induced target degradation. It was fitted to CRBN-modulator data from Sperling et al. and Yamanaka et al. — [Primary paper, Methods](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125); [Sperling et al.](https://doi.org/10.1182/blood.2019000789); [Yamanaka et al.](https://doi.org/10.1038/s41467-023-40385-9)
- The authors released MATLAB models and fitting code, including a “Combined Degrader Model,” at `jcd111/MG_Manuscript_Model`; the repository states that its degrader model is fitted to BRET, AlphaScreen, iMRM, and HiBiT data. — [Author repository](https://github.com/jcd111/MG_Manuscript_Model)
- A directly complementary 2025 paper by Couturier et al. provides a thermodynamic framework for extracting glue binding and cooperativity parameters from a biochemical concentration matrix using FRET/fluorescence-polarization data; Chae et al. then extend that kind of characterization to downstream cellular effects. — [Couturier et al., JACS 2025](https://doi.org/10.1021/jacs.5c09166); [Chae et al., Discussion](https://doi.org/10.1021/acs.jcim.5c03125)

### Inferences
- The user’s phrase almost certainly refers to Chae et al. rather than to a structure-prediction neural network: the title is an exact semantic match, the year is 2026, and the paper contains the requested CRBN/degradation kinetic machinery. — [Primary paper](https://doi.org/10.1021/acs.jcim.5c03125)
- For this PPIL4 project, the Chae model is best treated as a mechanistic pharmacology layer placed **after** structural docking and experimental ternary-complex characterization, not as another docking/ranking score. — [Primary paper](https://doi.org/10.1021/acs.jcim.5c03125); [Couturier et al.](https://doi.org/10.1021/jacs.5c09166)

### Gaps
- No evidence was found that Chae et al. fitted PPIL4-specific parameters. Their CRBN examples concern established IMiD neosubstrates, so CRBN–PPIL4 requires new calibration rather than direct transfer of fitted IKZF/GSPT1 values. — [Primary paper](https://doi.org/10.1021/acs.jcim.5c03125)

## What equations and parameters does the model require?

### Takeaway
It is a coupled mass-action ODE model, not a deep-learning predictor. Its minimal CRBN–PPIL4 version needs four reversible binding steps plus target synthesis, basal turnover, and an effective ternary-complex-to-degradation step.

### Cited Findings
- Let `E = CRBN/CRL4^CRBN effector`, `P = PPIL4`, and `G = glue`. The four reversible reactions are `E + G <-> EG` (`K_D1`), `EG + P <-> EGP` (`K_D2`), `E + P <-> EP` (`K_D3`), and `EP + G <-> EGP` (`K_D4`). Cooperativity satisfies the thermodynamic cycle `alpha = K_D3/K_D2 = K_D1/K_D4`; thus `alpha > 1` is positive cooperativity under this convention. — [Chae et al., CRBN model](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125); [Couturier et al.](https://doi.org/10.1021/jacs.5c09166)
- For each reversible reaction `x`, Chae et al. use `r_forward,x = k_fx[A][B]` and `r_reverse,x = k_fx K_Dx[C]`, so `k_off,x = k_on,x K_Dx`. Their published implementation assumes diffusion-limited association and uses the same association-rate constant for the binding reactions, although the framework can be generalized if measured `k_on/k_off` values exist. — [Chae et al., equations 22–23 and Methods](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125); [Author code](https://github.com/jcd111/MG_Manuscript_Model)
- A compact PPIL4 balance is `d[P]/dt = v_P - gamma[P] - r_f(E+P) + r_r(EP) - r_f(EG+P) + r_r(EGP)`. The ternary species additionally loses material through `k_cat[EGP]`; the model recycles free CRBN and glue after this event. The complete published system contains six ODEs for `P, E, G, EP, EG, EGP`. — [Chae et al., equations 24–29](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- `v_P` is basal PPIL4 production (concentration/time), `gamma` is basal PPIL4 degradation (1/time), and `k_cat` is the **effective combined** ubiquitination-plus-degradation rate from the ternary complex (1/time). Chae et al. choose `v_P` so the pre-dose system is at steady state and set initial glue concentration equal to the administered dose in their assay simulations. — [Chae et al., equations 24–34](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- The model can include a competing neosubstrate `I` through free `I`, competing ternary complex `EGI`, its abundance/production/turnover, its `K_D5`, and its effective degradation rate `k_cat*`; Chae et al. predict strongest interference when a competitor is abundant, binds well, but is degraded slowly or not at all. — [Chae et al., equations 39–44 and Figure 8](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)

### Inferences
- A first PPIL4 implementation can fit equilibrium `K_D1`, `K_D3`, and `alpha` (hence `K_D2` and `K_D4`) and use a lumped `k_cat`; only a later model should separate ubiquitination, polyubiquitin-chain elongation, ternary dissociation, proteasome delivery, and proteolysis. — [Chae et al., limitations](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- If kinetic data are available, replacing the shared-association-rate assumption with reaction-specific `k_on` and `k_off` values would allow ternary residence time to be tested separately from equilibrium cooperativity. — [Chae et al., Methods](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)

### Gaps
- The published Chae model lumps ubiquitination and proteasomal destruction into one `k_cat`; it therefore cannot identify which of those downstream steps limits PPIL4 loss without additional states and measurements. — [Chae et al., limitations](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- The article notes parameter correlations/identifiability limitations because many validation assays report relative rather than absolute abundances; PPIL4-specific fits need uncertainty intervals and profile/ensemble fitting, not single best-fit numbers. — [Chae et al., Discussion and Supporting Information description](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)

## What experimental inputs would a CRBN–PPIL4 model need?

### Takeaway
The project already has the right biological system and positive control, but candidate-specific docking results are not model inputs. The highest-value next data are concentration-matrix ternary recruitment, time-resolved cellular PPIL4 degradation, protein abundances, basal PPIL4 half-life, and ubiquitination kinetics.

### Cited Findings
- The PPIL4 discovery paper measured CRBN–compound–PPIL4 recruitment by TR-FRET, obtained cryo-EM structures of CRBN–FPFT-2216–PPIL4, demonstrated in-vitro PPIL4 ubiquitylation, and measured cellular PPIL4 abundance. Its focused screen reported Z6466608628 with apparent recruitment EC50 0.34 µM versus 1.05 µM for FPFT-2216 in that experiment, but cellular PPIL4 downregulation was selective and modest/incomplete even at 10 µM. — [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12289921/)
- The same PPIL4 study used 20 or 50 nM biotinylated neosubstrate, 2 nM terbium-streptavidin, and 200 nM GFP–CRBN–DDB1ΔB for TR-FRET ternary formation, providing an assay starting point but not, by itself, a full thermodynamic parameterization. — [Primary PPIL4 paper, Methods](https://www.nature.com/articles/s41467-025-62099-w)
- Couturier et al. show that a two-dimensional biochemical matrix varying component concentrations, globally fitted to the complete equilibrium model, can recover glue-binding and cooperativity parameters rather than only an apparent EC50. — [Couturier et al., JACS 2025](https://doi.org/10.1021/jacs.5c09166)
- Chae et al. fit ternary-proximity readouts (BRET/AlphaScreen) and target-abundance readouts (iMRM/HiBiT) jointly; their outputs are the fold increase in `EP + EGP` and the fraction of total target remaining relative to baseline. — [Chae et al., equations 37–38](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- Sperling et al. experimentally showed that CRBN-modulator activity depends on drug-induced ligase–substrate interaction strength, CRBN abundance, and competing-substrate expression; excess nonessential substrate can inhibit degradation of another substrate by competing for a limiting ligase pool. — [Sperling et al., Blood 2019](https://doi.org/10.1182/blood.2019000789)

### Inferences
- Minimum PPIL4 data package for quantitative fitting:
  1. `K_D1` or `k_on/k_off` for compound–CRBN/DDB1 binding (SPR, BLI, ITC, FP/competition assay).
  2. A CRBN × PPIL4 × glue concentration matrix for ternary recruitment, globally fitted to obtain `K_D3`, `K_D2`, and `alpha`; run FPFT-2216 and Z6466608628 as anchors and `negctrl_no_ppil4_arm` as a mechanistic negative control.
  3. Time-resolved PPIL4 degradation over multiple doses (preferably endogenous or calibrated HiBiT; Western blot can validate but is less suitable as the only quantitative fit source) to estimate `D_max`, `DC50/EC50`, lag, and effective `k_cat`.
  4. Basal PPIL4 half-life, giving `gamma = ln(2)/t_1/2`, plus calibrated starting PPIL4 and CRBN/CRL4 abundances in each cell line.
  5. In-vitro CRL4^CRBN ubiquitination time courses to test whether incomplete cellular loss reflects weak recruitment versus inefficient ubiquitin transfer/processing.
  6. Permeability/free intracellular exposure if biochemical and cellular potencies are to be linked quantitatively. — [Chae et al.](https://doi.org/10.1021/acs.jcim.5c03125); [Couturier et al.](https://doi.org/10.1021/jacs.5c09166); [PPIL4 study](https://doi.org/10.1038/s41467-025-62099-w)
- Because spliceosome-bound PPIL4 may occlude the CRBN-recruitment surface, a model extension with free and spliceosome-bound PPIL4 states could test whether an inaccessible PPIL4 pool explains the degradation plateau reported for Z6466608628. — [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w)

### Gaps
- The published PPIL4 paper reports apparent TR-FRET recruitment EC50 values, not all four microscopic dissociation constants or PPIL4-specific `k_on/k_off`; these must be measured or jointly inferred. — [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w)
- Absolute CRBN and PPIL4 abundances, PPIL4 basal half-life, and a PPIL4-specific effective `k_cat` were not identified in the cited sources; these remain project-specific measurements.

## What useful outputs could the model provide for this project?

### Takeaway
After calibration, the model can turn biochemical recruitment and cell-degradation measurements into interpretable bottlenecks and prospective dose/time predictions. Before calibration, it is still useful for assay design and sensitivity analysis, but not for ranking the current docked candidates by degradation potency.

### Cited Findings
- Chae et al. find that `K_D2`—PPIL4 binding to the preformed glue–CRBN complex in this application—drives the maximum ternary/degradation effect, while stronger glue–CRBN binding and higher cooperativity improve potency; strong binding alone does not overcome limits imposed by effective degradation rate or CRBN abundance. — [Chae et al., Results/Discussion](https://doi.org/10.1021/acs.jcim.5c03125)
- Their simulations predict that degrader performance is ultimately limited by the effective ternary-complex degradation rate relative to the target’s natural degradation rate, and by effector abundance; degraders tolerate a broader range of effector:target ratios than noncatalytic stabilizers but still fail when CRBN becomes limiting. — [Chae et al., Figures 6–7](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- The model can generate concentration- and time-dependent ternary-complex abundance, fraction of PPIL4 remaining, `D_max`, an apparent half-maximal degradation concentration, and local sensitivity to binding, abundance, and turnover parameters. — [Chae et al., equations 37–38 and sensitivity analysis](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- In the real PPIL4 system, strong in-vitro recruitment/ubiquitylation coexists with only modest incomplete cellular degradation, showing why a model that separates recruitment parameters from an effective downstream degradation parameter is directly relevant. — [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w)

### Inferences
- High-value project outputs would be: (i) whether each candidate is recruitment-limited (`K_D2/alpha`), CRBN-binding-limited (`K_D1`), exposure-limited, or processing-limited (`k_cat`); (ii) predicted dose/time windows for HiBiT and immunoblot experiments; (iii) expected `D_max` versus PPIL4/CRBN abundance across cell lines; (iv) how much improvement in `alpha` or `K_D1` is likely to matter before a `k_cat`/abundance ceiling is reached; and (v) whether competitor neosubstrates can account for weak or cell-line-dependent degradation. — [Chae et al.](https://doi.org/10.1021/acs.jcim.5c03125); [Sperling et al.](https://doi.org/10.1182/blood.2019000789)
- The best initial use is prospective experimental design: simulate plausible parameter ranges, identify which measurement most reduces uncertainty, collect that measurement, then refit. This adds mechanistic depth without claiming unmeasured affinity or efficacy. — [Chae et al., sensitivity/limitations](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- A rigorous ISEF comparison could ask whether the calibrated model predicts held-out dose–time points or a held-out candidate better than apparent TR-FRET EC50 alone; this is a testable contribution rather than a narrative overlay. — [Couturier et al.](https://doi.org/10.1021/jacs.5c09166); [Chae et al.](https://doi.org/10.1021/acs.jcim.5c03125)

### Gaps
- The model is not yet a validated PPIL4 predictor. Prospective utility must be demonstrated by fitting some PPIL4 measurements and predicting held-out conditions or compounds.

## What can it not infer from docking, and how should hook effect and kinetics be interpreted?

### Takeaway
HADDOCK/Vina scores can support structural plausibility and pose selection only. They are not `K_D`, cooperativity, residence time, ubiquitination rate, cellular exposure, `D_max`, or `DC50`; substituting them into the ODE model would create false precision.

### Cited Findings
- Chae et al. explicitly describe their required inputs as dissociation constants, protein expression levels, and degradation rates, and state that their current framework omits systemic PK, distribution, metabolism, and elimination. — [Chae et al., Discussion](https://doi.org/10.1021/acs.jcim.5c03125)
- The model collapses ubiquitination and subsequent degradation into a single first-order `k_cat`; it therefore cannot distinguish lysine-geometry/ubiquitin-transfer failure from proteasomal-processing failure without extending the state diagram. — [Chae et al., limitations](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125)
- The PPIL4 paper’s own docking interpretation for Z6466608628 is explicitly a potential structural explanation, whereas recruitment, ubiquitylation, and cellular depletion were established with separate biochemical/cellular experiments. — [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w)
- Classical high-dose hook effects are most characteristic of bifunctional degraders that independently saturate both proteins. Typical monovalent molecular glues preferentially bind one partner and often show saturating rather than bell-shaped recruitment; a complete three-body model can nevertheless represent alternative paths or measurable weak binding to the second partner, so the concentration range should still be tested empirically. — [Couturier et al.](https://doi.org/10.1021/jacs.5c09166); [Douglass et al., ternary-complex equations](https://doi.org/10.1074/jbc.RA120.014715)
- For ternary systems with two appreciable binary ligand-binding paths, exact mass-action solutions predict bell-shaped ternary-complex curves, with positive cooperativity broadening the useful concentration range; in the limiting molecular-glue case where one binary affinity approaches zero, the response becomes a conventional saturating two-step process. — [Douglass et al., JBC 2020](https://doi.org/10.1074/jbc.RA120.014715)

### Inferences
- The project should never map HADDOCK score differences to `K_D` or `alpha` without a separately validated calibration set, and even such a calibration would not identify `k_cat`, PPIL4 half-life, CRBN abundance, or intracellular free concentration. — [Chae et al.](https://doi.org/10.1021/acs.jcim.5c03125); [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w)
- A dose range extending well above the expected recruitment midpoint remains useful: absence of a hook supports the one-anchor molecular-glue mechanism; a declining high-dose signal would indicate either a second binding path, assay interference/aggregation, or other nonminimal behavior requiring investigation. — [Couturier et al.](https://doi.org/10.1021/jacs.5c09166); [Douglass et al.](https://doi.org/10.1074/jbc.RA120.014715)
- Structural models are most useful upstream to define mutations and controls—e.g., PPIL4 G-loop residues, CRBN pocket/interface residues, or the removed PPIL4-facing arm—that can experimentally test which reaction step changed. — [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w)

### Gaps
- No primary study was found that converts HADDOCK scores for this PPIL4 series into calibrated thermodynamic or kinetic parameters. Any such conversion would currently be unsupported.
- The Chae framework does not include explicit pharmacokinetics, cell permeability, spliceosome sequestration, multistep ubiquitin-chain formation, proteasomal queueing, or downstream RAS/spliceosome phenotypes; those require additional model states and data. — [Chae et al., Discussion/limitations](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03125); [Primary PPIL4 paper](https://doi.org/10.1038/s41467-025-62099-w)

## Which primary studies are the strongest precedents?

### Takeaway
The strongest chain of evidence is: exact ternary-equilibrium mathematics → CRBN substrate/competition experiments → PPIL4-specific recruitment/ubiquitylation/degradation → molecular-glue thermodynamic characterization → the 2026 cellular kinetic model.

### Cited Findings
- **Douglass et al., JBC 2020:** exact analytical solutions for ternary-complex formation, cooperativity, and hook behavior, applicable to high-cooperativity molecular-glue limits. — [Primary model](https://doi.org/10.1074/jbc.RA120.014715)
- **Sperling et al., Blood 2019:** CRBN substrate affinity, ligase abundance, competitive neosubstrates, and degradation kinetics determine cellular activity/resistance of thalidomide analogues. — [Primary experimental study](https://doi.org/10.1182/blood.2019000789)
- **Yamanaka et al., Nature Communications 2023:** CRBN-binding lenalidomide derivatives alter neosubstrate selectivity and provide matched AlphaScreen/HiBiT-like data used by the 2026 model. — [Primary experimental study](https://doi.org/10.1038/s41467-023-40385-9)
- **Couturier et al., JACS 2025:** a complete molecular-glue equilibrium framework that extracts binding/cooperativity values from biochemical matrix data. — [Primary model/experimental validation](https://doi.org/10.1021/jacs.5c09166)
- **Baek et al., Nature Communications 2025:** the project-defining CRBN–PPIL4 study establishes direct recruitment, 9DWV/structural evidence, in-vitro ubiquitylation, selective but incomplete cellular PPIL4 depletion, and Z6466608628 as an improved recruiter. — [Primary PPIL4 study](https://doi.org/10.1038/s41467-025-62099-w)
- **Chae et al., JCIM 2026:** the most directly relevant end-to-end biophysical kinetic framework, adding cellular target production/turnover, effective glue-induced degradation, abundance effects, and competitive neosubstrates. — [Primary model](https://doi.org/10.1021/acs.jcim.5c03125); [author code](https://github.com/jcd111/MG_Manuscript_Model)
- **Izaguirre et al., bioRxiv 2025 (GlueMap):** a separate structural/pharmacodynamic modeling preprint combining MD, binding free energy, ubiquitination geometry/metrics, and machine learning for molecular-glue discovery; useful as a structural complement but not a substitute for measured PPIL4 kinetic parameters, and all authors disclosed employment/shareholding in Atommap. — [Primary preprint](https://doi.org/10.1101/2025.01.13.632817)

### Inferences
- For an evidence-bound PPIL4 project, the most defensible integration is: stage-08 docking for pose/triage → PPIL4-specific ternary matrix for `K_D/alpha` → ubiquitination and HiBiT time courses for effective `k_cat`/degradation → Chae-model prediction of held-out conditions. — [Couturier et al.](https://doi.org/10.1021/jacs.5c09166); [PPIL4 study](https://doi.org/10.1038/s41467-025-62099-w); [Chae et al.](https://doi.org/10.1021/acs.jcim.5c03125)

### Gaps
- GlueMap remains a preprint in the source found and has declared commercial conflicts; its reported performance should not be treated as independent validation for PPIL4. — [GlueMap preprint](https://doi.org/10.1101/2025.01.13.632817)
