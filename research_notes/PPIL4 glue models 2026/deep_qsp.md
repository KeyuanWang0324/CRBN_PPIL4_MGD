# 2026 deep-learning QSP and hybrid neural–mechanistic modeling

## Is “the 2026 deep-learning QSP model” a specific model or a general class?

### Takeaway
There is an exact 2026 paper matching that phrase: Goo et al., *Proteolysis-targeting Chimera efficacy prediction using a deep-learning–QSP model* (Journal of Cheminformatics, 12 January 2026). It is a specific PROTAC “Deep Hook” implementation, not a universal QSP model and not a model developed for molecular glues.

### Cited Findings
- Goo et al. combine DeepCalici-predicted binary dissociation constants with a mechanistic Hook model and a small DNN that scales Hook-model parameters to predict PROTAC `DC50` and `Dmax` — [peer-reviewed paper, DOI 10.1186/s13321-026-01152-2](https://doi.org/10.1186/s13321-026-01152-2).
- The curated dataset contains 444 in-vitro measurements, 272 PROTACs, 53 protein targets, two E3 ligases (CRBN and VHL), 95 warheads, 210 linkers and 13 E3 ligands; data were split 80:20 by individual measurement rather than held-out compound, scaffold or target — [paper full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888238/).
- Test performance was `R² = 0.547` for `DC50`, but `R² = -0.180` for `Dmax`; the authors attribute the poor `Dmax` generalization to unstructured/unstandardized context such as cell-line E3 abundance, proteasome availability, assay duration and detection method — [paper full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888238/).
- The released implementation expects PROTAC-, linker-, warhead- and E3-ligand descriptors plus warhead/E3-ligand binding energies, `KD,E`, `KD,P`, and target turnover `kdeg,P`; the repository documents local Docker/GPU reproduction and new-row prediction — [official author repository](https://github.com/swgoo/protac_deep_qsp).
- The model code implements closed-form Hook equations for `Dmax` and `DC50`, with a feed-forward network outputting multiplicative coefficients on `E0`, cooperativity, catalytic degradation, target turnover and the two binary `Kd` values — [official `model.py`](https://github.com/swgoo/protac_deep_qsp/blob/main/model.py).
- “DeepQSP” can also refer to an unrelated neural classifier of quorum-sensing peptides, so the name alone is ambiguous outside the targeted-protein-degradation context — [DeepQSP paper](https://doi.org/10.3389/fgene.2024.1301175).

### Inferences
- The user’s phrase almost certainly refers to Goo et al.’s 2026 PROTAC paper because its title is an exact match.
- The reported random measurement split is a weak test of out-of-chemotype or new-target generalization: repeated observations and shared PROTAC/warhead/linker chemotypes can occur on both sides. A PPIL4 application should require grouped splits by compound and scaffold, and ideally a held-out target test, before treating the model as prospective.
- The model cannot be transferred to PPIL4 molecular glues by simply entering a SMILES string. Molecular glues do not have separate warhead/linker/E3-ligand components, and the relevant ternary recruitment mechanism is not equivalent to two independent PROTAC binary affinities plus a single cooperativity correction.

### Gaps
- The article does not report prospective experimental validation on newly synthesized compounds or a held-out molecular-glue dataset.
- No evidence was found that the released model was trained on PPIL4, on PPIL4–CRBN molecular glues, or on monovalent molecular glues generally.

## Can the 2026 PROTAC Deep Hook model be useful for the PPIL4–CRBN molecular-glue project?

### Takeaway
It is useful mainly as a methodological reference and reproducibility/negative-control exercise, not as a defensible PPIL4 ranking engine. The project should adapt a molecular-glue-specific kinetic ODE model and generate PPIL4-specific time-course data before considering any neural augmentation.

### Cited Findings
- The 2026 model itself states that full QSP models are impractical in early design without extensive parameterization data, and its workaround depends on PROTAC-specific structural decomposition and predicted binary affinities — [Goo et al.](https://doi.org/10.1186/s13321-026-01152-2).
- A 2026 primary molecular-glue paper developed a unified kinetic framework covering glue stabilizers and CRBN glue degraders, explicitly adding the degradation step for the CRBN case; this mechanism is directly closer to PPIL4–CRBN than a bifunctional PROTAC Hook model — [Zaidman et al., *Biophysical Modeling Elucidates Mechanistic Principles for Rational Molecular Glue Design*](https://pubmed.ncbi.nlm.nih.gov/41973479/).
- The established PROTAC mechanistic-PD framework separates ternary-complex/target engagement, intrinsic protein turnover and catalytic degradation, and uses time-course degradation to estimate how degradation and recovery evolve — [Haid et al., 2023](https://doi.org/10.3390/ph16010098).
- In a separate open PROTAC activity model, training data include `pDC50`, `Dmax`, E3 type, target amino-acid sequence and cell type; performance fell from 0.865 ROC-AUC in the easiest test to 0.604 ROC-AUC when generalizing to novel targets, illustrating how strongly apparent performance depends on the split — [Ribes et al., 2024 preprint](https://arxiv.org/abs/2406.02637).

### Inferences
- Defensible current use: reproduce the published model locally, run known CRBN PROTAC controls, and document that its feature schema and training domain do not match PPIL4 glue candidates. This can support a Methods/Limitations discussion, but its outputs should not be merged with the locked PPIL4 ligand-inclusive docking rank.
- A custom PPIL4 model should start with states such as free glue, CRBN–glue, PPIL4–CRBN–glue ternary complex, ubiquitinated PPIL4 and total PPIL4. Candidate-specific quantities should enter only where measured or explicitly uncertain: intracellular free drug, glue–CRBN binding, ternary recruitment/cooperativity, ubiquitination/catalytic rate, and PPIL4 basal synthesis/turnover.
- This model would answer “what kinetic parameter ranges are consistent with observed PPIL4 degradation?” and “which experiment most reduces uncertainty?”, not “which docked pose is biologically active?”. It must remain downstream of and separate from structural docking.

### Gaps
- Public PPIL4–CRBN glue time courses, ternary-complex affinities/cooperativity values, PPIL4 half-life in the intended cell line and intracellular exposure measurements were not found in the searched literature.
- The transferability of the 2026 molecular-glue biophysical framework to PPIL4 specifically must be tested experimentally; a mechanistic resemblance is not validation.

## What data would be needed?

### Takeaway
For a small ISEF study, prioritize a compact, replicated concentration-by-time degradation dataset plus key controls and fix as many parameters as possible from direct measurements. Deep learning becomes credible only after assembling hundreds of independent, standardized compound–condition observations, not three docked candidates.

### Cited Findings
- Goo et al.’s neural-QSP model used 444 measurements from 272 compounds yet still failed to generalize `Dmax`, with the paper specifically identifying cell line, E3 abundance, proteasome availability, assay duration and detection method as missing context — [Goo et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888238/).
- The model obtained target turnover from protein half-life (`kdeg = ln(2)/t½`), substituted an average when unavailable, and initialized unavailable `E0`, `kcat` and cooperativity values to one, demonstrating which biological quantities must otherwise be measured or treated as uncertain — [Goo et al. methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888238/).
- Hybrid pharmacometric–ML reporting guidance identifies insufficient benchmarking, error propagation, feature-stability analysis, ablation, external validation and uncertainty quantification as recurring deficiencies, and recommends explicit estimands, convergence checks, diagnostics, sensitivity analysis, validation and verification — [Fochesato et al., 2026](https://doi.org/10.1002/psp4.70113).
- Hybrid neural ODEs can lose mechanistic parameter identifiability because flexible neural components compensate for mechanistic terms; Giampiccolo et al. therefore add explicit identifiability analysis and confidence intervals after fitting — [Giampiccolo et al., 2024](https://doi.org/10.1038/s41540-024-00460-3).

### Inferences
- Minimum useful PPIL4 experimental package: (1) vehicle and multiple log-spaced concentrations for each candidate; (2) PPIL4 protein remaining at several early, intermediate and late times with biological replicates; (3) washout/recovery measurements; (4) basal PPIL4 half-life/turnover; (5) CRBN and PPIL4 abundance in the chosen cell line; (6) CRBN-binding and ternary-recruitment/cooperativity data where feasible; (7) intracellular unbound exposure or at least stability/permeability measurements; and (8) CRBN-loss/competition, proteasome-inhibition and no-PPIL4-arm controls.
- Add matched assay metadata: cell line, passage/batch, treatment duration, detection method, replicate structure and censoring/limit-of-quantification. Without these, `Dmax` variation may be mistaken for chemistry.
- If only a few compounds are available, fit a hierarchical mechanistic model sharing basal PPIL4/CRBN parameters across compounds while allowing a small number of candidate-specific parameters. Avoid training a DNN on rows created by technical repeats; the independent unit is the compound/condition, not each measurement well.
- For any future ML dataset, validation must be group-based: leave out whole compounds, then scaffolds/chemotypes, and preferably targets or cell lines. Hold a prospective blinded set untouched until the analysis is frozen.

### Gaps
- There is no universal minimum sample size for neural-QSP models; it depends on feature dimension, task and effective number of independent compounds. The 444-measurement PROTAC example should be treated as a cautionary lower-bound case, not a guarantee of adequacy.

## What tasks are appropriate for deep learning or hybrid QSP?

### Takeaway
Appropriate roles are parameter inference, mechanistic-model emulation, virtual-population generation, dose/schedule exploration and uncertainty analysis after data exist. Structural candidate ranking and claims of PPIL4 degradation efficacy are not appropriate outputs from QSP alone.

### Cited Findings
- Machine-learning emulators have been used to accelerate virtual-population inference from an immuno-oncology QSP model, where repeated full-model simulation was the bottleneck — [Pawłowski et al., 2024](https://doi.org/10.1016/j.ymeth.2023.12.006).
- Simulation-based inference has been used to learn distributions over parameterizations of a large rheumatoid-arthritis QSP model for individual virtual patients, explicitly representing multiple parameter sets compatible with sparse/noisy patient data — [Paul et al., 2024](https://doi.org/10.3389/fsysb.2024.1444912).
- A 2025 neural emulator of a cardiac cellular QSP model achieved large simulation speedups but showed larger mismatches when experimental tissue data differed from the single-cell simulations used for training, illustrating domain-shift risk — [Grandits et al., 2025](https://doi.org/10.1016/j.vascn.2025.107850).
- Neural-ODE PK modeling has been evaluated for predicting new dosing regimens, while hybrid pharmacometric guidance categorizes ML roles as covariate learning, unknown-dynamics components and fast surrogate/emulator models — [Lu et al., 2021](https://doi.org/10.1016/j.isci.2021.102804); [Fochesato et al., 2026](https://doi.org/10.1002/psp4.70113).

### Inferences
- Near-term PPIL4 tasks: estimate kinetic parameters with profile-likelihood/Bayesian uncertainty; run global sensitivity analysis; design the next concentration/time point; compare candidate-vs-negative-control degradation dynamics; and explore whether a hook-shaped response is supported by data.
- Later tasks, only with PK/exposure and richer biological data: emulate expensive simulations, propagate cell-to-cell or cell-line variability, create virtual populations, and optimize schedules under explicitly stated efficacy/toxicity objectives.
- A surrogate is unnecessary if the compact PPIL4 ODE runs in milliseconds; it adds approximation error without solving a real bottleneck.

### Gaps
- No recent primary paper was found applying neural ODEs, simulation-based inference or virtual populations specifically to molecular-glue degraders.

## Should a conventional mechanistic ODE/QSP model come first?

### Takeaway
Yes. For the current PPIL4 project, a small molecular-glue-specific ODE with prespecified observables, sensitivity/identifiability analysis and uncertainty bounds is preferable. Add a neural component only after a reproducible residual pattern and enough independent data justify it.

### Cited Findings
- Mechanistic PROTAC models already describe target turnover and catalytic degradation with compact equations and can connect target degradation to downstream PD without requiring a neural network — [Haid et al., 2023](https://doi.org/10.3390/ph16010098).
- Hybrid neural ODE research identifies scarce/noisy data and partial observability as sources of non-identifiability and shows that neural flexibility can further compromise mechanistic parameter recovery — [Giampiccolo et al., 2024](https://doi.org/10.1038/s41540-024-00460-3).
- Current hybrid-model recommendations require comparison against simpler baselines, ablation of the ML component, uncertainty propagation and external validation, rather than assuming hybrid complexity improves inference — [Fochesato et al., 2026](https://doi.org/10.1002/psp4.70113).

### Inferences
- A staged plan is most defensible: **Stage 1**, pre-register a minimal glue–CRBN–PPIL4 kinetic model; **Stage 2**, collect replicated concentration/time courses and controls; **Stage 3**, fit only identifiable parameters and report intervals/predictive checks; **Stage 4**, validate on a held-out experiment; **Stage 5**, add a neural residual/covariate function only if it improves held-out prediction and leaves mechanistic parameters identifiable.
- With only docking scores and three candidate/control results, neither conventional QSP nor deep QSP can estimate degradation efficacy. The scientifically honest output is an experimental-design/sensitivity model, not a potency prediction.

### Gaps
- The exact ODE structure and parameter count should be chosen after the planned wet-lab readouts, assay times and controls are fixed; otherwise the model may contain states that are not observable.
