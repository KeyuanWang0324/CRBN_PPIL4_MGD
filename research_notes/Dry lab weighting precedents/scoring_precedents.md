# Evidence-Calibrated Multidimensional Scoring Precedents

## Which projects and papers are closest to the proposed framework?

### Takeaway
The affinity-prediction core of the idea is an established field: empirical and machine-learning scoring functions already learn relationships between structural/interaction descriptors and experimental affinity. The closest literal precedent is ΔVinaRF20, which learns a correction to the Vina score from 20 additional descriptors; ECIF plus ligand descriptors and several multimodal neural models are also close. A project framed only as “combine docking and structural scores with learned weights to predict affinity” would therefore have low methodological novelty.

### Cited Findings
- **RF-Score (2010)** used random forests to map counts of protein–ligand atom-type contacts to experimental binding affinity. It was trained on 1,105 PDBbind refined-set complexes and tested on the 195-complex core set, attaining reported Pearson R = 0.776 and RMSE = 1.58. This is a foundational example of replacing a hand-written scoring formula with weights/relationships learned from experimental data. — [Original paper](https://academic.oup.com/bioinformatics/article/26/9/1169/199938)
- **ΔVinaRF20 (2017)** is the closest direct precedent to the proposed “weighted combination”: it predicts a random-forest correction to the AutoDock Vina affinity estimate using 20 features (10 Vina-derived interaction terms and 10 buried-SASA terms). The training set contained 3,336 crystal structures plus 3,322 computer-generated decoys; it was evaluated for scoring, ranking, docking, and screening power on CASF-2007 and CASF-2013. It is a post-scoring method, not a docking sampler. — [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5140681/)
- **Pafnucy (2018)** represents a complex on a 3D grid using atom type, hybridization, heavy-atom/heteroatom bonds, hydrophobicity, aromaticity, hydrogen-bond donor/acceptor status, ring membership, partial charge, and protein/ligand identity, then learns affinity with a 3D CNN. It trained on 11,906 PDBbind v2016 complexes, used 1,000 randomly selected refined-set complexes for validation, and tested on the PDBbind core set and 73 independent Astex complexes. — [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6198856/)
- **OnionNet (2019)** learns experimental affinity from protein–ligand element-pair contact counts divided into multiple distance shells, thereby encoding local and nonlocal structure without a conventional docking-energy term. It was tested on CASF-2013, the PDBbind v2016 core set, and docking-generated poses; reported R/RMSE on the v2016 core set were 0.816/1.278 pK units. — [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6776976/)
- **AK-Score (2020)** ensembles multiple 3D convolutional neural networks operating on voxelized protein–ligand structures. It evaluated both CASF-2016 and a temporal-style holdout of structures added to PDBbind 2018 but absent from PDBbind 2016; reported correlation with experiment was 0.827 versus 0.616 for AutoDock Vina and 0.650 for X-Score on its benchmark comparison. — [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7697539/)
- **ECIF (2021)** uses counts of extended-connectivity-defined protein/ligand atom-type pairs and gradient-boosted trees; the best model combined ECIF structural interaction features with ligand descriptors. On CASF-2016, reported Pearson R was 0.857 for ECIF alone and 0.866 with ligand descriptors. This is a close precedent for fusing structural and ligand-level dimensions. — [Original paper](https://academic.oup.com/bioinformatics/article/37/10/1376/5998664)
- **Structure-Based Deep Fusion Inference (2021)** explicitly fused heterogeneous 3D model representations for affinity prediction, trained/validated on PDBbind 2016, and created a PDBbind-2019 holdout using both a temporal cutoff and 3D structure-based clustering to assess novel protein–ligand complexes. — [Original paper](https://pubs.acs.org/doi/10.1021/acs.jcim.0c01306)
- **HAC-Net (2023)** fuses a channel-attention 3D CNN with two attention-based graph convolutional networks. Its evaluation included splits designed to maximize differences in protein structure, protein sequence, or ligand fingerprints, plus similarity-controlled cross-validation, making it especially relevant to rigorous generalization testing. — [Original paper](https://doi.org/10.1021/acs.jcim.3c00251)
- **Consensus scoring is related but not equivalent.** It combines rankings or scores from multiple docking/scoring programs, often after normalization, primarily to improve virtual-screening enrichment. The original consensus-scoring study analyzed two docking methods and 13 scoring functions; later EFO work learned target-specific linear combinations optimized for enrichment. These methods do not necessarily learn calibrated experimental affinity or cellular performance. — [Original consensus-scoring paper](https://pubs.acs.org/doi/10.1021/jm990352k); [EFO linear-combination study](https://doi.org/10.3390/ijms20092060)
- **Target-specific ML scoring is also established.** A current protocol explicitly covers preparation of target-specific experimental data, train/test partitioning, descriptor generation, and evaluation of ML scoring functions for discriminating actives from inactives. — [Nature Protocols](https://www.nature.com/articles/s41596-023-00885-w)

### Inferences
- The project should not claim that the act of combining docking, interaction, and ligand descriptors with learned weights is new; RF-Score, ΔVinaRF20, ECIF and deep-fusion models occupy that space.
- A simple linear weighted score could still be educational and interpretable, but its contribution would be implementation, curation, interpretability, or validation rather than a new algorithmic principle.
- The most defensible prototype is a single-target or target-family model with one precisely defined experimental endpoint. This reduces biological and assay heterogeneity while leaving enough room to test whether each feature family adds value.

### Gaps
- Most published models do not expose a small, human-readable set of stable “weights”; random forests and neural networks learn nonlinear interactions. A sparse linear or monotonic model would be more interpretable but might sacrifice accuracy.
- Published benchmarks rarely establish that improved affinity prediction translates into cellular potency, selectivity, toxicity, pharmacokinetics, or clinical performance.

## How have similar models been validated, and what should count as convincing validation?

### Takeaway
Random or conventional PDBbind core-set tests are useful for comparison but can substantially overestimate generalization. Strong evidence requires separation by ligand chemistry and target/pocket similarity, a temporal or publication holdout, a fully external dataset, and ideally a blind prospective test.

### Cited Findings
- CASF-2016 evaluates four distinct capabilities: scoring power (correlation with affinity), ranking power (ordering ligands for the same target), docking power (selecting native-like poses), and screening power (recognizing target binders). High performance on one is not proof of high performance on the others. — [CASF-2016 paper](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00545)
- The 2017 similarity analysis of ML scoring functions showed that protein structural and sequence similarity between train and test sets strongly affected results; after removing highly similar training proteins, the tested ML methods no longer outperformed conventional scoring functions. — [Original study](https://pubmed.ncbi.nlm.nih.gov/28358210/)
- LP-PDBbind found cross-contamination by similar proteins and ligands between standard PDBbind subsets and argued that even a pure time split does not eliminate leakage. It rebuilt train/validation/test sets using ligand and protein similarity controls and added BDB2020+, an independent post-2020 BindingDB/PDB test set. — [LP-PDBbind paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10462179.1/)
- In LP-PDBbind testing, original ML models showed larger generalization problems than their conventional core-set results implied; retraining on cleaned, similarity-controlled splits improved 3D methods on the independent BDB2020+ set. The paper also evaluated same-target ranking on 40 Mpro and 23 EGFR complexes. — [LP-PDBbind results](https://pmc.ncbi.nlm.nih.gov/articles/PMC10462179.1/)
- Pafnucy removed exact overlaps and used an external Astex set, but its validation selection was random and exact de-duplication does not exclude close protein or ligand analogues across splits. — [Pafnucy methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6198856/); [later similarity-controlled critique](https://pubmed.ncbi.nlm.nih.gov/28358210/)
- D3R Grand Challenges provide blind, prospective pose and affinity-ranking tests. In Grand Challenge 3, machine-learning/structure hybrids did not outperform other approaches overall; pose accuracy did not correlate well with affinity-prediction accuracy, and supplying crystal poses often produced little or no ranking improvement. — [D3R Grand Challenge 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6472484/)
- A D3R analysis reported best prospective affinity-ranking correlations around Spearman rho 0.5 in the studied systems and showed that receptor selection could substantially change rankings, illustrating that even sophisticated retrospective performance does not guarantee robust use. — [Prospective D3R validation](https://pmc.ncbi.nlm.nih.gov/articles/PMC5771500/)

### Inferences
- Recommended evaluation hierarchy: (1) nested cross-validation grouped jointly by Bemis–Murcko scaffold and target/pocket cluster; (2) leave-one-publication or leave-one-assay-series-out testing; (3) time-split external testing; and (4) a preregistered prospective wet-lab test.
- For a single-target prototype, a ligand-scaffold split plus a chronological/publication holdout is more relevant than a random split. For a multi-target model, both ligand and protein/pocket similarity must be controlled.
- The model should be compared with at least four baselines: docking score alone, structure/contact features alone, ligand descriptors alone, and a simple similarity/QSAR baseline. Without those ablations, apparent gains cannot be attributed to multimodal integration.
- Uncertainty intervals and an applicability-domain flag should accompany every prediction; a point estimate outside the chemical/target domain invites false precision.

### Gaps
- A scaffold split alone is not a complete leakage safeguard because different scaffolds may still be highly similar; explicit similarity thresholds or clustering should supplement scaffold labels.
- Public blind challenges test only selected targets and endpoints, so prospective performance can remain target-dependent.

## What methodological pitfalls are documented?

### Takeaway
The dominant risks are leakage, assay-label heterogeneity, learning ligand or target identity instead of interactions, conflating pose/affinity/activity, and combining incompatible docking scores. These risks can make a sophisticated model look accurate while offering little value on genuinely new compounds.

### Cited Findings
- PDBbind contains Kd, Ki, and IC50-derived labels, but IC50 is not directly convertible to binding free energy without assay and inhibition-mechanism information. LP-PDBbind notes that IC50 for the same protein–ligand system can vary by up to an order of magnitude across assays and defines a cleaner level retaining exact Kd or Ki values. — [LP-PDBbind data-quality analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10462179.1/)
- BindingDB records solution composition, pH, and temperature because they affect measured affinities. This supports keeping assay metadata rather than treating literature values as context-free ground truth. — [BindingDB database paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC1751547/)
- A statistical analysis of mixed IC50 data found comparability depends on matched assay type, cell line, substrate, and binding context, underscoring the danger of pooling heterogeneous literature measurements. — [IC50 comparability study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3628986/)
- Latent-bias studies have shown that protein-only and ligand-only models can obtain substantial apparent performance on popular affinity datasets, demonstrating that high benchmark correlation need not mean the model learned the physical interaction. — [Original latent-bias study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9850481/)
- Docking pose quality and affinity-ranking quality are different tasks: D3R GC3 found little correlation between the accuracy of poses used and the resulting affinity predictions. — [D3R GC3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6472484/)
- Consensus methods must normalize heterogeneous scores because docking programs use different units and ranges; common transformations include ranks, min–max scaling, and z-scores. This means raw outputs should not simply be averaged or interpreted as commensurate physical energies. — [Consensus-scoring analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC9941253/)
- Traditional affinity models are not models of whole-drug performance. The experimental endpoint must remain explicit: affinity, biochemical inhibition, cellular potency, degradation, permeability, selectivity, and toxicity represent different outcomes and causal stages. D3R separately evaluates affinity ranking and active/inactive classification, and CASF separates scoring, ranking, docking, and screening. — [D3R GC3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6472484/); [CASF-2016](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00545)

### Inferences
- All computational features should be regenerated with one locked preparation/docking protocol. Published docking scores from different programs, versions, structures, grids, protonation states, and exhaustiveness settings should not be pooled as if they were one feature scale.
- Kd, Ki, IC50, cellular EC50, and degradation endpoints should not be placed in one response column. Either train separate endpoint-specific models or use a hierarchical/multitask model that retains assay identity and context.
- Repeated measurements should be preserved with provenance and uncertainty, not silently averaged. Paper, laboratory, assay format, species, construct, pH, temperature, substrate concentration, cell line, and exposure time can be modeled as covariates or used to define comparable strata.
- Model interpretation should use ablation, permutation importance or SHAP only after leakage-controlled validation; feature importance from a leaked model can explain the shortcut rather than the mechanism.

### Gaps
- The literature does not provide a universally transferable conversion from docking/structural scores to cellular efficacy.
- Many papers report only conventional CASF/PDBbind validation, so direct comparison of strict out-of-distribution performance remains incomplete.

## What is already established, and where could this project still be creative and significant?

### Takeaway
As a new affinity-scoring algorithm, the project has low-to-moderate creativity because many close precedents exist. It can become meaningfully creative as a transparent, assay-aware, uncertainty-calibrated framework that predicts several clearly separated endpoints, quantifies transfer failure, and is tested prospectively rather than merely achieving another retrospective correlation.

### Cited Findings
- Empirical scoring functions already consist of interaction descriptors calibrated against experimental affinity using regression; modern ML scoring functions generalize this idea with nonlinear models. — [Review of empirical scoring functions](https://pmc.ncbi.nlm.nih.gov/articles/PMC6165880/)
- ΔVinaRF20 already combines a docking score with additional interaction and SASA descriptors, while ECIF already combines protein–ligand interaction counts with ligand descriptors. These substantially overlap the proposed technical core. — [ΔVinaRF20](https://pmc.ncbi.nlm.nih.gov/articles/PMC5140681/); [ECIF](https://academic.oup.com/bioinformatics/article/37/10/1376/5998664)
- Blind D3R results show that improved retrospective metrics and correct-looking poses do not automatically deliver accurate affinity rankings, so rigorous external or prospective validation is itself scientifically important. — [D3R GC3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6472484/)
- The availability of leak-proof splits and independent post-2020 tests makes it possible to evaluate whether an apparently strong method actually transfers to dissimilar protein–ligand systems. — [LP-PDBbind](https://pmc.ncbi.nlm.nih.gov/articles/PMC10462179.1/)

### Inferences
- **Creativity rating if left broad:** approximately **3/10**. “Learn weights for docking, structural, and ligand features from published wet-lab values” essentially describes an empirical/ML scoring function.
- **Creativity rating with a strong redesign:** approximately **6–7/10** for a student research project. The distinctive package would be: interpretable feature families; assay-context modeling; separate affinity, biochemical, cellular, and possibly degradation heads; calibrated uncertainty/applicability domain; leakage-resistant splits; and a preregistered prospective test.
- **Scientific significance:** potentially high for candidate triage and reproducibility, but only moderate if the output is another single retrospective “overall drug score.” A credible negative result—showing which dry-lab dimensions fail to transfer—would itself be valuable.
- A strong project question is: **“Under leakage-resistant and assay-aware validation, do structure-derived features add predictive information beyond ligand similarity and conventional docking for a defined experimental endpoint?”** This directly tests scientific value instead of presuming that more features will help.
- A second, more original question is: **“Can a calibrated multi-endpoint model predict affinity, biochemical activity, and cellular activity separately while learning where structure-derived evidence ceases to transfer?”** The output should be a profile with uncertainty, not an opaque universal drug-success number.
- The prototype should begin with one data-rich target or coherent target family and one primary endpoint (preferably Kd or Ki). Additional endpoints should be added only where matched compound-level measurements and assay metadata exist.

### Gaps
- No source located in this search demonstrates a broadly validated model that converts docking/structure features into a universally meaningful, multidimensional measure of real-world “drug performance.” This gap is partly scientific opportunity and partly a sign that the requested endpoint is causally too broad.
- The feasibility and novelty of a specific prototype will depend on the chosen target, the number of matched compounds, access to assay metadata, and whether a genuinely unseen prospective test set can be obtained.
