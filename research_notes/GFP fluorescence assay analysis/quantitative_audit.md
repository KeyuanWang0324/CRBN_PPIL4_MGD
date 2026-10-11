# Quantitative audit of the GFP fluorescence plate

## Do the workbook calculations reconcile?

### Takeaway
Yes. An independent reconstruction from the four `Long Data` sheets reproduced every value in `Per-Well Normalized` to floating-point precision (maximum absolute difference `2.22e-16`) and reproduced all means and sample SDs in `Final Results`. The numerical pipeline is internally correct for the documented method, but the final workbook is not dynamically linked to the four raw workbooks.

### Cited Findings
- Within each read, the workbook method uses the median of `A2-A6, H2-H6, B1-G1` as the 293FT-GFP background and the median of `A7-A11, H7-H11, B12-D12` as the HeLa-GFP background; it then subtracts that background, divides by the mean corrected signal in `B2-G2` or `B11-G11`, averages each physical treatment well across the four reads, and finally reports mean and `STDEV.S` across four treatment wells. Evidence: `2026-10-06_Final_Normalized_Fluorescence_Results_293FT-GFP_HeLa-GFP.xlsx`, `Method!A3:C14`; the four processed workbooks, `Long Data!A1:N97`.
- The six final rows contain only `AVERAGE` and `STDEV.S` formulas over the hardcoded `Per-Well Normalized` values. The per-well values themselves do not link back to the raw files. Evidence: final workbook, `Final Results!A7:H12` and `Per-Well Normalized!A1:G13`.
- E12-G12 were excluded in all calculations. They are present but flagged in IMG_1123 and IMG_1128, and blank/omitted in IMG_1124 and IMG_1126. Evidence: the four processed workbooks, `Long Data`, records E12-G12.
- The reconstructed pooled-control results are:

| Cell line | Treatment | 5 µM (% control, mean ± SD) | 10 µM (% control, mean ± SD) |
|---|---|---:|---:|
| 293FT-GFP | Z64 | 69.9 ± 9.8 | 86.3 ± 6.3 |
| 293FT-GFP | TM | 97.3 ± 5.1 | 93.8 ± 5.0 |
| 293FT-GFP | Negative control | 95.2 ± 6.7 | 116.4 ± 11.4 |
| HeLa-GFP | Z64 | 75.7 ± 6.1 | 94.9 ± 18.7 |
| HeLa-GFP | TM | 103.6 ± 22.1 | 86.0 ± 4.6 |
| HeLa-GFP | Negative control | 114.8 ± 25.7 | 112.4 ± 31.3 |

### Inferences
- The workbook contains no arithmetic error in the stated pipeline.
- Reproducibility would improve if a master analysis workbook or script calculated normalized values directly from the four source tables. At present, a corrected raw value would not propagate automatically to the final workbook.

### Gaps
- The raw instrument export was not available, so transcription accuracy from the images into the four processed workbooks could not be independently verified.
- The precise reader settings, gain mode, excitation/emission wavelengths, shaking, focal height, and read direction are not documented in these workbooks.

## How reproducible are the four sequential reads, and are they independent replicates?

### Takeaway
The plate's spatial ranking is reproducible, but the absolute RFU scale changes dramatically between scans. The four scans are repeated measurements of the same wells, not four independent experiments, so they can reduce readout noise but cannot increase biological sample size.

### Cited Findings
- Read times were 20:19:49, 20:22:33, 20:28:17, and 20:30:15 on 2026-10-06; recorded temperature fell from 37.0°C to 33.7°C. Evidence: each processed workbook, `Plate Matrix!A4:B5` and `Long Data!L:N`.
- Across the 93 mutually usable non-bubble wells, raw RFU medians were 7,007; 25,544; 23,537; and 7,857. Pairwise Pearson correlations were nevertheless 0.981-0.988 (Spearman 0.953-0.966), indicating strong preservation of relative well ranking despite roughly 3.6-fold scale shifts. Evidence: independent calculation from the four `Long Data!A1:N97` tables.
- After background and pooled-control normalization, pairwise Pearson correlations across the 48 treatment wells were 0.806-0.905. The median within-well CV across the four reads was 7.2% (IQR 5.0-10.0%; maximum 19.2%). Evidence: independent calculation from the four `Long Data` tables.
- Across-read CVs of condition means were 1.6-3.9% for 293FT-GFP and 3.5-9.1% for HeLa-GFP. These values describe repeat-read stability only.
- Repeated measurements of the same experimental units do not increase `N`; treating them as independent inflates evidence and is pseudoreplication. [Lazic et al., PLOS Biology](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2005282)

### Inferences
- The large raw scale jumps are most consistent with an undocumented reader/gain or measurement-state change; the workbooks alone cannot identify the cause. Per-read normalization is essential and raw RFUs from different reads should not be pooled directly.
- Averaging the four reads for each physical well is a reasonable descriptive noise-reduction step. The effective biological experiment count remains one plate/run.

### Gaps
- It is unknown whether auto-gain or other settings changed between scans.
- It is unknown whether the same cell preparation, passage, and treatment preparation were used for all wells, although the layout strongly indicates one experimental run.

## Which findings survive pooled-control versus same-row-control normalization?

### Takeaway
Only a few directional findings are stable. Z64 is below control at both doses in 293FT-GFP and at 5 µM in HeLa-GFP under both methods. The apparent HeLa TM 10 µM decrease and several HeLa negative-control elevations are normalization-dependent. There is no robust monotonic dose-response pattern.

### Cited Findings
- Same-row normalization used the no-drug well in the same plate row (`B2/B11` for Z64 5 µM, `C2/C11` for Z64 10 µM, etc.), calculated separately in each read after the same background subtraction.

| Cell line | Treatment | Dose | Pooled-control mean ± SD (%) | Same-row mean ± SD (%) | Assessment |
|---|---|---:|---:|---:|---|
| 293FT-GFP | Z64 | 5 µM | 69.9 ± 9.8 | 83.0 ± 11.8 | Direction robust; magnitude control-sensitive |
| 293FT-GFP | Z64 | 10 µM | 86.3 ± 6.3 | 86.8 ± 6.3 | Direction and magnitude robust |
| 293FT-GFP | TM | 5 µM | 97.3 ± 5.1 | 84.3 ± 4.3 | Not robust |
| 293FT-GFP | TM | 10 µM | 93.8 ± 5.0 | 92.1 ± 5.0 | Small directional decrease under both |
| 293FT-GFP | Negative control | 5 µM | 95.2 ± 6.7 | 94.6 ± 6.7 | Near control under both |
| 293FT-GFP | Negative control | 10 µM | 116.4 ± 11.4 | 119.2 ± 11.6 | Elevated under both |
| HeLa-GFP | Z64 | 5 µM | 75.7 ± 6.1 | 84.3 ± 6.6 | Direction robust; magnitude control-sensitive |
| HeLa-GFP | Z64 | 10 µM | 94.9 ± 18.7 | 79.3 ± 15.7 | Below control, but magnitude and dose ordering fragile |
| HeLa-GFP | TM | 5 µM | 103.6 ± 22.1 | 109.5 ± 23.2 | No supported decrease |
| HeLa-GFP | TM | 10 µM | 86.0 ± 4.6 | 108.0 ± 5.7 | Direction reverses; not robust |
| HeLa-GFP | Negative control | 5 µM | 114.8 ± 25.7 | 96.1 ± 21.6 | Direction changes; not robust |
| HeLa-GFP | Negative control | 10 µM | 112.4 ± 31.3 | 122.1 ± 34.2 | Elevated mean but highly unstable |

- For Z64, dose ordering is not robust: in 293FT-GFP the large pooled difference between 5 and 10 µM nearly disappears with row matching; in HeLa-GFP, row matching reverses which dose has the lower mean.
- Fluorescence assays can produce false decreases through absorbance/quenching and false increases through compound autofluorescence; counter-screens are needed to distinguish optical interference from biological activity. [NIH Assay Guidance Manual](https://www.ncbi.nlm.nih.gov/books/NBK343429/?report=printable)

### Inferences
- The best current descriptive signal is Z64-associated lower fluorescence, particularly 293FT-GFP at 5 µM and HeLa-GFP at 5 µM. This remains a preliminary within-plate observation, not evidence of PPIL4 degradation, toxicity, or potency.
- The 293FT-GFP negative-control increase at 10 µM could reflect biology or compound-dependent optical interference. The experiment cannot separate those mechanisms.

### Gaps
- There are no compound-only/no-cell wells, GFP-negative parental cells, or independent cell-count/viability measurements to test autofluorescence, quenching, or differences in cell number.
- The biological meaning of total-well GFP fluorescence is not documented. It is unclear whether GFP is a target-linked reporter, a constitutive marker, or a proxy for cell number.

## Are there plate-position effects, outliers, heteroscedasticity, or QC failures?

### Takeaway
Yes. Treatment identity and dose are completely confounded with plate row, and the no-drug controls show persistent row differences large enough to alter conclusions. HeLa-GFP also contains influential high wells and markedly larger technical variance than 293FT-GFP.

### Cited Findings
- Averaged across reads after within-read normalization, 293FT-GFP control wells ranged from 84.6% (B2) to 115.5% (D2), a 30.9-percentage-point span. HeLa-GFP controls ranged from 81.1% (E11) to 120.9% (C11), a 39.8-point span. Evidence: independent calculation from the four raw workbooks.
- Because every condition occupies exactly one row, this row variation is inseparable from treatment/dose without matched-row normalization or a randomized/block-balanced layout.
- The 293FT-GFP Z64 5 µM wells increase monotonically from B3 to B6: 57.4%, 67.6%, 74.6%, and 80.2% of pooled control. Across all six 293FT-GFP conditions, the first treatment column position averages 84.5%, versus 95.8%, 97.8%, and 94.6% for positions 2-4. This suggests a column-position or dispensing/seeding component.
- HeLa-GFP has several influential wells. For negative control 5 µM, F7 is 152.8% while the other three are 96.8-105.2%; deleting F7 changes the mean from 114.8% to 102.1%. For negative control 10 µM, G8 is 157.8%; deleting it changes the mean from 112.4% to 97.2%. These wells are not statistical proof of error, but the treatment claims are not robust to them.
- HeLa-GFP TM 5 µM is split into two low wells (82.2%, 86.9%) and two high wells (122.5%, 122.8%), consistent with a positional pattern rather than a homogeneous treatment effect.
- Across physical treatment wells, pooled-normalization CVs were 5.2-14.1% for 293FT-GFP, but 5.4-27.9% for HeLa-GFP. The most variable HeLa groups were negative control 10 µM (27.9%), negative control 5 µM (22.4%), TM 5 µM (21.3%), and Z64 10 µM (19.7%).
- Per-read control CV was 6.0-15.0% for 293FT-GFP and 10.5-27.0% for HeLa-GFP. The NIH Assay Guidance Manual uses CV ≤20% plus absence of material spatial effects as common plate-uniformity acceptance criteria and recommends explicit spatial-uniformity studies; it also warns that row, column, and edge effects can arise from incubation, dispensing, and readers. [HTS Assay Validation](https://www.ncbi.nlm.nih.gov/sites/books/NBK83783/?report=classic); [Microplate recommended practices](https://www.ncbi.nlm.nih.gov/books/NBK558077/)
- B12-D12 are labeled cell-free but are consistently 27-31% brighter than the adjacent HeLa-side background wells. Because the analysis uses a median of 13 background wells, excluding B12-D12 changes final HeLa means by at most about 0.9 percentage points, so this anomaly is real but not materially driving the conclusions.

### Inferences
- Row-matched controls are a useful sensitivity analysis and likely more appropriate than a single pooled mean for this layout, but one control well per row is insufficient to estimate a row-specific control distribution.
- Standard one-way tests on the four wells per condition would confound drug/dose with row and would not solve the lack of independent experiments.
- A formal Z-prime value is not available because the plate lacks well-replicated positive/max and negative/min control states. A compound labeled “negative control” cannot substitute for an assay-state positive control when it itself changes fluorescence.

### Gaps
- The layout does not allow separation of drug effect from row, dispensing order, cell-density gradients, or local reader effects.
- The source images do not establish whether F7, G8, D8, D9, or other influential wells had bubbles, uneven cell density, or pipetting errors.

## Is statistical inference justified?

### Takeaway
No confirmatory biological inference is justified from this dataset alone. Descriptive summaries of the one plate are valid, but p-values or confidence intervals that treat wells or scans as independent biological replicates would be pseudoreplication.

### Cited Findings
- There is one plate/run, four wells per condition, and four sequential measurements of those same wells. The four wells quantify within-plate technical variability; the repeated scans quantify readout repeatability.
- Genuine replication requires the treatment to be applied to multiple independent experimental units. Repeated measurements of the same units do not increase sample size, and statistical analysis cannot repair this design after data collection. [Lazic et al., PLOS Biology](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2005282)
- Assay-validation guidance separates within-run variation from between-run variation and recommends replicate experiments plus plate-uniformity studies to characterize assay reproducibility. [NIH Assay Guidance Manual](https://www.ncbi.nlm.nih.gov/sites/books/NBK83783/?report=classic)

### Inferences
- Report these data as a pilot plate with `n = 4 technical wells per condition on one experimental run`; state that each well was averaged across four sequential reads.
- Do not report a biological p-value, DC50, Dmax, or confirmed hit from these two concentrations.
- If a purely descriptive uncertainty interval is shown, label it explicitly as within-plate technical variation, not biological uncertainty.

### Gaps
- There are no independent experiment days, cell preparations/passages, treatment preparations, or plates from which to estimate between-run biological and procedural variance.

## What calculations and experimental changes are justified next?

### Takeaway
The next analysis should preserve raw reads, summarize at the physical-well level, and then use independent experiment means as the inferential units. The repeat experiment must break the row-treatment confounding and add controls for cell number and optical interference.

### Cited Findings
- NIH assay guidance recommends plate-uniformity and replicate-experiment studies, inspection for row/column/edge drift, and on-plate controls. [HTS Assay Validation](https://www.ncbi.nlm.nih.gov/sites/books/NBK83783/?report=classic)
- Positional effects can arise from evaporation, temperature gradients, dispensing sequence, and the reader; uniformity plates and controls distributed across rows/columns are recommended diagnostic tools. [Microplate selection and recommended practices](https://www.ncbi.nlm.nih.gov/books/NBK558077/)
- Fluorescence interference can be assessed with counter-assays designed to detect autofluorescence and quenching. [NIH Assay Guidance Manual](https://www.ncbi.nlm.nih.gov/books/NBK343429/?report=printable)

### Inferences
- Minimum defensible repeat design:
  1. Run at least three independent experiments on different days, preferably using separately prepared cell suspensions and drug dilutions; define the experiment/run as `N`.
  2. Randomize or block each treatment/dose across multiple rows and columns. Place multiple vehicle/no-drug controls in every block or row rather than one control per row.
  3. Include a known positive assay control and enough max/min controls to assess dynamic range and Z-prime.
  4. Use a broader dose series (for example 6-8 concentrations with appropriate spacing) before attempting a dose-response fit.
  5. Fix and record reader gain, wavelengths, focal height, temperature/equilibration, shaking, and scan order; use one prespecified endpoint read unless repeated-read kinetics are part of the design.
  6. Add compound-only/no-cell wells at every dose, GFP-negative parental cells, and a parallel cell-count or viability readout. If imaging is available, analyze GFP per cell as well as total-well GFP.
  7. Inspect and document bubbles before reading; retain explicit exclusion flags and prespecified exclusion rules.
- Recommended analysis hierarchy for the repeat study: background correct within plate; normalize to block/row-matched vehicle controls; average technical reads within each well; summarize technical wells within each independent experiment; then perform statistical analysis across independent experiment means with experiment/block represented in the model.
- For the current pilot, add median and leave-one-well-out sensitivity summaries. These clearly show which HeLa means are driven by one well.

### Gaps
- The correct positive control and the expected direction of fluorescence change depend on the reporter biology and must be defined before the next plate is designed.
- Power/sample-size planning requires an expected biologically meaningful effect and an estimate of between-experiment variance, neither of which is available from this single run.
