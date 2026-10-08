# Z64 fluorescence and Western-blot statistics and correlation

- Fluorescence experiment date: 2026-10-06
- Western-blot experiment date: 2026-10-06
- Statistical analysis date: 2026-10-08
- Compound: Z6466608628 (Z64)
- Cell models: 293FT-GFP, HeLa-GFP, and SY5Y

## Scope and inferential boundary

The fluorescence experiment contains one plate with four technical treatment wells per condition. The plate was read four times sequentially. The sequential reads are repeated measurements of the same physical wells and are not independent biological replicates. Treatment and dose are also confounded with plate row.

The Western blot contains one lane per dose in each cell model and no independent replicate blots. Consequently, valid biological significance cannot be calculated for either experiment. The fluorescence p-values below describe within-plate technical separation only. The Western-blot values are approximate densitometry and have no valid p-values.

## Fluorescence significance calculations

For each physical well, the four sequential normalized readings were averaged first. A two-sided Welch test then compared the four treatment wells with the six pooled vehicle wells for the corresponding cell line. Benjamini-Hochberg correction was applied across all 12 treatment-versus-control comparisons.

| Cell line | Treatment | Dose | Mean ± technical SD (% vehicle) | Difference from vehicle (percentage points) | Welch p | BH q | Result after correction |
|---|---|---:|---:|---:|---:|---:|---|
| 293FT-GFP | Z64 | 5 µM | 69.9 ± 9.8 | -30.1 | 0.0025 | 0.030 | Significant within this plate |
| 293FT-GFP | Z64 | 10 µM | 86.3 ± 6.3 | -13.7 | 0.0280 | 0.112 | Not significant after correction |
| 293FT-GFP | TM | 5 µM | 97.3 ± 5.1 | -2.7 | 0.5827 | 0.699 | Not significant |
| 293FT-GFP | TM | 10 µM | 93.8 ± 5.0 | -6.2 | 0.2297 | 0.459 | Not significant |
| 293FT-GFP | Negative control | 5 µM | 95.2 ± 6.7 | -4.8 | 0.3871 | 0.581 | Not significant |
| 293FT-GFP | Negative control | 10 µM | 116.4 ± 11.4 | +16.4 | 0.0577 | 0.173 | Not significant |
| HeLa-GFP | Z64 | 5 µM | 75.7 ± 6.1 | -24.3 | 0.0135 | 0.081 | Not significant after correction |
| HeLa-GFP | Z64 | 10 µM | 94.9 ± 18.7 | -5.1 | 0.6723 | 0.733 | Not significant |
| HeLa-GFP | TM | 5 µM | 103.6 ± 22.1 | +3.6 | 0.7900 | 0.790 | Not significant |
| HeLa-GFP | TM | 10 µM | 86.0 ± 4.6 | -14.0 | 0.0942 | 0.226 | Not significant |
| HeLa-GFP | Negative control | 5 µM | 114.8 ± 25.7 | +14.8 | 0.3573 | 0.581 | Not significant |
| HeLa-GFP | Negative control | 10 µM | 112.4 ± 31.3 | +12.4 | 0.5068 | 0.676 | Not significant |

Only 293FT-GFP treated with 5 µM Z64 remains below vehicle after correction across all 12 comparisons. This is a technical-well result from one confounded plate and must not be described as biological significance.

### Z64 5 µM versus 10 µM

| Cell line | Difference, 5 minus 10 µM | Unadjusted p | BH q across the two cell-line comparisons | Interpretation |
|---|---:|---:|---:|---|
| 293FT-GFP | -16.4 percentage points | 0.0366 | 0.073 | Not significant after correction |
| HeLa-GFP | -19.2 percentage points | 0.1300 | 0.130 | Not significant |

The fluorescence experiment therefore does not provide corrected evidence that 5 and 10 µM Z64 differ.

## Approximate Western-blot densitometry

The presentation-embedded images were converted to grayscale. Fixed lane and band regions were background corrected. Each PPIL4 or PPIL4-GFP value was normalized to its corresponding GAPDH lane and then expressed relative to vehicle. These estimates are not publication-grade because several bands contain saturated pixels, GAPDH varies between lanes, and raw imager files were unavailable.

| Panel and band | Vehicle | 2.5 µM | 5 µM | 10 µM |
|---|---:|---:|---:|---:|
| HeLa PPIL4-GFP, anti-PPIL4 | 100% | ~67% | ~71% | ~64% |
| HeLa endogenous PPIL4 | 100% | ~81% | ~91% | ~90% |
| 293FT PPIL4-GFP, anti-PPIL4 | 100% | ~82% | ~85% | ~61% |
| 293FT endogenous PPIL4 | 100% | ~80% | ~87% | ~56% |
| 293FT PPIL4-GFP, anti-GFP | 100% | ~75% | ~73% | ~105% |
| SY5Y endogenous PPIL4 | 100% | ~41% | ~22% | ~19% |

No Western-blot p-value is valid because there is one lane per condition and no replicate blot. The anti-GFP and anti-PPIL4 measurements of the 293FT fusion protein also give opposite high-dose patterns, which argues against interpreting the apparent 10 µM anti-GFP rebound as established PPIL4 stabilization.

## Correlation calculations

### Matched fluorescence versus Western blot

The matched observations were HeLa-GFP and 293FT-GFP at 5 and 10 µM Z64, giving four observations. Because two doses are nested within each of only two cell models, these observations are not fully independent.

| Western-blot measurement | Pearson r | Pearson p | Spearman rho | Spearman p | Interpretation |
|---|---:|---:|---:|---:|---|
| PPIL4-GFP detected by anti-PPIL4 | -0.846 | 0.154 | -0.800 | 0.200 | Strong negative direction but not significant; opposite the expected positive assay concordance |
| Endogenous PPIL4 detected by anti-PPIL4 | -0.222 | 0.778 | 0.000 | 1.000 | No detectable relationship |

The present measurements do not demonstrate that lower whole-well GFP fluorescence tracks lower PPIL4-GFP protein abundance.

### Dose correlations

| Measurement | Pearson r with dose | p | Interpretation |
|---|---:|---:|---|
| 293FT-GFP fluorescence, vehicle/5/10 µM | -0.453 | 0.700 | No dose correlation |
| HeLa-GFP fluorescence, vehicle/5/10 µM | -0.199 | 0.872 | No dose correlation |
| HeLa PPIL4-GFP, anti-PPIL4, vehicle/2.5/5/10 µM | -0.753 | 0.247 | Descriptive downward trend only |
| 293FT PPIL4-GFP, anti-PPIL4, vehicle/2.5/5/10 µM | -0.949 | 0.051 | Strong descriptive trend, but one lane per dose |
| 293FT endogenous PPIL4, vehicle/2.5/5/10 µM | -0.928 | 0.072 | Strong descriptive trend, but one lane per dose |
| SY5Y endogenous PPIL4, vehicle/2.5/5/10 µM | -0.817 | 0.183 | Monotonic rank decrease, but one lane per dose |
| 293FT PPIL4-GFP, anti-GFP, vehicle/2.5/5/10 µM | +0.264 | 0.736 | No dose relationship because the 10 µM signal rebounds |

## Conclusion

The only fluorescence comparison surviving correction is 293FT-GFP with 5 µM Z64, and even that comparison represents technical separation within one confounded plate. Western-blot statistical significance cannot be calculated. The cross-assay correlation is negative and non-significant, so fluorescence is not validated as a surrogate for PPIL4-GFP abundance. Independent experiments with unsaturated blot images, balanced plate layouts, and at least three biological replicates are required before biological significance, dose response, or assay correlation can be claimed.
