# Z6466608628 positive-control Western blot interpretation

- Experiment date: 2026-10-06
- Initial interpretation date: 2026-10-07
- Date metadata added: 2026-10-08
- Statistical and correlation follow-up: `reports/2026-10-08_Z64_Fluorescence_Western_Blot_Statistics_and_Correlation.md`

## File identity

- Source blot: `2026-10-06_Z6466608628_Positive_Control_Western_Blot_Original.pptx`
- Compound: Z6466608628 (Z64)
- Doses shown: vehicle/no-drug (`-`), 2.5 µM, 5 µM, and 10 µM
- Cell models shown: HeLa PPIL4-GFP, 293FT PPIL4-GFP, and SY5Y
- Blot labels shown: GFP, PPIL4, and GAPDH
- Treatment duration, replicate number, antibody identities, exposure settings, and raw image format were not provided on the slide.

## Main conclusion

Z64 is a literature-supported PPIL4 positive-control compound, but this blot only partially validates the control locally. The SY5Y panel qualitatively supports depletion of endogenous PPIL4 after Z64 treatment. The HeLa PPIL4-GFP panel does not show convincing dose-dependent loss of PPIL4-GFP, and the 293FT PPIL4-GFP panel is inconclusive because any change is small relative to band saturation and loading-control variation.

The result should therefore be reported as **qualitative evidence that Z64 can reduce endogenous PPIL4 in SY5Y under the tested conditions**, not as quantitative proof that Z64 degrades PPIL4-GFP in HeLa or 293FT cells.

## Panel-by-panel interpretation

| Panel | Bands evaluated | Visual observation | Interpretation |
|---|---|---|---|
| SY5Y | Endogenous PPIL4 near 60-65 kDa; GAPDH near 35 kDa | The PPIL4 band is weaker in the Z64-treated lanes than in the vehicle lane. GAPDH is present but varies across lanes and becomes stronger in the higher-dose lanes. The 10 µM PPIL4 region appears split into two nearby bands. | This is the strongest local positive-control evidence, but the magnitude cannot be estimated reliably without background-corrected densitometry and replicate blots. The complete predefined PPIL4 band region should be quantified consistently at 10 µM. |
| HeLa PPIL4-GFP | PPIL4-GFP near 90-100 kDa; endogenous PPIL4 near 60-65 kDa; GFP; GAPDH | The GFP signal is weak. Neither the putative PPIL4-GFP band nor the endogenous PPIL4 band shows a clear monotonic decrease from 2.5 to 10 µM. GAPDH is comparatively even. | No convincing Z64-dependent degradation is demonstrated in this panel. Weak GFP detection also limits evaluation of the tagged protein. |
| 293FT PPIL4-GFP | PPIL4-GFP near 90-100 kDa; endogenous PPIL4 near 60-65 kDa; GFP; GAPDH | There may be a small reduction at 10 µM, but the main PPIL4 bands are very intense and GAPDH varies between lanes. A clear monotonic dose response is not visible. | Inconclusive. Saturation and uneven loading/transfer can conceal or exaggerate modest changes. |

## Correct control comparison

For each cell model, every treated lane must be compared with the vehicle/no-drug (`-`) lane from the **same cell model and blot panel**. The vehicle lane is not averaged across HeLa, 293FT, and SY5Y.

The normalized value for each lane should be calculated as:

`normalized PPIL4 = background-corrected PPIL4 integrated density / background-corrected GAPDH integrated density`

The vehicle value within each panel should then be set to 1.00 or 100%. For the engineered lines, the approximately 90-100 kDa PPIL4-GFP band and the approximately 60-65 kDa endogenous PPIL4 band should be quantified separately.

Terminology for reporting:

- Z64 is the **positive-control compound**.
- The `-` lane is the **vehicle/no-drug control**.
- The SY5Y Z64 dose series is the panel that most clearly shows the expected positive-control response.
- Z64 should not be described as a validated positive control for the HeLa PPIL4-GFP or 293FT PPIL4-GFP assay until the response is reproduced and quantified in those cell lines.

## Relationship to the fluorescence plate result

The fluorescence experiment showed lower background-corrected whole-well GFP signal with Z64, most clearly at 5 µM, but lacked a robust monotonic dose response. This Western blot does not establish that the fluorescence decrease was caused by PPIL4-GFP degradation:

- HeLa PPIL4-GFP shows no clear reduction of the tagged band.
- 293FT PPIL4-GFP shows, at most, a small and technically uncertain reduction.
- SY5Y supports endogenous PPIL4 depletion but does not directly validate the GFP readout in the engineered cell lines.

The fluorescence result should remain described as a **preliminary Z64-associated decrease in whole-well GFP fluorescence**, while the Western blot supplies qualitative evidence for endogenous PPIL4 depletion only in SY5Y.

## Literature context

Baek et al. identified Z6466608628 as a selective PPIL4 recruiter with a biochemical TR-FRET recruitment EC50 of 0.34 µM. That EC50 is not a cellular degradation potency value. The study reported selective but modest PPIL4 downregulation in MOLT4 cells and incomplete cellular degradation even at 10 µM. Therefore, partial depletion or cell-line-dependent activity is consistent with the published evidence.

Primary source: Baek et al., *Nature Communications* (2025), “Unveiling the hidden interactome of CRBN molecular glues.” https://www.nature.com/articles/s41467-025-62099-w

## Minimum requirements before quantitative claims

1. Obtain the original unsaturated blot images rather than relying only on a slide screenshot.
2. Perform background-corrected densitometry with a predefined band region.
3. Normalize PPIL4 or PPIL4-GFP to GAPDH or, preferably, a validated stable loading control/total-protein measurement.
4. Repeat the complete dose series in at least three independent biological experiments.
5. Report treatment duration, vehicle concentration, antibody identity, exposure settings, and exact biological replicate count.
6. Confirm that the fusion band and endogenous band are identified correctly, ideally using parental cells or a tag-specific control.
7. For mechanism, test whether PPIL4 loss requires CRBN and whether proteasome or cullin-neddylation inhibition rescues the signal.

## Reporting statement

> Z6466608628 treatment was associated with qualitative depletion of endogenous PPIL4 in SY5Y cells. No clear dose-dependent depletion of PPIL4-GFP was observed in HeLa PPIL4-GFP cells, while the 293FT PPIL4-GFP result remained inconclusive because of strong target bands and variable loading-control signal. These data support Z6466608628 as a literature-anchored positive-control compound but do not yet validate it as a reproducible positive control for the engineered GFP assays.
