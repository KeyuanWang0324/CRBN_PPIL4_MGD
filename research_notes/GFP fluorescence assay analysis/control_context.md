# Biological and experimental context for Z6466608628, TM-04-064-02, and the GFP fluorescence assay

## What is established for Z6466608628 (Z64)?

### Takeaway
Z6466608628 is the only compound in this experiment for which the cited primary paper directly demonstrates PPIL4 recruitment, in-vitro ubiquitylation, and cellular PPIL4 loss. Its reported EC50 of 0.34 µM is a **biochemical TR-FRET recruitment EC50**, not a cellular PPIL4 degradation DC50, viability IC50, or efficacy value. Cellular degradation was modest and incomplete, and it was established in MOLT4 cells rather than GFP-engineered 293FT or HeLa cells.

### Cited Findings
- Baek et al. screened more than 6,000 IMiD analogues by TR-FRET using GFP-CRBN-DDB1ΔB and biotinylated PPIL4. Z6466608628 produced a recruitment EC50 of 0.34 µM; the paper explicitly describes this assay as measurement of compound-induced PPIL4 recruitment to CRBN. This is therefore a biochemical ternary-complex value, not a cellular degradation or viability value. — [Baek et al., Nature Communications (2025)](https://www.nature.com/articles/s41467-025-62099-w)
- The TR-FRET assay used 20 or 50 nM biotinylated neosubstrate, 2 nM terbium-labelled streptavidin, 200 nM eGFP-CRBN-DDB1ΔB, a 1-hour room-temperature incubation, and the 520/490-nm ratio. The eGFP in this purified-protein assay labels CRBN and is not equivalent to whole-well GFP fluorescence from engineered living cells. — [Baek et al., Methods](https://www.nature.com/articles/s41467-025-62099-w)
- In Kelly-cell lysate IP-MS, Z6466608628 selectively recruited PPIL4 together with PPIL4's known partner DHX40. In MOLT4 whole-cell proteomics, 1 µM Z6466608628 for 5 hours caused modest PPIL4 downregulation; the experiment used n=2 compound samples and n=4 DMSO controls. — [Baek et al., Fig. 5 and Methods](https://www.nature.com/articles/s41467-025-62099-w)
- The official Supplementary Data 5 reports PPIL4 log2 fold-change = -0.5289921 in the MOLT4 whole-cell experiment (approximately 0.69-fold abundance, or a 31% decrease), nominal P=6.61×10^-5 and adjusted P=0.0146. In the Kelly lysate IP-MS experiment, PPIL4 had log2 enrichment = 4.066273 (approximately 16.7-fold), nominal P=6.44×10^-5 and adjusted P=0.0196. These are different assays and must not be conflated. — [Official Supplementary Data 5](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-62099-w/MediaObjects/41467_2025_62099_MOESM7_ESM.xlsx)
- The authors confirmed dose-dependent PPIL4 depletion by western blot over 0.1, 1, 5, and 10 µM Z6466608628 and demonstrated time-dependent in-vitro PPIL4 ubiquitylation. They nevertheless state that cellular PPIL4 degradation remained incomplete even at 10 µM. — [Baek et al., Supplementary Fig. 8 and main text](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-62099-w/MediaObjects/41467_2025_62099_MOESM1_ESM.pdf)
- The published experimental ternary structure is not a Z6466608628 structure. PDB 9DWV is a 3.50-Å cryo-EM structure of CRBN-DDB1, the PPIL4 RRM domain, and **FPFT-2216**. The Z6466608628 pose in Fig. 5g is a docking model constrained by experimental SAR. — [RCSB PDB 9DWV](https://www.rcsb.org/structure/9DWV); [Baek et al., Fig. 5](https://www.nature.com/articles/s41467-025-62099-w)

### Inferences
- Z6466608628 is a defensible **literature-anchored positive-control candidate** for PPIL4 recruitment/depletion, but only after its identity and activity are verified in the actual stock and in each local cell model.
- A 5 or 10 µM treatment is above the biochemical recruitment EC50, but this does not predict the magnitude of cellular degradation; permeability, CRBN/PPIL4 abundance, PPIL4 complex occupancy, time, and degradation machinery can all change the response.
- The observed non-monotonic whole-well GFP response (stronger reduction at 5 µM than 10 µM) cannot be interpreted as a PPIL4 degradation dose response. The published PPIL4 western blot was monotonic but incomplete, whereas the local readout is generic fluorescence rather than PPIL4 protein.

### Gaps
- The paper and supplement do not clearly state the cell line and exposure time for the Supplementary Fig. 8a western blot; this should not be guessed.
- No primary evidence was located establishing Z6466608628-driven PPIL4 degradation in 293FT-GFP or HeLa-GFP cells.
- No cellular PPIL4 DC50 or Dmax is reported. The 0.34 µM number must not be labelled DC50, IC50, viability potency, or PPIL4 degradation potency.

## What is established for TM-04-064-02 (TM)?

### Takeaway
No peer-reviewed primary source located in this search identifies a compound by the exact name **TM-04-064-02** or provides a traceable PPIL4 recruitment/degradation experiment under that name. Commercial pages cite Baek et al., but the paper uses public drug names and internal compound codes rather than TM-04-064-02. Therefore TM cannot presently be treated as a validated PPIL4 positive control.

### Cited Findings
- Baek et al. profiled 20 CRBN-directed compounds at 1 µM in MOLT4 and Kelly cell lysates using activity-impaired recombinant CRBN-DDB1ΔB and IP-MS. This assay identifies compound-dependent recruitment and, by design, prevents CUL4 engagement and productive ubiquitylation; recruitment in this screen is not itself cellular degradation. — [Baek et al., Fig. 2 and Methods](https://www.nature.com/articles/s41467-025-62099-w)
- The paper emphasizes that recruitment and degradation are distinct: some proteins were recruited to CRBN without productive degradation, and the workflow was intended to find interactions missed by endpoint degradation screens. — [Baek et al., Discussion](https://www.nature.com/articles/s41467-025-62099-w)
- Commercial listings describe TM-04-064-02 as inducing CRBN association of IKZF1/3, PPIL4, ETF1, and PDE6D and cite Baek et al.; however, this exact compound name is absent from the article text and official Supplementary Data 4/5 examined here. These listings are secondary/vendor assertions, not independent primary validation. — [MedChemExpress product page](https://www.medchemexpress.com/tm-04-064-02.html); [Baek et al. primary paper](https://www.nature.com/articles/s41467-025-62099-w)

### Inferences
- TM-04-064-02 may be a vendor-assigned name for one of the paper's internal compounds, but that mapping cannot be assumed without a structure, catalog/lot documentation, or an author/vendor cross-reference.
- Even if the identity is mapped to a compound that enriched PPIL4 in lysate IP-MS, that result would establish recruitment under an artificial lysate assay, not PPIL4 degradation in living 293FT-GFP or HeLa-GFP cells.
- TM should be labelled **test compound / unvalidated comparator**, not “PPIL4 positive control,” until direct PPIL4 protein loss and CRBN dependence are shown locally.

### Gaps
- The exact chemical identity, CAS number, vendor catalog/lot, purity, and correspondence between “TM-04-064-02” and the internal compound identifiers in Baek et al. are unresolved.
- No traceable primary cellular PPIL4 degradation value, dose-response, time course, or assay system was found for the exact TM-04-064-02 name.
- The “neg” compound identity is also necessary before concluding that it is mechanistically inert; a no-drug/DMSO vehicle control is the only clearly interpretable baseline without that information.

## Can published activity be transferred to GFP-engineered 293FT and HeLa cells?

### Takeaway
No. Published Z6466608628 evidence supports a mechanism in purified proteins, Kelly lysate, and MOLT4 cells. It does not validate a generic GFP fluorescence endpoint in 293FT-GFP or HeLa-GFP cells as a PPIL4-degradation readout. Unless GFP is fused directly and functionally to endogenous PPIL4, whole-well GFP mainly reflects reporter abundance plus cell number/state and optical interference.

### Cited Findings
- The published whole-cell Z6466608628 experiment was performed in MOLT4 cells at 1 µM for 5 hours and directly quantified the proteome; the local experiment instead uses GFP-engineered adherent 293FT and HeLa cells and a whole-well fluorescence endpoint. — [Baek et al., Methods](https://www.nature.com/articles/s41467-025-62099-w)
- A primary GFP viability study found that GFP fluorescence tracked total protein more closely than all viability measures and that GFP behaved differently under different toxic mechanisms; the authors concluded that a multiparametric assay is necessary for precise viability interpretation. — [Juhász et al., *Sensors* (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6100089/)
- Test compounds can interfere with GFP-channel measurements through autofluorescence. A primary microplate study used fluorescence polarization specifically to discriminate GFP from compound autofluorescence. — [Knight et al., *Mutagenesis* (2002), PubMed](https://pubmed.ncbi.nlm.nih.gov/12062116/)
- Baek et al. themselves rejected one apparent TR-FRET screening hit after determining it was autofluorescent, demonstrating that optical artifacts are relevant even in the discovery assay from which Z6466608628 emerged. — [Baek et al., Fig. 5 results](https://www.nature.com/articles/s41467-025-62099-w)

### Inferences
- The processed table can support only a statement such as “treatment-associated change in background-corrected whole-well GFP fluorescence relative to the cell-line-matched no-drug mean.” It cannot establish PPIL4 degradation, CRBN dependence, molecular-glue activity, cytotoxicity, or cell proliferation by itself.
- The decrease seen with 5 µM Z64 in both cell lines is hypothesis-generating. Its failure to strengthen at 10 µM, the increase produced by the negative compound in some conditions, and larger HeLa dispersion favor cautious interpretation and artifact/cell-state counterchecks.
- Because the four images are repeated reads of the same plate rather than independent treatments, they improve measurement precision but do not provide four biological replicates. Statistical inference should use independently cultured/treatment plates as the experimental unit.

### Gaps
- It is not documented whether GFP is constitutive, promoter-driven, or fused to PPIL4; this distinction is essential.
- Exposure duration, GFP excitation/emission settings, gain, plate-reader linearity, confluence at treatment, compound-only fluorescence, and non-GFP parental-cell fluorescence are needed to interpret the endpoint.
- The number of independent cultures/plates is not established; within-plate wells alone do not demonstrate biological reproducibility.

## What orthogonal evidence is required to claim PPIL4 degradation and CRBN dependence?

### Takeaway
At minimum, the claim requires direct endogenous PPIL4 quantification, proof that protein loss depends on CRBN/CRL4 and the proteasome, and separation of target loss from cell-number, transcriptional, and optical effects. Recruitment, ubiquitination, and selectivity assays make the mechanism substantially stronger.

### Cited Findings
- Baek et al. used a layered chain of evidence for Z6466608628: purified-protein TR-FRET recruitment, lysate IP-MS recruitment/selectivity, whole-cell proteomics, PPIL4 western blot dose response, and in-vitro ubiquitylation. The paper also demonstrates that recruitment alone need not yield degradation. — [Baek et al.](https://www.nature.com/articles/s41467-025-62099-w)
- A recent primary CRBN molecular-glue study established mechanism using target immunoblotting, CRBN-WT versus CRBN-knockout cells, rescue experiments, proteasome inhibition (MG132), CRL inhibition (MLN4924), ternary-complex assays, ubiquitination, and proteomic selectivity. This illustrates the controls needed to separate protein loss from a nonspecific phenotype. — [Molecular glue degraders of HuR, *Nature* (2026)](https://www.nature.com/articles/s41586-026-10613-5)

### Inferences
- **Direct degradation endpoint:** quantify endogenous PPIL4 by western blot with a validated antibody, normalized to total-protein stain or stable loading controls, across a dose series and time course. Use at least three independent biological experiments. Targeted MS or global quantitative proteomics is a strong orthogonal confirmation and reveals off-target degradation.
- **CRBN dependence:** repeat the PPIL4 endpoint in CRBN-knockout/knockdown cells and restore activity by wild-type CRBN rescue. A CRBN-binding-defective compound analogue or CRBN-competition control would further support direct CRBN engagement.
- **UPS mechanism:** demonstrate rescue of PPIL4 abundance by proteasome inhibition and by inhibition of cullin neddylation/CRL activity, using inhibitor-only and toxicity controls. Show PPIL4 ubiquitination under compound treatment where feasible.
- **Recruitment:** demonstrate a compound-dependent PPIL4–CRBN ternary complex using TR-FRET, NanoBRET, co-immunoprecipitation, or another direct proximity assay. Recruitment should not be presented as degradation.
- **Exclude transcriptional effects:** measure PPIL4 mRNA by RT-qPCR; stable mRNA with falling PPIL4 protein supports post-translational loss.
- **Separate GFP/cell-state effects:** pair GFP with direct cell counts or nuclei counts, an orthogonal viability assay, and microscopy/flow cytometry. Include compound-only wells at the GFP wavelengths, parental non-GFP cells, vehicle-matched controls, and equal final DMSO. Bubbles and out-of-range detector values must be excluded by predefined QC rules.
- **Positive-control language:** Z6466608628 can be called the literature-anchored PPIL4 positive-control candidate. It becomes an assay positive control only after reproducibly decreasing endogenous PPIL4 in that assay and cell line. TM-04-064-02 should remain an unvalidated comparator until its identity and mechanism are confirmed.

### Gaps
- Definitive inhibitor concentrations and pretreatment times should be optimized for these exact cells rather than copied across cell lines; inhibitor toxicity can itself distort GFP and protein abundance.
- A validated CRBN-knockout/rescue 293FT-GFP and HeLa-GFP pair is not currently documented.
- Without endogenous PPIL4 measurements, the current fluorescence table cannot answer whether either compound degraded PPIL4.
