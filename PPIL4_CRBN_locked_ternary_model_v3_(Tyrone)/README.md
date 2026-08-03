# PPIL4–CRBN Locked Ternary Docking Model v3 (Tyrone)

This self-contained folder is Tyrone's locked 30-candidate PPIL4–CRBN ternary
docking model. See `MODEL_CARD.md` for its scope, contents, and current
validation status.

Prepare the complete candidate queue from the repository root with:

```bash
python3 'PPIL4_CRBN_locked_ternary_model_v3_(Tyrone)/11_run_locked_ppil4_ternary_campaign_(Tyrone).py' --prepare
```

Candidate execution remains intentionally blocked until the FPFT-2216
reference is independently reproduced and `REFERENCE_FPFT.json` is approved as
specified below.

This directory implements `PROTOCOL_LOCK.md`. Do not edit an individual file
without creating a new protocol version and regenerating `MANIFEST.sha256`.
The necessary version-specific placement of the seed is documented in
`IMPLEMENTATION_NOTE.md`.

## Required order

1. Verify every file against `MANIFEST.sha256`.
2. Prepare each ligand with the bundled PRODRG, retaining its renamed atoms.
3. Substitute `<candidate>` in a copy of `guided_ppil4_haddock.cfg`; do not
   change any other setting. `ncores` may differ between machines.
4. Run the pinned `FPFT-2216.cfg` first with `ncores = 1` and write the measured values to
   `REFERENCE_FPFT.json` using `REFERENCE_FPFT.template.json` as the schema.
5. Do not run or rank candidates until the second setup reproduces that reference.
6. Emit candidate results using `results_schema.csv` and the exact contents of
   `PROTOCOL_VERSION.txt` as `protocol_version`.

## Locked candidate universe

Candidate runs are restricted to the 30 rows in
`candidates_30_final_corrected_ranked.csv`. They are the 30 non-control rows
from the pinned `sources/final_corrected_31_ranked.csv` file in its existing
corrected HADDOCK order. `locked_candidate_rank` records the contiguous queue
order, while `final_corrected_haddock_rank` preserves the source rank.

Do not re-sort, substitute, or extend this list within this protocol version.
The source file's 31st row type, FPFT-2216, is already the pinned positive
control and is deliberately excluded from the 30-candidate queue. The three
negative controls are also controls, not candidate-ranking members.

## Locked result calculation

Use `RANKING_RULES.md` for the exact all-model definition of `mean_best_10`,
tie handling, insufficient-model behavior, and final candidate ordering.
Use `aggregate_minimized_models.py` on the final post-`emref`
`caprieval/capri_ss.tsv`; it implements that definition and emits the locked
result schema.

## Fixed pocket derivation

`CRBN_receptor_thalidomide_Ryan.pdb` is receptor-only despite its filename. The
CRBN active residues were therefore derived once as chain-A protein residues
with any atom within 5.0 Å of EF2 residue 501 in the pinned bound complex
`inputs/CRBN-Thalidomide-SALL4_Ryan_bound_reference.pdb`. The receptor-only PDB
remains the docking input. PPIL4 active residues are locked to 249, 250, 273,
275, 276, 277, 278, and 279. Passive residues were computed once with
HADDOCK3 2026.7.0 `passive_from_active`.

`ambig_FIXED.tbl` is generated once from the two fixed active/passive files and
the fixed ligand residue 900. It must be reused byte-for-byte for every ligand.

## Locked DockQ reference

Every CAPRI stage uses
`inputs/reference_9DWV_CRBN_A_PPIL4_B_FPFT_C.pdb`. It is the deposited 9DWV
ternary structure remapped to the locked chain convention. DockQ is calculated
for CRBN chain A against PPIL4 chain B. Candidate ligand chain C is excluded
from the DockQ atom matching, so DockQ measures recovery of the protein-protein
ternary geometry rather than chemical identity with FPFT-2216.
