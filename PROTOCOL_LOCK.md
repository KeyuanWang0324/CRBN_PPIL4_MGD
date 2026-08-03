# PROTOCOL_LOCK — the one docking protocol you both run

Ryan and Tyrone: this file exists so your two rankings become **comparable**. It locks
**one** protocol. Commit the executable artifacts (cfg, restraints, inputs, seed, versions)
— not a prose description — and both run *those files*.

> **Implementation status:** wired into `08_haddock3_ternary_with_ligand_(Ryan).py`
> (`ppil4_lock_v1`). Two errors in the original draft were found and corrected while wiring
> it — see **§8**. The settings below are the corrected, as-implemented ones.

## Rule 0 — share the artifact, not a description
A HADDOCK score is only comparable **within one exact protocol**. Share and pin the
**executable artifacts** below, commit them, and both run them. Your **wrapper/launcher may
differ** (that's even useful as a cross-check); the **cfg / restraints / inputs / seed must not.**

---

## 1. The locked config
Ryan's `08` protocol with two changes marked `### LOCK`. `08` generates a per-candidate cfg
from this exact template (only the `<candidate>` ligand paths differ):

```toml
run_dir = "run1"
ncores = 8                     # speed only — does NOT change scores (but see §5: ncores=1 for byte-exact)

molecules = [
    "CRBN_receptor_thalidomide_Ryan.pdb",        # chain A
    "<candidate>_prodrg_fixed.pdb",              # chain C, resnum 900 (PRODRG, per-candidate)
    "PPIL4_chainB.pdb"                           # chain B
]

[topoaa]
ligand_top_fname   = "<candidate>_prodrg.top"    # PRODRG explicit — NOT autotoppar
ligand_param_fname = "<candidate>_prodrg.param"

[rigidbody]
ambig_fname        = "ambig_FIXED.tbl"           ### LOCK #2: one fixed restraint file for ALL candidates
ligand_top_fname   = "<candidate>_prodrg.top"
ligand_param_fname = "<candidate>_prodrg.param"
sampling           = 1000
iniseed            = 42                           ### LOCK #1: pin the CNS seed (param is 'iniseed', see §8)

[caprieval]
[seletop]
select = 200
[flexref]
ambig_fname        = "ambig_FIXED.tbl"
ligand_top_fname   = "<candidate>_prodrg.top"
ligand_param_fname = "<candidate>_prodrg.param"
[caprieval]
[emref]
ambig_fname        = "ambig_FIXED.tbl"
ligand_top_fname   = "<candidate>_prodrg.top"
ligand_param_fname = "<candidate>_prodrg.param"
[caprieval]
[clustfcc]
[caprieval]
```

**Ligand topology = PRODRG, explicit** on every module. Do **not** use `autotoppar = true`
(it renames ligand atoms and desyncs them from the coordinates — documented in 08).

---

## 2. The fixed restraints (`ambig_FIXED.tbl`) — LOCK #2, the important one
Identical for every candidate, so the only thing that changes between runs is the molecule.
Previously 08 derived the CRBN-side active residues from *each candidate's own Vina pose*,
leaking Vina noise into the ranking. Now fixed:

- **PPIL4 side** (`ppil4_active_FIXED.txt`): chain B active `249, 250, 273, 275, 276, 277, 278, 279`
  (real RRM interface from 9DWV / FPFT-2216).
- **CRBN side** (`crbn_active_FIXED.txt`): the glutarimide/degron pocket —
  `350, 351, 352, 353, 357, 377, 378, 379, 380, 381, 386, 400, 402`. Derived once as chain A
  residues within 5 Å of the co-crystal thalidomide (see §8 for the file-source correction).
- `ambig_FIXED.tbl` is built once from these against the ligand (resnum 900) and committed.
  Because the ligand is always resnum 900, this is genuinely **one file for all candidates and controls.**

---

## 3. Fixed inputs + versions
| Artifact | Pin |
|---|---|
| CRBN receptor | `CRBN_receptor_thalidomide_Ryan.pdb` (chain A) |
| PPIL4 | `PPIL4_chainB.pdb` (chain B, committed) |
| Ligand starting pose | **Ryan's `04` recipe** — Vina into the CRBN thalidomide pocket — same for candidates AND controls |
| HADDOCK3 | **2026.7.0** (`haddock3 --version`) |
| PRODRG | bundled `haddock/prodrg/prodrg_<arch>` |

---

## 4. What may differ / must not
| May differ | Must be identical |
|---|---|
| launcher / parser / plotting | `ambig_FIXED.tbl`, `*_active_FIXED.txt` |
| which AI wrote the wrapper | cfg modules/params, `iniseed`, input PDBs, ligand-pose recipe, HADDOCK3 version, PRODRG |

---

## 5. Acceptance test (golden reference) — do BEFORE any candidate
1. Run **FPFT-2216** through the locked cfg. Record `mean_best_10, best_haddock, sd_best_10,
   cluster_id, n_models, dockq_vs_9dwv` → commit `REFERENCE_FPFT.json`.
2. The other runs the **same FPFT-2216** first and must reproduce it:
   - `iniseed` pinned + same version + **`ncores = 1` → exact**;
   - multi-core → allow tolerance, but **same cluster** and DockQ within ~0.05.
3. **Only if FPFT reproduces** do you trust the run and move to candidates.

---

## 6. Output schema (one CSV, fixed columns)
```
molecule, mean_best_10, best_haddock, sd_best_10, dockq_vs_9dwv, cluster_id, n_models, protocol_version
```
`mean_best_10/best_haddock/sd_best_10` = over the 10 best single models (`capri_ss.tsv`).
`dockq_vs_9dwv` = from the `09` validation. `protocol_version` = `ppil4_lock_v1`.

---

## 7. Divergent wrappers as a cross-check
Both point at the same locked cfg + shared control set (FPFT + 2–3 common molecules).
Numbers match → unified *and* independently reproducible (strong result). Numbers differ →
a real bug in inputs/version/seed, caught on 3 molecules instead of after a full run.

---

## 8. Implementation corrections (found while wiring `ppil4_lock_v1`)
1. **Seed param is `iniseed`, not `seed`.** HADDOCK3 2026.7.0 has no top-level global `seed`;
   the sampling seed is `iniseed` (default 917) on `[rigidbody]`. A top-level `seed = 42` is
   silently ignored. Also: even with `iniseed` pinned, multi-core runs are **not** byte-identical
   (CNS orders parallel models nondeterministically) — exact reproduction needs `ncores = 1`.
2. **The degron pocket came from the co-crystal, not the receptor.** `CRBN_receptor_thalidomide_Ryan.pdb`
   is apo (Zn only, no thalidomide), so "within 5 Å of the bound thalidomide" in that file is empty.
   The fixed CRBN residues were derived from the co-crystal `CRBN-Thalidomide-SALL4_(Ryan).pdb`
   (ligand EF2) and confirmed present in the apo receptor with matching numbering.
