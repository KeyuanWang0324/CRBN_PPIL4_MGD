# PROTOCOL_LOCK — the one docking protocol you both run

Ryan and Tyrone: this file exists so your two rankings become **comparable**. It locks
**one** protocol. Commit the executable artifacts (cfg, restraints, inputs, seed, versions)
— not a prose description — and both run *those files*.

> **Implementation status:** wired into `08_haddock3_ternary_with_ligand_(Ryan).py`
> (`ppil4_lock_v1`) and scored by `09_score_vs_fpft2216_reference_(Ryan).py`. All 18
> candidates and all 6 controls have been run under it; the golden reference
> (`REFERENCE_FPFT.json`, §5) is emitted and committed. Errors found while wiring and
> running it are in **§9**. The settings below are the corrected, as-implemented ones.

## Rule 0 — share the artifact, not a description
A HADDOCK score is only comparable **within one exact protocol**. Share and pin the
**executable artifacts** below, commit them, and both run them. Your **wrapper/launcher may
differ** (that's even useful as a cross-check); the **cfg / restraints / inputs / seed must not.**

---

## 1. The locked config
Ryan's `08` protocol with two changes marked `### LOCK`. `08` generates a per-candidate cfg
from this exact template (only the `<candidate>` ligand paths and `run_dir` differ):

```toml
run_dir = "<run_root>/haddock3_ternary_with_ligand_run/<candidate>/run1"   # see §3
ncores = <cpu_count - 1>       # speed only — does NOT change scores (but see §5: ncores=1 for byte-exact)

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
iniseed            = 42                           ### LOCK #1: pin the CNS seed (param is 'iniseed', see §9)

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

The final module is `9_caprieval`, and every downstream consumer (07, `00_validate_..._08`,
09, and 08's own resume check) reads `9_caprieval/capri_ss.tsv` by that name. Don't add or
remove modules — that renumbers the step directories and silently breaks all four.

---

## 2. The fixed restraints (`ambig_FIXED.tbl`) — LOCK #2, the important one
Identical for every candidate and control, so the only thing that changes between runs is
the molecule. Previously 08 derived the CRBN-side active residues from *each candidate's own
Vina pose*, leaking Vina noise into the ranking. Now fixed:

- **PPIL4 side** (`ppil4_active_FIXED.txt`): chain B active `249, 250, 273, 275, 276, 277, 278, 279`
  (real RRM interface from 9DWV / FPFT-2216); line 2 of that file is the passive shell.
- **CRBN side** (`crbn_active_FIXED.txt`): the glutarimide/degron pocket —
  `350, 351, 352, 353, 357, 377, 378, 379, 380, 381, 386, 400, 402`. Derived once as chain A
  residues within 5 Å of the co-crystal thalidomide (see §9 for the file-source correction);
  line 2 is the passive shell.
- `ambig_FIXED.tbl` is built once from these against the ligand (resnum 900) and committed:
  **23 `assign` statements** — 13 CRBN(A)→ligand, 8 PPIL4(B)→ligand, plus the two reciprocal
  ligand→set restraints. Because the ligand is always resnum 900, this is genuinely **one
  file for all candidates and controls.**
- **There is no direct CRBN↔PPIL4 restraint.** The ternary is bridged *only* through the
  ligand — which is the point: the protein–protein interface has to be earned by the
  molecule, not imposed by the restraints.

08 refuses to start if `ambig_FIXED.tbl`, `PPIL4_chainB.pdb`, or
`CRBN_receptor_thalidomide_Ryan.pdb` is missing.

---

## 3. Fixed inputs + versions
| Artifact | Pin |
|---|---|
| CRBN receptor | `CRBN_receptor_thalidomide_Ryan.pdb` (chain A) |
| PPIL4 | `PPIL4_chainB.pdb` (chain B, committed) |
| Ligand starting pose | **Ryan's `04` recipe** — Vina into the CRBN thalidomide pocket — same for candidates AND controls |
| HADDOCK3 | **2026.7.0** (`haddock3 --version`) |
| PRODRG | bundled `haddock/prodrg/prodrg_<arch>` |
| 9DWV reference | `reference_structures/FPFT-2216_9DWV_reference_(Ryan).pdb` — chains B(CRBN)→A, C(PPIL4)→B; DDB1, ligand and Zn dropped |

### 3a. Where the run tree lives — never inside a synced folder
`run_dir` is **not** in the project. HADDOCK3 writes thousands of small files per candidate
and hands off between steps through a large `io.json`; a file-sync daemon (iCloud on
`~/Desktop` here) races it and swaps `io.json` mid-read, so `caprieval` sees an empty
topology and dies with `find_ff -> models[0].topology[0]` IndexError. That killed **all 18
candidates** in one overnight run, while the earlier controls — run while sync was idle —
survived.

The run root therefore defaults to `~/haddock_runs/` (override with `$HADDOCK_RUN_ROOT`) and
must be outside any synced directory. Only the disposable run tree moves out; the
deliverables (§6) stay in the project. This changes no number, but ignoring it costs the
whole batch.

---

## 4. What may differ / must not
| May differ | Must be identical |
|---|---|
| launcher / parser / plotting | `ambig_FIXED.tbl`, `*_active_FIXED.txt` |
| `ncores`, run root, machine | cfg modules/params, `iniseed`, input PDBs, ligand-pose recipe, HADDOCK3 version, PRODRG |
| which AI wrote the wrapper | the module sequence and its step numbering (final = `9_caprieval`) |

---

## 5. Acceptance test (golden reference) — do BEFORE any candidate
1. Run **FPFT-2216** through the locked cfg, then run `09` — it writes `REFERENCE_FPFT.json`
   with the locked-schema stats plus the full CAPRI metric set. **Committed, current values:**

   | | |
   |---|---|
   | `mean_best_10` | **-111.382** |
   | `best_haddock` | **-130.503** |
   | `sd_best_10` | 7.296 |
   | `cluster_id` / `n_models` | 9 / 161 |
   | `dockq_vs_9dwv` | **0.376** (CAPRI *acceptable*, ≥ 0.23) |
   | `irmsd` / `fnat` / `lrmsd` / `ilrmsd` | 2.073 / 0.647 / 21.31 / 4.245 |

2. The other runs the **same FPFT-2216** first and must reproduce it:
   - `iniseed` pinned + same version + **`ncores = 1` → exact**;
   - multi-core → allow tolerance, but **same cluster** and DockQ within ~0.05.
3. **Only if FPFT reproduces** do you trust the run and move to candidates.

09 prints PASS/FAIL against the 0.23 CAPRI bar automatically. **This check currently passes**:
handed the real drug, the protocol recovers roughly the real CRBN–PPIL4 architecture. That
validates the *methodology* — it says nothing about any individual candidate.

---

## 6. Output schema (two CSVs, same fixed columns)
```
molecule, mean_best_10, best_haddock, sd_best_10, dockq_vs_9dwv, cluster_id, n_models, protocol_version
```
- `08_final_ternary_with_ligand_results_(Ryan).csv` — candidates, best `TOP_N_KEEP` by `mean_best_10`.
- `08_controls_results_(Ryan).csv` — the calibration set (§7), kept out of the candidate ranking.

`mean_best_10/best_haddock/sd_best_10` = over the 10 best single models
(`9_caprieval/capri_ss.tsv`). `cluster_id` = the best model's cluster. `protocol_version` =
`ppil4_lock_v1`, stamped on every row so old-protocol rows can't quietly mix in.

**Run order: 08, then 09.** 08 always writes `dockq_vs_9dwv` blank (its own caprieval dockq is
self-referential — the reference there is the run's own best model); 09 fills it for every
molecule in *both* CSVs by re-scoring 08's existing rank-1 models against 9DWV (~1 s each, no
re-docking). One asymmetry to know: the controls CSV *is* its own ledger, so its dockq values
survive an 08 re-run — the candidates' ledger lives in the run tree with the column blank, so
**re-running 08 blanks `dockq_vs_9dwv` in the finalists CSV. Re-run 09 after any 08 run.**

### `dockq_vs_9dwv` is a diagnostic, NOT a ranking column
The column is now filled for all 24 molecules, and the sweep that filled it **falsified its
use as a descriptor.** Measured under `ppil4_lock_v1`:

- FPFT-2216 tops all 24 at **0.376** — the §5 methodology check passes, cleanly.
- But `negctrl_no_ppil4_arm` (0.294) and `negctrl_no_crbn_ppil4_arm` (0.290) rank **2nd and
  3rd, above every candidate**. The latter is the fragment with both functional handles
  stripped and the worst HADDOCK score of anything docked (-92.19).
- `Z6466608628_positive_control`, the strongest positive control by HADDOCK score, scores 0.061.
- Spearman vs `mean_best_10` is +0.055, vs ligand heavy-atom count +0.066 — so it is neither
  a restatement of the score nor a size artifact.
- The distribution is bimodal (irmsd ~2–3 Å vs ~5–9 Å): it mostly records whether *this one
  rank-1 model out of ~160* landed in the native-like basin. A sampling outcome, not a
  property of the molecule.

It is also **not independent** in the first place: the PPIL4 restraints are themselves derived
from 9DWV's contact residues, so PPIL4 lands near the real interface almost by construction.

**Therefore: rank on `mean_best_10`. Never prefer or reject a molecule on `dockq_vs_9dwv`.**
Use it only to spot which models are worth opening in PyMOL. 09's `main()` carries a
calibration guard that recomputes best-negative-vs-candidates every run and prints this
warning automatically, so it can't be quietly forgotten. The lock-era buy list's structure
axis is built from BSA / cluster convergence / AIR / `sd_best_10` instead — none of which
touch this column.

---

## 7. The control set (calibration) — run these, they set the bar
Docked through the **same** locked protocol, written to their own CSV:

| Control | Role | `mean_best_10` | `dockq_vs_9dwv` |
|---|---|---|---|
| `FPFT_2216_positive_control` | real drug; golden reference | -111.382 | 0.376 |
| `Z6466608628_positive_control` | isoindolinone-piperazine (Enamine) | **-125.610** | 0.061 |
| `TM-04-064-02_positive_control` | adenine-linked thalidomide | -109.961 | 0.143 |
| `negctrl_no_glutarimide` | degron-dead | **-112.202** | 0.029 |
| `negctrl_no_crbn_ppil4_arm` | both handles stripped | -92.187 | 0.290 |
| `negctrl_no_ppil4_arm` | no PPIL4 arm | -105.097 | 0.294 |

Two numbers here do real work:

- **The noise floor is ~-112.** `negctrl_no_glutarimide` (-112.202) *outscores* FPFT-2216
  (-111.382), so anything above ~-112 is not distinguishable from a degron-dead negative and
  "beats FPFT-2216" is a meaningless bar.
- **The usable bar is Z6466608628 (-125.61)**, the strongest control — that is what a
  candidate has to beat to mean anything.

A control whose 04 Vina pose doesn't exist yet is skipped with a message (08 can't make poses
— no Vina in the haddock venv), so run 04 with `MANUAL_CANDIDATE` for it first. Both wrappers
must run the full control set: it is the only calibration this protocol has.

**FPFT-2216 is docked but never bought.** It is not commercially available — sourcing it would
require custom synthesis — so it is deliberately absent from the buy list and the purchase
sheet, which carry the other five controls. That is not an oversight to "fix": it is the
golden reference (§5) and must keep being docked and scored. The same caution applies in
reverse — appearing on the buy list does not imply a molecule is orderable; the `RW_*` rows
are novel analogs generated by 01.

---

## 8. Divergent wrappers as a cross-check
Both point at the same locked cfg + shared control set (FPFT + 2–3 common molecules).
Numbers match → unified *and* independently reproducible (strong result). Numbers differ →
a real bug in inputs/version/seed, caught on 3 molecules instead of after a full run.

### 8a. The shared molecules — implemented
Three of Tyrone's guided-screen picks are now docked through **Ryan's** wrapper under this
lock, so both screens have numbers for the same molecules. They run through the same
04 → 08 → 09 path as everything else and are named `tyrone_<his number>`:

| name | Tyrone's | SMILES | also in 01's library as |
|---|---|---|---|
| `tyrone_118` | Candidate 118 | `COc1ccc2c(c1C#N)C(=O)N(C1CCC(=O)NC1=O)C2=O` | `cand_1981` |
| `tyrone_45` | Candidate 45 | `Cc1ccc2c(c1Br)C(=O)N(C1CCC(=O)NC1=O)C2=O` | `cand_467` |
| `tyrone_6` | Candidate 6 | `CCc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O` | `cand_85` |

All three are **already in Ryan's own 01-generated library** under a `cand_*` name (identical
SMILES) — they simply never made 03/06's cut, so none was ever docked in 08. There is no
double-counting in the finalists. They are re-docked under the `tyrone_*` names rather than
reusing the twins' 04 poses: Vina's search is stochastic and the twins' on-disk poses no
longer match the affinities recorded for them in 04's summary (that file predates its RESUME
flag). Re-docking gives one self-consistent pose → affinity → 08 chain.

Each stage keeps them in its **own** CSV, never the funnel's ranked output — they are not from
03/06's pool, and folding them in would silently change what "top 20 of the funnel" means:

| stage | list | output |
|---|---|---|
| 04 | `CROSSCHECK` | `04_crosscheck_vina_scores_(Ryan).csv` |
| 05 (legacy) | `CROSSCHECK` + `RUN_CROSSCHECK_ONLY` | `05_crosscheck_ternary_scores_(Ryan).csv` |
| 08 | `CROSSCHECK` | `08_crosscheck_results_(Ryan).csv` (locked schema, §6) |
| 09 | — | fills `dockq_vs_9dwv` in that CSV too |
| 11 | — | `final_buy_list_crosscheck_(Ryan).csv` (buy-list schema) |

The buy list renames these to `TZ_<his number>` (and Ryan's to `RW_<n>`) for readability, and
keeps the pipeline name in a `run_name` column. **Everything in this document, and every run
directory and CSV outside the buy list, uses the pipeline names** — `tyrone_118`, not `TZ_118`.

**Result — all three land below the noise floor:**

| molecule | `mean_best_10` | `dockq_vs_9dwv` | vs the bars |
|---|---|---|---|
| `tyrone_6` | -108.564 | 0.188 | above the ~-112 noise floor |
| `tyrone_118` | -104.011 | 0.281 | above the ~-112 noise floor |
| `tyrone_45` | -98.973 | 0.062 | above the ~-112 noise floor |

All three score **worse than every one of Ryan's 18** and worse than `negctrl_no_glutarimide`
(-112.202), so under this protocol they are not distinguishable from a degron-dead negative.
Read that as a **statement about the score, not a verdict on the molecules**: all three are
bare phthalimides with no PPIL4-facing arm, and §6/README already record that the HADDOCK
score is size- and interface-dominated. A molecule with no arm to bury cannot produce a large
interface, so this is the expected direction. What it does establish is that Tyrone's guided
ranking and Ryan's `mean_best_10` **do not agree** on these three — which is exactly the
disagreement a cross-check exists to surface, and it is not explained by inputs/seed/version
(the FPFT golden reference reproduced byte-identically in the same session).

`tyrone_118` and `tyrone_6` also give the cleanest available demonstration of why 05 was
retired: run through 05 they return **byte-identical** scores (-144.986 / dockq 0.277 /
fnat 0.314) — to each other *and* to their `cand_*` twins — because their 04 CRBN contact
sets are identical and 05 never sees the ligand. (`cand_779` and `cand_1078` do the same.)

### 8b. Do NOT merge Tyrone's `completed_results_with_controls` CSV — measured, not assumed
That file (30 candidates + controls) is real data, but its `mean_best_10` **cannot be ranked
against v1 numbers**. It carries `protocol_version = ppil4-lock-v3-45ee93aea3ab`. Three
independent checks, in increasing order of directness:

1. **The §5 acceptance test fails.** Same FPFT-2216 through both wrappers:

   | | Ryan v1 | Tyrone v3 |
   |---|---|---|
   | `mean_best_10` | -111.382 | -115.719 |
   | `dockq_vs_9dwv` | **0.376** | **0.042** |
   | cluster | 9 | 5 |

   §5 allows DockQ within ~0.05 and requires the same cluster. This is 0.334 off, in a
   different cluster. Under §5's own rule, v3 candidate runs are not trusted against v1.

2. **Measured on three molecules run under BOTH protocols** — §8a's cross-check set is
   exactly Tyrone's Candidates 118, 45 and 6:

   | molecule | Ryan v1 | Tyrone v3 | delta |
   |---|---|---|---|
   | Candidate 118 | -104.011 | -119.781 | -15.77 |
   | Candidate 45 | -98.973 | -119.045 | -20.07 |
   | Candidate 6 | -108.564 | -116.267 | -7.70 |

   v3 scores the *same molecule* **14.5 units more negative on average, sd 6.3**. The spread
   is as large as the effect, so this is **not** a constant offset that could be subtracted
   out — there is no calibration that rescues a merge.

3. **The pools are different size classes.** Ryan's 18 are 30–34 heavy atoms and every one
   carries the phenylpiperazine PPIL4-facing arm. All 30 of Tyrone's are 20–24 heavy atoms
   with no arm — **zero overlap**. The HADDOCK score is size- and interface-dominated
   (README, Limitations), so the two pools are not on one scale even before the protocol
   difference. (Also: **26 of his 30 are already in 01's library** under a `cand_*` name —
   they were generated and scored by 02/03 and didn't make the cut.)

A raw merge on `mean_best_10` would have put 4 of his molecules into a top-20 on that
~14.5-unit artifact. **The only valid route in is re-docking under `ppil4_lock_v1`** (§8a).
Dock more of them and they enter the buy list automatically via 11.

---

## 9. Implementation corrections (found while wiring and running `ppil4_lock_v1`)
1. **Seed param is `iniseed`, not `seed`.** HADDOCK3 2026.7.0 has no top-level global `seed`;
   the sampling seed is `iniseed` (default 917) on `[rigidbody]`. A top-level `seed = 42` is
   silently ignored. Also: even with `iniseed` pinned, multi-core runs are **not** byte-identical
   (CNS orders parallel models nondeterministically) — exact reproduction needs `ncores = 1`.
2. **The degron pocket came from the co-crystal, not the receptor.** `CRBN_receptor_thalidomide_Ryan.pdb`
   is apo (Zn only, no thalidomide), so "within 5 Å of the bound thalidomide" in that file is empty.
   The fixed CRBN residues were derived from the co-crystal `CRBN-Thalidomide-SALL4_(Ryan).pdb`
   (ligand EF2) and confirmed present in the apo receptor with matching numbering.
3. **`dockq_vs_9dwv` was specified but nothing ever wrote it.** 08 hardcodes it blank and 09
   only printed its value, so the column stayed empty however many times either was run. 09 now
   persists it for every molecule — and doing that is what revealed the column fails its own
   calibration (§6), which is why the schema keeps it but the ranking doesn't use it.
4. **The run tree must leave the synced Desktop** (§3a). Not a scoring change, but it cost a
   full 18-candidate overnight batch before it was diagnosed.
5. **Resume by detecting finished docks, not ledger names.** A crashed or interrupted run
   leaves a `-` ledger row and a partial run dir; trusting the name skips it forever. 08 now
   treats a candidate as done iff `9_caprieval/capri_ss.tsv` exists with ≥1 scored model, and
   recovers results from disk when the ledger has lost them.
6. **`run_dir` must be ASCII, not just unsynced.** HADDOCK3 validates `run_dir` against
   `[a-zA-Z0-9._-/\]` and hard-fails on anything else. This project's own path is now
   non-ASCII (`.../校外/ISEF`), so *any* run tree under the project root cannot start —
   `ConfigurationError: The 'run_dir' parameter can only have [...] characters`. 08 was
   already immune (its tree is under `~/haddock_runs` for the iCloud reason in §3a); 05 was
   not, and was moved to the same root. Same fix, second independent reason.
7. **The buy list has no generator script.** `final_buy_list_lock_v1_(Ryan).csv` was assembled
   by hand, so two of its derived columns cannot be recomputed from source: `structure_score_z`
   and `rule_a_upstream_crbn_pose`. 11 re-derives both from written-down definitions
   (§6's named structure axis; "top-half Vina rank AND pose overlap ≥ 0.70") — the
   reconstruction matches all 10 committed `rule_a` values and reproduces the committed
   `structure_score_z` **ordering** but not its values (max deviation 0.42). `dock_z`,
   `vina_rank_of_18` and `rule_c` do reproduce exactly, 10/10. **Do not mix the two
   `structure_score_z` columns** — recompute the 18 with 11's `--recompute-18` if you want one
   consistent column.
