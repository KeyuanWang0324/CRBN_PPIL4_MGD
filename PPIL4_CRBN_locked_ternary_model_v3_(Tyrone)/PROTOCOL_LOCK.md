# PROTOCOL_LOCK — the one docking protocol you both run

Ryan and Tyrone: this file exists so your two rankings become **comparable**. Right now they aren't, because you ran two different protocols. This locks **one**.

## Rule 0 — share the artifact, not a description

A HADDOCK score is only comparable **within one exact protocol**. A prose description (even this file) under-determines the run: sampling, restraints, seed, scoring all have to match to the setting. If each of you hands a description to a different AI and lets it regenerate the code, you get **different code → different protocol → uncomparable scores again.**

So the thing you share and pin is the **executable artifact** below — the config, the restraint files, the inputs, the versions, the seed — **not** a paragraph telling an AI to rebuild them. Commit these files to the repo; both of you run *those files*, not your own regenerated copies.

Your **wrapper is allowed to differ** (see the last section) — that's even useful. The **cfg / restraints / inputs / seed are not.**

---

## 1. The locked config (`guided_ppil4_haddock.cfg`)

This is Ryan's `08` protocol (the more rigorous one), with **two changes** marked `### LOCK`. Commit this exact file; do not regenerate it.

```toml
run_dir = "run_ppil4_ternary"
ncores = 8                     # speed only — may differ per machine, does NOT change scores
seed = 42                      ### LOCK #1: pin the seed so the run is reproducible (see §5)

molecules = [
    "CRBN_receptor_thalidomide_Ryan.pdb",   # chain A
    "ligand_prep/<candidate>_prodrg_fixed.pdb",  # chain C, resnum 900 (per-candidate, deterministic)
    "PPIL4_chainB.pdb"                        # chain B
]

[topoaa]
ligand_top_fname   = "ligand_prep/<candidate>_prodrg.top"     # PRODRG explicit — NOT autotoppar
ligand_param_fname = "ligand_prep/<candidate>_prodrg.param"

[rigidbody]
ambig_fname        = "ambig_FIXED.tbl"       ### LOCK #2: one fixed restraint file for ALL candidates
ligand_top_fname   = "ligand_prep/<candidate>_prodrg.top"
ligand_param_fname = "ligand_prep/<candidate>_prodrg.param"
sampling           = 1000

[caprieval]

[seletop]
select = 200

[flexref]
ambig_fname        = "ambig_FIXED.tbl"
ligand_top_fname   = "ligand_prep/<candidate>_prodrg.top"
ligand_param_fname = "ligand_prep/<candidate>_prodrg.param"

[caprieval]

[emref]
ambig_fname        = "ambig_FIXED.tbl"
ligand_top_fname   = "ligand_prep/<candidate>_prodrg.top"
ligand_param_fname = "ligand_prep/<candidate>_prodrg.param"

[caprieval]

[clustfcc]

[caprieval]
```

**Ligand topology = PRODRG, explicit** (`ligand_top_fname`/`ligand_param_fname` on every module). Do **not** use `[topoaa] autotoppar = true` — it renames the ligand atoms and desyncs them from the coordinates (Ryan documented this; it's the bug in the `autotoppar` path). Run PRODRG once per candidate, keep its renamed atoms, and pass the `.top`/`.param` explicitly.

---

## 2. The fixed restraints (`ambig_FIXED.tbl`) — LOCK #2, the important one

The restraint set must be **identical for every candidate**, so the only thing that changes between runs is the molecule. Right now Ryan's script derives the **CRBN-side** active residues from *each candidate's own Vina pose* — that leaks the Vina step's noise into the ranking. Fix it to one set.

- **PPIL4 side (already fixed — keep it):** chain B active residues `249, 250, 273, 275, 276, 277, 278, 279` (the real RRM interface from 9DWV / FPFT-2216).
- **CRBN side (change to fixed):** the glutarimide/degron pocket — the same pocket for every candidate, because every candidate binds CRBN through its glutarimide. Derive it **once** from the co-crystal thalidomide in `CRBN_receptor_thalidomide_Ryan.pdb` (chain A residues within 5 Å of the bound thalidomide), write it to `crbn_active_FIXED.txt`, and **commit that file**. Do not recompute per candidate.
- Build `ambig_FIXED.tbl` **once** from `crbn_active_FIXED` + `ppil4_active_FIXED` (each `active_passive_to_ambig` against the ligand, resnum 900), commit it, and point all three modules at it.

With both sides fixed, `ambig_FIXED.tbl` is the same file for all candidates and both of you.

---

## 3. Fixed inputs + versions (commit / record these)

| Artifact | What to pin |
|---|---|
| CRBN receptor | `CRBN_receptor_thalidomide_Ryan.pdb` (chain A, thalidomide bound) |
| PPIL4 | `PPIL4_alphafold_(Ryan).pdb` → chain B (`PPIL4_chainB.pdb`) |
| Ligand starting pose | Vina pose transformed into the 9DWV CRBN frame (04's recipe) — **same recipe for candidates AND controls** |
| HADDOCK3 | record `haddock3 --version` in the repo; both use the same version |
| PRODRG | the bundled `haddock/prodrg/prodrg_<arch>` binary |

---

## 4. What may differ, and what may not

| ✅ May differ (your own code — even good to differ) | ❌ Must be byte-identical |
|---|---|
| launcher / driver script | `guided_ppil4_haddock.cfg` |
| CSV parsing, tabulation | `ambig_FIXED.tbl` + the two `_active_FIXED.txt` |
| plotting, progress bars | input PDBs + ligand-pose recipe |
| which AI wrote the wrapper | HADDOCK3 version, seed, PRODRG |

If it feeds the cfg the same inputs, it produces the same numbers — regardless of who (or which AI) wrote it.

---

## 5. Acceptance test (golden reference) — do this BEFORE any candidate

You need one molecule with a known expected result under the **locked** protocol, so you can prove your setup matches before spending hours on candidates. Use the positive control.

1. One of you runs **FPFT-2216** through the locked cfg. Record `mean_best_10`, `best_haddock`, `sd_best_10`, cluster id, and DockQ-vs-9DWV (from the `00` validation step). Commit it as `REFERENCE_FPFT.json`.
2. The other runs **the same FPFT-2216** first. Your numbers must reproduce the reference:
   - seed pinned + same version + `ncores = 1` → should match **exactly**;
   - with multi-core, allow a small tolerance and also confirm the **same cluster** and DockQ within ~0.05.
3. **Only if FPFT reproduces** do you trust the run and move to candidates. If it doesn't, something in inputs/version/seed still differs — fix it now, on this one cheap molecule, not after a 20-candidate run.

(Rough expectation from earlier runs: FPFT `mean_best_10 ≈ −470`. But the committed `REFERENCE_FPFT.json` — not this number — is the thing you must match.)

---

## 6. Output schema (one CSV, fixed columns)

Both of you emit exactly these columns, so the two files concatenate and sort together:

```
molecule, mean_best_10, best_haddock, sd_best_10, dockq_vs_9dwv, cluster_id, n_models, protocol_version
```

`protocol_version` = a short tag/commit hash of this locked bundle, so a stray old-protocol row can never silently mix in.

---

## 7. The good kind of "different code" — use it as a cross-check

Divergent wrappers are a **feature** if you use them right: you each write your own launcher/parser, both point at the **same locked cfg**, and both run the shared control set (FPFT + 2–3 common molecules). Then:

- **Numbers match** → the protocol is genuinely unified *and* independently reproducible. Strong result — put it in the write-up.
- **Numbers differ** → a real bug in inputs/version/seed, and you caught it on a 3-molecule set instead of after the full run.

Divergence at the wrapper level that still agrees on numbers = validation. Divergence in the cfg/restraints = the bug. Keep the first, kill the second.
