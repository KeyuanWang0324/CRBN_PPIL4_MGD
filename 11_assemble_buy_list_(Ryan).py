"""
Assemble the final buy list from every molecule docked under ppil4_lock_v1 --
Ryan's 18 funnel candidates plus the cross-check molecules from Tyrone's screen
(PROTOCOL_LOCK.md section 8a) -- ranked by 08's HADDOCK score, in one fixed
column schema.

WHY THIS SCRIPT EXISTS. Every number in the buy list comes from a different
stage, and the original final_buy_list_lock_v1_(Ryan).csv was assembled by
hand. That is fine for a one-off, but it means a molecule cannot be added
without redoing that assembly, and the derived columns cannot be recomputed
identically. This writes those definitions down and applies them uniformly.

=============================================================================
ONE PROTOCOL ONLY -- WHY TYRONE'S completed_results_with_controls CSV IS NOT
MERGED IN HERE
=============================================================================
That file is real data and it is NOT ignored -- but its scores cannot be
ranked against these. It was produced under `ppil4-lock-v3-45ee93aea3ab`,
not `ppil4_lock_v1`, and PROTOCOL_LOCK.md Rule 0 is that a HADDOCK score is
only comparable within one exact protocol. Three independent checks say the
difference is real and large:

  1. Section 5's acceptance test FAILS between the two wrappers. Same
     FPFT-2216: DockQ 0.376 (v1) vs 0.042 (v3) against a ~0.05 tolerance,
     landing in different clusters (9 vs 5).
  2. Measured directly on three molecules run under BOTH protocols
     (Tyrone's Candidates 118, 45, 6 = the cross-check set here):
     v3 scores the SAME molecule 7.7 / 15.8 / 20.1 units more negative,
     mean -14.5, sd 6.3. The spread is as large as the effect, so this is
     not a constant offset that could be subtracted out.
  3. The two pools are different size classes, and the HADDOCK score is
     size/interface-dominated (README, Limitations): Ryan's 18 are 30-34
     heavy atoms and all carry the phenylpiperazine PPIL4-facing arm;
     all 30 of Tyrone's are 20-24 heavy atoms with no arm. Zero overlap.

So Tyrone's molecules enter this list only by being re-docked under
ppil4_lock_v1 (04 -> 08 -> 09, see PROTOCOL_LOCK.md section 8a), which is
what the cross-check set is. Merging his v3 numbers directly would promote
molecules on a ~14.5-unit protocol artifact. Dock more of them under v1 and
they will appear here automatically.

WHERE EACH COLUMN COMES FROM
  measured, read straight from the pipeline's own outputs:
    haddock_mean_best10/haddock_best_model/sd_best_10/n_models
                                    <- 08's two locked-schema result CSVs
    bsa/air/vdw_mean_best10, top_cluster_fraction
                                    <- 08's run tree, <mol>/run1/9_caprieval/capri_ss.tsv
    vina_crbn_affinity/vina_combined_affinity/degron_contact_overlap
                                    <- 04's screening + cross-check CSVs
    p_crbn_glue                     <- 02_crbn_glue_scores_for_03_(Ryan).csv
    legacy_05_dockq/legacy_05_fnat  <- 05's funnel + cross-check CSVs

  derived here: dock_z, structure_score_z, the three ranks, rules A/B/C,
  rules_passed, and (for molecules not already in the committed buy list)
  decision, selection_axis and reason.

DERIVED COLUMNS -- read before comparing against the committed buy list.
  dock_z            EXACT. -(mean_best_10 - mean)/sd over Ryan's 18.
                    Reproduces all 10 committed rows 10/10.
  vina_rank_of_18   EXACT. Ranked by crbn_affinity (NOT combined) -- the
                    committed column is monotonic in crbn_affinity for all
                    10 rows and is not monotonic in combined. 10/10.
  rule_c            EXACT. mean_best_10 < -125.610, the Z6466608628 bar
                    (PROTOCOL_LOCK.md section 7). 10/10.
  structure_score_z RE-DERIVED, not exact. The committed column's formula is
                    not in the repo. PROTOCOL_LOCK.md section 6 names the
                    axis -- "BSA / cluster convergence / AIR / sd_best_10" --
                    so this is the equal-weight mean of those four z-scores
                    (BSA and cluster fraction as-is; AIR and sd negated).
                    Reproduces the committed ORDERING closely but not its
                    values (max deviation 0.42). This script recomputes it
                    for every row so the column is internally consistent --
                    which means the 10 committed rows' values change here.
                    That is deliberate: one function, one column.
  rule_a            RE-DERIVED, not verified. "04's Vina pose is a credible
                    CRBN-pocket pose": top-half Vina rank among the 18 AND
                    pose overlap >= 0.70. Reproduces all 10 committed
                    yes/no values, but it is a 2-parameter rule fitted to 10
                    points, not a definition recovered from source.

WHAT IS PRESERVED RATHER THAN REGENERATED. For a molecule already in the
committed buy list, `decision`, `selection_axis`, `chemotype` and `reason`
are carried over VERBATIM. Those encode Ryan's judgment -- in particular the
ADD rows were hand-picked to cover an axis extreme regardless of how many
rules they pass, and regenerating them would silently discard that. New rows
get a rule-derived decision (3 rules -> BUY, 2 -> CONSIDER, 0-1 -> NO-BUY)
and a generated reason. --report prints every row whose recomputed
rules_passed disagrees with its carried-over decision.

Run with the SYSTEM python (needs rdkit), after 04, 05, 08 and 09:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "11_assemble_buy_list_(Ryan).py"
    # --top N    how many candidates to keep (default 20)
    # --report   show deltas vs the committed buy list, and the full ranking
"""
import argparse
import collections
import csv
import os
import statistics
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

FINALISTS_CSV = os.path.join(SCRIPT_DIR, "08_final_ternary_with_ligand_results_(Ryan).csv")
CROSSCHECK_CSV = os.path.join(SCRIPT_DIR, "08_crosscheck_results_(Ryan).csv")
COMMITTED_BUY_LIST = os.path.join(SCRIPT_DIR, "final_buy_list_lock_v1_(Ryan).csv")
VINA_CROSSCHECK_CSV = os.path.join(SCRIPT_DIR, "04_crosscheck_vina_scores_(Ryan).csv")
VINA_FUNNEL_CSV = os.path.join(SCRIPT_DIR, "04_vina_screening_scores_for_05_(Ryan).csv")
GLUE_SCORES_CSV = os.path.join(SCRIPT_DIR, "02_crbn_glue_scores_for_03_(Ryan).csv")
LEGACY_05_CSV = os.path.join(SCRIPT_DIR, "05_ternary_docking_scores_for_06_(Ryan).csv")
LEGACY_05_CROSSCHECK_CSV = os.path.join(SCRIPT_DIR, "05_crosscheck_ternary_scores_(Ryan).csv")
# 05 only ever docked the top TOP_FRACTION (20%) by Vina, so two 08 finalists
# (cand_779, cand_1078) had no legacy_05_* value; this is their backfill run.
LEGACY_05_BACKFILL_CSV = os.path.join(SCRIPT_DIR, "05_backfill_ternary_scores_(Ryan).csv")

OUT_CSV = os.path.join(SCRIPT_DIR, "final_buy_list_lock_v1_top20_(Ryan).csv")
OUT_CROSSCHECK_CSV = os.path.join(SCRIPT_DIR, "final_buy_list_crosscheck_(Ryan).csv")
# The order sheet is built by 12, which reads this script's buy list directly --
# one owner, so the two cannot drift apart.

RUN_DIR_BASES = [
    os.path.join(os.path.expanduser("~"), "haddock_runs", "haddock3_ternary_with_ligand_run"),
    os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run"),
]
VINA_OUT_DIR = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_novel_candidate")
CRBN_ACTIVE_FIXED = os.path.join(SCRIPT_DIR, "crbn_active_FIXED.txt")

PROTOCOL_VERSION = "ppil4_lock_v1"
RYAN_POOL = "ryan18_lock_v1"
CROSS_POOL = "tyrone_crosscheck_lock_v1"
Z6466608628_BAR = -125.610      # PROTOCOL_LOCK.md section 7: the usable bar
NOISE_FLOOR = -112.202          # negctrl_no_glutarimide -- below this is meaningless
RULE_A_MIN_POSE_OVERLAP = 0.70

# Cross-check molecules: Tyrone's, re-docked under this protocol. tyrone_id is
# his numbering; ryan_twin is the identical-SMILES molecule already in 01's
# library (which was never docked in 08, so nothing is double-counted).
CROSSCHECK_META = {
    "tyrone_118": {"tyrone_id": "Candidate 118", "ryan_twin": "cand_1981",
                   "chemotype": "cyano/methoxy phthalimide (di-oxo), no PPIL4 arm"},
    "tyrone_45": {"tyrone_id": "Candidate 45", "ryan_twin": "cand_467",
                  "chemotype": "bromo/methyl phthalimide (di-oxo), no PPIL4 arm"},
    "tyrone_6": {"tyrone_id": "Candidate 6", "ryan_twin": "cand_85",
                 "chemotype": "ethyl phthalimide (di-oxo), no PPIL4 arm"},
}

# Prose chemotypes for the docked candidates that were never in the committed
# buy list, in the same style as its hand-written column. Read off the SMILES:
# "mono-oxo" = isoindolinone (one ring C=O), "di-oxo" = phthalimide (two).
NEW_CHEMOTYPES = {
    "cand_2392": "carboxylic-acid isoindolinone (mono-oxo), phenylpiperazine arm",
    "cand_234": "cyano isoindolinone (mono-oxo), phenylpiperazine arm",
    "cand_445": "chloro phthalimide (di-oxo), phenylpiperazine arm",
    "cand_779": "bromo phthalimide (di-oxo), phenylpiperazine arm",
    "cand_1092": "chloro isoindolinone (mono-oxo), phenylpiperazine arm",
    "cand_1078": "fluoro phthalimide (di-oxo), phenylpiperazine arm",
    "cand_72": "unsubstituted isoindolinone (mono-oxo), phenylpiperazine arm",
    "cand_672": "methyl isoindolinone (mono-oxo), phenylpiperazine arm",
}

COLUMNS = [
    "name", "run_name", "source_pool", "role", "decision", "selection_axis", "smiles", "chemotype",
    "has_glutarimide", "haddock_mean_best10", "haddock_best_model", "sd_best_10",
    "dock_rank_of_18", "dock_z", "structure_score_z", "structure_rank_of_18",
    "bsa_mean_best10", "top_cluster_fraction", "air_mean_best10", "vdw_mean_best10",
    "n_models", "vina_crbn_affinity", "vina_combined_affinity", "vina_rank_of_18",
    "degron_contact_overlap", "p_crbn_glue", "legacy_05_dockq", "legacy_05_fnat",
    "rule_a_upstream_crbn_pose", "rule_b_structure_z_ge_0", "rule_c_beats_Z6466608628",
    "rules_passed", "protocol_version", "reason",
]

STRUCTURE_AXIS = [
    ("bsa_mean_best10", 1),       # bigger buried interface is better
    ("top_cluster_fraction", 1),  # more of the sampling converged on one mode
    ("air_mean_best10", -1),      # less restraint violation is better
    ("sd_best_10", -1),           # tighter top-10 spread is better
]


def display_name(run_name):
    """The name a human reads on the buy list, keyed to whose screen the
    molecule came from: Ryan's `cand_708` -> `RW_708`, Tyrone's `tyrone_6` ->
    `TZ_6`. Controls keep their own descriptive names -- they are not
    candidates and belong to neither screen.

    This is a DISPLAY name only. Every run directory, every 08/09 CSV and every
    upstream file still uses the pipeline name, so the buy list also carries a
    `run_name` column -- without it a row could not be traced back to the run
    that produced it."""
    if run_name.startswith("cand_"):
        return "RW_" + run_name[len("cand_"):]
    if run_name.startswith("tyrone_"):
        return "TZ_" + run_name[len("tyrone_"):]
    return run_name


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def capri_ss_path(molecule):
    for base in RUN_DIR_BASES:
        path = os.path.join(base, molecule, "run1", "9_caprieval", "capri_ss.tsv")
        if os.path.exists(path):
            return path
    return None


def structure_metrics(molecule):
    """The structure-axis measurements for one molecule, over the SAME 10
    best-scoring single models the locked schema's score stats use, plus
    cluster convergence over the whole run.

    top_cluster_fraction is the LARGEST cluster's share of all scored models --
    how much of the sampling converged on one binding mode -- not the rank-1
    model's own cluster share. Verified against the committed buy list
    (cand_1987 -> 0.805, cand_1959 -> 0.655)."""
    path = capri_ss_path(molecule)
    if path is None:
        return None
    with open(path, newline="") as f:
        recs = list(csv.DictReader(f, delimiter="\t"))
    if not recs:
        return None
    recs.sort(key=lambda r: float(r["score"]))
    best10 = recs[:10]
    counts = collections.Counter(r["cluster_id"] for r in recs)
    return {
        "bsa_mean_best10": round(statistics.fmean(float(r["bsa"]) for r in best10), 1),
        "air_mean_best10": round(statistics.fmean(float(r["air"]) for r in best10), 1),
        "vdw_mean_best10": round(statistics.fmean(float(r["vdw"]) for r in best10), 2),
        "top_cluster_fraction": round(max(counts.values()) / len(recs), 3),
        "n_models": len(recs),
    }


def z_params(pool, key, sign=1):
    """(mean, sd, sign) for one metric over the reference pool. sign=-1 where
    lower is better, so every resulting z reads 'higher = better'.

    Returned as parameters rather than finished scores so a molecule OUTSIDE
    the reference pool can be scored against the same fixed distribution. That
    is the point: a cross-check molecule's z must mean 'this far from Ryan's
    18', and must not move when another cross-check molecule is added."""
    values = [m[key] for m in pool.values()]
    return statistics.fmean(values), statistics.stdev(values), sign


def apply_z(params, value):
    mean, sd, sign = params
    return sign * (value - mean) / sd


def structure_axis_params(pool):
    return [z_params(pool, key, sign) for key, sign in STRUCTURE_AXIS]


def structure_score(params, metrics):
    """PROTOCOL_LOCK.md section 6's structure axis, equal weight. See DERIVED
    COLUMNS: a written-down reconstruction, not the committed column's formula."""
    return statistics.fmean(
        apply_z(p, metrics[key]) for p, (key, _) in zip(params, STRUCTURE_AXIS)
    )


def insertion_rank(value, reference_values, higher_is_better):
    """The rank `value` takes against reference_values (1 = best). For a
    molecule inside the reference pool this is its own rank; for one outside
    it is the rank it WOULD take, so the *_rank_of_18 columns stay meaningful."""
    if higher_is_better:
        return sum(1 for v in reference_values if v > value) + 1
    return sum(1 for v in reference_values if v < value) + 1


def crbn_pose_contacts(molecule):
    """The CRBN residues 04's selected Vina pose touches -- written by 04 into
    <molecule>/crbn_contacts.txt."""
    path = os.path.join(VINA_OUT_DIR, molecule, "crbn_contacts.txt")
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        first = f.readline()
    return {int(x) for x in first.split()} if first.strip() else set()


def has_glutarimide(smiles):
    """The CRBN degron itself: the 2,6-dioxopiperidine (glutarimide) ring.
    A molecule without it has no CRBN handle -- that is what
    negctrl_no_glutarimide tests."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return ""
    return mol.HasSubstructMatch(Chem.MolFromSmarts("O=C1NC(=O)CCC1"))


def main():
    parser = argparse.ArgumentParser(description="Assemble the ppil4_lock_v1 buy list.")
    parser.add_argument("--top", type=int, default=20, help="candidates to keep (default 20)")
    parser.add_argument("--report", action="store_true",
                        help="print the full ranking and every delta vs the committed buy list")
    args = parser.parse_args()

    finalists = read_csv(FINALISTS_CSV)
    if not finalists:
        sys.exit(f"{os.path.basename(FINALISTS_CSV)} not found or empty -- run 08 first.")
    crosscheck = [r for r in read_csv(CROSSCHECK_CSV) if r["mean_best_10"] not in ("", "-", None)]

    # --- Reference distribution: Ryan's 18. Every z and rank is measured
    # against these and only these, so adding cross-check molecules never
    # shifts an existing row's numbers.
    pool = {}
    for row in finalists:
        metrics = structure_metrics(row["molecule"])
        if metrics is None:
            print(f"NOTE: no run dir for {row['molecule']} -- excluded from the reference distribution.")
            continue
        metrics["mean_best_10"] = float(row["mean_best_10"])
        metrics["sd_best_10"] = float(row["sd_best_10"])
        pool[row["molecule"]] = metrics
    if len(pool) < 3:
        sys.exit(f"Only {len(pool)} reference candidates have run dirs -- too few for a z-score.")

    struct_params = structure_axis_params(pool)
    dock_params = z_params(pool, "mean_best_10", -1)
    struct_z_ref = sorted((structure_score(struct_params, m) for m in pool.values()), reverse=True)
    dock_ref = [m["mean_best_10"] for m in pool.values()]

    vina_all = {r["name"]: r for r in read_csv(VINA_FUNNEL_CSV)}
    vina_all.update({r["name"]: r for r in read_csv(VINA_CROSSCHECK_CSV)})
    vina_ref = [float(vina_all[n]["crbn_affinity"]) for n in pool if n in vina_all]

    glue_by_name = {r["name"]: r["p_crbn_glue"] for r in read_csv(GLUE_SCORES_CSV)}
    smiles_by_name = {r["name"]: r["smiles"] for r in read_csv(GLUE_SCORES_CSV)}
    legacy05 = {r["name"]: r for r in read_csv(LEGACY_05_CSV)}
    legacy05.update({r["name"]: r for r in read_csv(LEGACY_05_CROSSCHECK_CSV)})
    legacy05.update({r["name"]: r for r in read_csv(LEGACY_05_BACKFILL_CSV)})
    degron_pocket = {int(x) for x in open(CRBN_ACTIVE_FIXED).readline().split()}
    committed = {r["name"]: r for r in read_csv(COMMITTED_BUY_LIST)}

    rows = []
    for row in finalists + crosscheck:
        name = row["molecule"]
        metrics = structure_metrics(name)
        if metrics is None:
            continue
        metrics["mean_best_10"] = float(row["mean_best_10"])
        metrics["sd_best_10"] = float(row["sd_best_10"])
        meta = CROSSCHECK_META.get(name, {})
        is_cross = name in CROSSCHECK_META

        struct_z = structure_score(struct_params, metrics)
        dock_z = apply_z(dock_params, metrics["mean_best_10"])
        vina_row = vina_all.get(name)
        if vina_row is None:
            sys.exit(f"No Vina row for {name} -- run 04 first.")
        vina_crbn = float(vina_row["crbn_affinity"])
        vina_rank = insertion_rank(vina_crbn, vina_ref, higher_is_better=False)
        # A cross-check molecule's own SMILES lives under its identical-SMILES twin.
        lookup = meta.get("ryan_twin", name)
        smiles = smiles_by_name.get(lookup, "")

        rule_a = vina_rank <= (len(vina_ref) + 1) // 2 and float(vina_row["overlap"]) >= RULE_A_MIN_POSE_OVERLAP
        rule_b = struct_z >= 0
        rule_c = metrics["mean_best_10"] < Z6466608628_BAR
        rules_passed = sum((rule_a, rule_b, rule_c))

        dock_rank = insertion_rank(metrics["mean_best_10"], dock_ref, higher_is_better=False)
        struct_rank = insertion_rank(struct_z, struct_z_ref, higher_is_better=True)

        prior = committed.get(name)
        if prior and prior["role"] == "candidate":
            decision, axis = prior["decision"], prior["selection_axis"]
            chemotype, reason = prior["chemotype"], prior["reason"]
        else:
            decision = {3: "BUY", 2: "CONSIDER"}.get(rules_passed, "NO-BUY")
            axis = ("both" if dock_z >= 0 and struct_z >= 0 else
                    "docking" if dock_z >= 0 else
                    "structure" if struct_z >= 0 else "neither")
            chemotype = meta.get("chemotype") or NEW_CHEMOTYPES.get(name, "")
            beats = ("beats the Z6466608628 glue bar" if rule_c else
                     f"above the {NOISE_FLOOR} noise floor -- not distinguishable from a "
                     "degron-dead negative" if metrics["mean_best_10"] > NOISE_FLOOR else
                     "below the noise floor but short of the Z6466608628 bar")
            hits = crbn_pose_contacts(name) & degron_pocket
            origin = (f"{meta['tyrone_id']} from Tyrone's guided screen, re-docked under "
                      f"{PROTOCOL_VERSION}. " if is_cross else "")
            tail = (f" Same molecule as Ryan's {meta['ryan_twin']} (identical SMILES), which was "
                    "never docked in 08." if is_cross else "")
            reason = (
                f"{origin}HADDOCK {metrics['mean_best_10']:.1f} ({beats}); dock #{dock_rank} of "
                f"{len(dock_ref) + (1 if is_cross else 0)}. Structure axis z={struct_z:+.2f} "
                f"(#{struct_rank}): BSA {metrics['bsa_mean_best10']:.0f}, "
                f"{metrics['top_cluster_fraction']:.0%} of models in one cluster, AIR "
                f"{metrics['air_mean_best10']:.0f}. Vina #{vina_rank} ({vina_crbn:.2f} CRBN), "
                f"pose hits {len(hits)}/{len(degron_pocket)} degron-pocket residues. "
                f"Passes {rules_passed}/3 rules -> {decision}.{tail}"
            )

        rows.append({
            "name": display_name(name),
            "run_name": name,
            "source_pool": CROSS_POOL if is_cross else RYAN_POOL,
            "role": "crosscheck" if is_cross else "candidate",
            "decision": decision, "selection_axis": axis, "smiles": smiles,
            "chemotype": chemotype,
            "has_glutarimide": has_glutarimide(smiles),
            "haddock_mean_best10": row["mean_best_10"],
            "haddock_best_model": row["best_haddock"],
            "sd_best_10": row["sd_best_10"],
            "dock_rank_of_18": dock_rank,
            "dock_z": round(dock_z, 3),
            "structure_score_z": round(struct_z, 3),
            "structure_rank_of_18": struct_rank,
            "bsa_mean_best10": metrics["bsa_mean_best10"],
            "top_cluster_fraction": metrics["top_cluster_fraction"],
            "air_mean_best10": metrics["air_mean_best10"],
            "vdw_mean_best10": metrics["vdw_mean_best10"],
            "n_models": metrics["n_models"],
            "vina_crbn_affinity": vina_crbn,
            "vina_combined_affinity": round(float(vina_row["combined_affinity"]), 3),
            "vina_rank_of_18": vina_rank,
            "degron_contact_overlap": round(float(vina_row["overlap"]), 3),
            "p_crbn_glue": glue_by_name.get(lookup, ""),
            "legacy_05_dockq": legacy05.get(name, {}).get("dockq", ""),
            "legacy_05_fnat": legacy05.get(name, {}).get("fnat", ""),
            "rule_a_upstream_crbn_pose": "yes" if rule_a else "no",
            "rule_b_structure_z_ge_0": "yes" if rule_b else "no",
            "rule_c_beats_Z6466608628": "yes" if rule_c else "no",
            "rules_passed": rules_passed,
            "protocol_version": PROTOCOL_VERSION,
            "reason": reason,
        })

    rows.sort(key=lambda r: float(r["haddock_mean_best10"]))
    kept, dropped = rows[:args.top], rows[args.top:]

    # Controls carried over verbatim -- they are the calibration set, not
    # candidates, and nothing here recomputes them. They keep their own names
    # (they belong to neither screen), so run_name == name.
    controls = [r for r in read_csv(COMMITTED_BUY_LIST) if r["role"] != "candidate"]
    for c in controls:
        c["run_name"] = c["name"]

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(kept)
        writer.writerows(controls)
    print(f"Wrote {OUT_CSV}\n  {len(kept)} candidates (of {len(rows)} docked under "
          f"{PROTOCOL_VERSION}) + {len(controls)} controls")
    if dropped:
        print(f"  below the cut: {', '.join(r['name'] for r in dropped)}")

    cross_rows = [r for r in rows if r["role"] == "crosscheck"]
    if cross_rows:
        with open(OUT_CROSSCHECK_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(cross_rows)
        print(f"Wrote {OUT_CROSSCHECK_CSV} ({len(cross_rows)} cross-check molecules)")

    print(f"\n{'#':>3} {'molecule':<12} {'haddock':>9} {'dock_z':>7} {'struct_z':>9} "
          f"{'rules':>5}  {'decision':<9} source")
    print("-" * 74)
    for i, r in enumerate(kept, 1):
        flag = "" if float(r["haddock_mean_best10"]) < NOISE_FLOOR else "  << above noise floor"
        print(f"{i:>3} {r['name']:<12} {float(r['haddock_mean_best10']):>9.3f} {r['dock_z']:>7.3f} "
              f"{r['structure_score_z']:>9.3f} {r['rules_passed']:>5}  {r['decision']:<9} "
              f"{r['source_pool']}{flag}")
    print(f"\nBars: Z6466608628 {Z6466608628_BAR} (rule C) | noise floor {NOISE_FLOOR} "
          "(negctrl_no_glutarimide).")

    if args.report:
        print("\n== Deltas vs the committed buy list ==")
        for r in kept:
            prior = committed.get(r["run_name"])
            if not prior or prior["role"] != "candidate":
                print(f"  {r['name']:<12} NEW ROW ({r['decision']}, {r['rules_passed']}/3 rules)")
                continue
            d = abs(float(r["structure_score_z"]) - float(prior["structure_score_z"]))
            note = f"structure_score_z {prior['structure_score_z']} -> {r['structure_score_z']} (d={d:.3f})"
            if r["rules_passed"] != int(prior["rules_passed"]):
                note += f" | rules_passed {prior['rules_passed']} -> {r['rules_passed']}"
            if r["decision"] != prior["decision"]:
                note += f" | decision {prior['decision']} -> {r['decision']}"
            print(f"  {r['name']:<12} {note}")
        mismatched = [r for r in kept
                      if committed.get(r["run_name"], {}).get("role") == "candidate"
                      and r["rules_passed"] != int(committed[r["run_name"]]["rules_passed"])]
        if mismatched:
            print(f"\n  {len(mismatched)} carried-over decision(s) no longer match their recomputed "
                  "rules_passed (the recomputed structure_score_z moved rule B). The decision is "
                  "kept as Ryan set it -- see WHAT IS PRESERVED in the docstring.")


if __name__ == "__main__":
    main()
