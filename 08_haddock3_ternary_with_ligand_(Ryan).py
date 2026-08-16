"""
HADDOCK3 3-body ternary docking that includes the actual drug molecule --
unlike 05/06 (which only ever dock CRBN against PPIL4, using the
candidate's Vina-derived contact residues as a proxy and discarding the
ligand's own atoms), this script docks THREE separate bodies together:
CRBN (chain A), the candidate ligand itself (chain C), and PPIL4 (chain
B). The ligand bridges the other two via real AIR restraints on both
sides, so its actual geometry/chemistry can influence whether/how the
ternary complex comes together -- not just a residue-list hint.

THIS IS NOW THE PIPELINE'S FINAL RANKING STAGE (see chat discussion):
06's own dockq turned out to be a bad ranking signal, not just something
needing a tiebreak -- since 05/06 never put the ligand in the CNS
topology, the ONLY thing that varies between candidates is which CRBN
residues happen to be "active" (from 04's Vina pose), and since every
candidate shares the same glutarimide-isoindolinone core, many candidates'
CRBN contact sets overlap heavily or are literally identical -- so 06 is
frequently solving the exact same rigid-body docking problem for
different candidates and (correctly) returns identical scores. That's not
noise to tiebreak around, it's 06's score being structurally insensitive
to candidate chemistry. This script actually puts the ligand in the CNS
force field, so its real chemistry can differentiate candidates -- so
FINAL ranking uses THIS script's own HADDOCK score (see TOP_N_KEEP below),
not dockq (self-referential, no real reference structure) and not 06's.

06 is kept upstream purely as a coarse, cheap-relative-to-this pre-filter:
a candidate completing 06's (ligand-free) run without HADDOCK3 erroring
out is weak evidence the CRBN-side restraint/setup is basically viable
before spending 45 min - 1.5 hr per candidate here. This script's
candidate pool is whatever's currently in 06's progress ledger with a
real result (06_complete_run_progress_(Ryan).csv) -- not 06's TOP_N_KEEP
finalists, since that cut is itself decided by the same insensitive dockq
this script exists to replace.

WHY THIS NEEDS ITS OWN SCRIPT (not just a flag on 06):
HADDOCK3's CNS engine needs a topology + force-field parameter file for any
non-standard residue (the ligand), which Vina/meeko's PDBQT route doesn't
produce. HADDOCK3 bundles PRODRG (haddock/prodrg/prodrg_<arch>) and can run
it automatically via `[topoaa] autotoppar = true` -- but testing that path
surfaced a real gap in HADDOCK3's own integration: PRODRG renames the
ligand's atoms to match its generated topology (e.g. "H1" -> "HG3"), and the
auto path never propagates that renaming back into the working PDB's
coordinates, so CNS aborts with "unknown coordinates for atom ...". (There
is literally a `# TODO: Check if the atom names have been changed!` comment
in haddock3's own libligand.py for this.) This script works around it by
running PRODRG once per candidate itself, keeping PRODRG's own renamed
coordinate file (DRGFIN.PDB) instead of the original input coordinates, and
passing the resulting topology/param explicitly via `ligand_top_fname`/
`ligand_param_fname` -- which must be set on EVERY module that runs its own
CNS minimization (topoaa, rigidbody, flexref, emref), not just topoaa, or
CNS aborts with "missing nonbonded Lennard-Jones parameters" instead.

This imports a few non-public haddock.libs functions (PRODRG binary
resolution, atom-name/NBONDS sanitization) rather than reimplementing them,
so behavior stays exactly consistent with whatever HADDOCK3 version is
installed -- but that does mean this script is more sensitive to a
haddock3 upgrade than 05/06 (which only ever shell out to the stable
`haddock3`/`haddock3-restraints` CLIs).

Validated with a 4-model smoke test (rigidbody -> caprieval -> seletop ->
flexref -> caprieval -> emref -> caprieval -> clustfcc -> caprieval) on
cand_5 before writing this: completed in 74s, top cluster n=3/4, dockq
0.627. At the funnel's full scale (sampling=1000, top 200 for
flexref/emref) each candidate costs ~45 min - 1.5 hr, same order as 06.

For validating this whole protocol's ability to reproduce a KNOWN real
ternary complex (as opposed to ranking library candidates), see
09_score_vs_fpft2216_reference_(Ryan).py -- it runs the real drug
(FPFT-2216) through this same script via CANDIDATE_NAME and checks the
result against the real PDB 9DWV structure. Don't use this script's own
dockq column (self-referential, no real reference) as evidence a
candidate is "good" -- only 09's positive control speaks to that.

Run in the haddock3 venv:
    source .venv-haddock3/bin/activate
    python3 "08_haddock3_ternary_with_ligand_(Ryan).py"
"""
import csv
import glob
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time

SCRIPT_START_TIME = time.time()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

if shutil.which("haddock3") is None:
    haddock_venv_bin = os.path.join(SCRIPT_DIR, ".venv-haddock3", "bin")
    haddock_python = os.path.join(haddock_venv_bin, "python3")
    env = os.environ.copy()
    env["PATH"] = haddock_venv_bin + os.pathsep + env.get("PATH", "")
    env["VIRTUAL_ENV"] = os.path.join(SCRIPT_DIR, ".venv-haddock3")
    os.execve(haddock_python, [haddock_python] + sys.argv, env)

from haddock.libs.libutil import get_prodrg_exec

VINA_OUT_DIR = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_novel_candidate")
# 06's progress ledger -- the candidate pool source. Deliberately NOT
# 06_final_ternary_results_(Ryan).csv (06's own TOP_N_KEEP finalists):
# that cut is decided by 06's dockq, the exact signal this script exists
# to replace. Any candidate with a real (non "-") result here -- i.e. it
# completed 06's ligand-free run without HADDOCK3 erroring -- is eligible.
SIX_PROGRESS_CSV = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_complete_run",
                                 "06_complete_run_progress_(Ryan).csv")
# HADDOCK3 writes THOUSANDS of tiny files per candidate and hands off between
# steps through a large io.json. If that happens inside an iCloud-synced folder
# (this project lives under ~/Desktop, which iCloud syncs on this machine), the
# iCloud file provider races HADDOCK3: it momentarily swaps the freshly-written
# io.json while caprieval is reading it, so caprieval sees an EMPTY topology and
# dies with `find_ff -> models[0].topology[0]` IndexError. That killed all 18
# candidates in the overnight run while the earlier controls (run when iCloud was
# idle) survived. Fix: run the actual docking OUTSIDE the synced Desktop. The
# home dir (~) is not synced here -- only Desktop/Documents were. The deliverables
# (FINALISTS_CSV, CONTROLS_CSV) still live in the project; only the disposable,
# regenerable run tree moves out. Override with $HADDOCK_RUN_ROOT if desired.
LOCAL_RUN_ROOT = os.environ.get("HADDOCK_RUN_ROOT") or os.path.expanduser("~/haddock_runs")
RUN_DIR_BASE = os.path.join(LOCAL_RUN_ROOT, "haddock3_ternary_with_ligand_run")
# This script's own running ledger -- every candidate attempted so far,
# win or lose. Read back in on every run so a later session only adds new
# candidates (06's progress ledger can grow over time as 06 keeps running;
# nothing here requires 06 to be "finished" first). Lives under
# docking_tmp/ (gitignored, out of the way) rather than the project root --
# it's internal bookkeeping for resume/dedup, not something to browse;
# FINALISTS_CSV below is the actual deliverable.
RESULTS_CSV = os.path.join(RUN_DIR_BASE, "08_ternary_with_ligand_progress_(Ryan).csv")
# The pipeline's actual final output -- the best TOP_N_KEEP of RESULTS_CSV
# by THIS script's own HADDOCK score (lower/more negative = better),
# rewritten after every candidate so it's always current, not gated on
# every 06 candidate finishing first (08 is far more expensive than 06;
# waiting for a full batch before showing any ranking isn't useful here).
FINALISTS_CSV = os.path.join(SCRIPT_DIR, "08_final_ternary_with_ligand_results_(Ryan).csv")

# Set this to a specific name to run only that one, bypassing the
# 06-sourced pool entirely -- e.g. a real known compound as a positive
# control (see 09_score_vs_fpft2216_reference_(Ryan).py) rather than a
# library candidate. Requires <name>/crbn_contacts.txt and
# CRBN_candidate_complex.pdb already present under VINA_OUT_DIR (run 04
# with MANUAL_CANDIDATE set to that name/SMILES first). Its result is
# NOT written to RESULTS_CSV/FINALISTS_CSV -- those are the funnel's
# shared ranked output; a one-off candidate isn't part of that ranking.
CANDIDATE_NAME = ""

# How many of RESULTS_CSV's candidates to keep in FINALISTS_CSV, ranked by
# HADDOCK score. Matches the scale 06/the rest of this pipeline has been
# using for a "final" list -- adjust freely, this isn't derived from
# anything upstream.
TOP_N_KEEP = 20

CRBN_RECEPTOR_ONLY_PDB = os.path.join(SCRIPT_DIR, "CRBN_receptor_thalidomide_Ryan.pdb")
PPIL4_SOURCE_PDB = os.path.join(SCRIPT_DIR, "PPIL4_alphafold_(Ryan).pdb")
PPIL4_PDB = os.path.join(SCRIPT_DIR, "PPIL4_chainB.pdb")   # LOCK: committed input (chain B), not regenerated per run

# --- Locked protocol (see PROTOCOL_LOCK.md) --------------------------------
# These are what make this comparable to Tyrone's runs. The restraints are one
# committed file identical for every candidate/control (so only the ligand varies),
# and the sampling seed is pinned.
AMBIG_FIXED = os.path.join(SCRIPT_DIR, "ambig_FIXED.tbl")   # LOCK #2: ONE restraint file for all candidates/controls
INISEED = 42               # LOCK #1: haddock3 'iniseed' (NOT a top-level 'seed'); ncores=1 -> byte-exact reproduction
PROTOCOL_VERSION = "ppil4_lock_v1"                          # stamped on every output row so old-protocol rows can't mix in

# --- Control molecules (calibration set) -----------------------------------
# Docked through the SAME locked protocol as candidates, but kept OUT of the
# candidate finalists -- written to their own CSV so you can see where candidates
# land relative to real glues (positive) and degron-dead / arm-less molecules
# (negative). Each needs its 04 Vina pose to already exist; 08 runs in the
# haddock venv (no Vina), so a control without a pose is SKIPPED with a message
# telling you to run 04 for it first. Negative-control SMILES for that 04 step:
#   negctrl_no_glutarimide     : CN1Cc2c(cccc2-c2ccc(N3CCNCC3)cc2)C1=O
#   negctrl_no_crbn_ppil4_arm  : COc1cscc1-c1c[nH]nn1
#   negctrl_no_ppil4_arm       : O=C1CCC(N2Cc3ccccc3C2=O)C(=O)N1
#   TM-04-064-02_positive_control : O=C(NCCCCCCCCNC1=CC=CC2=C1C(N(C2=O)C3CCC(NC3=O)=O)=O)CCN4C=NC5=C(N=CN=C54)N
RUN_CONTROLS = True   # funnel mode also docks CONTROLS below (into CONTROLS_CSV), before candidates
CONTROLS = [
    "FPFT_2216_positive_control",     # positive (also the golden reference); 04 pose exists
    "Z6466608628_positive_control",   # positive (isoindolinone-piperazine, Enamine Z6466608628); already docked -> skipped
    "TM-04-064-02_positive_control",  # positive (adenine-linked thalidomide); NEW -> run 04 for its pose first
    "negctrl_no_glutarimide",      # negative; run 04 for it first
    "negctrl_no_crbn_ppil4_arm",   # negative; run 04 for it first
    "negctrl_no_ppil4_arm",        # negative; run 04 for it first
]
CONTROLS_CSV = os.path.join(SCRIPT_DIR, "08_controls_results_(Ryan).csv")

# --- Cross-check molecules (PROTOCOL_LOCK.md section 8) --------------------
# Candidates picked by Tyrone's independent guided screen (10_*), docked here
# through Ryan's locked protocol so the two screens can be compared on the SAME
# molecules. Same locked cfg/restraints/seed as everything else -- the only
# thing that differs is where the molecule came from.
#
# Kept OUT of the candidate finalists (they aren't from 03/06's pool, and
# mixing them in would silently change what "top 20 of the funnel" means) and
# out of the controls CSV (they calibrate nothing) -- they get their own CSV,
# same locked schema. Named tyrone_<his number>; all three also exist in 01's
# library under a different name (tyrone_118 = cand_1981, tyrone_45 = cand_467,
# tyrone_6 = cand_85, identical SMILES), but none of those three was ever
# docked here, so there is no double-counting in the finalists.
# Each needs its 04 Vina pose first (04's CROSSCHECK list docks them).
RUN_CROSSCHECK = True
CROSSCHECK = [
    "tyrone_118",   # COc1ccc2c(c1C#N)C(=O)N(C1CCC(=O)NC1=O)C2=O   (cyano/methoxy phthalimide)
    "tyrone_45",    # Cc1ccc2c(c1Br)C(=O)N(C1CCC(=O)NC1=O)C2=O     (bromo/methyl phthalimide)
    "tyrone_6",     # CCc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O        (ethyl phthalimide)
]
CROSSCHECK_CSV = os.path.join(SCRIPT_DIR, "08_crosscheck_results_(Ryan).csv")

NCORES = max(1, (os.cpu_count() or 4) - 1)

STEP_PLAN = [
    (0, "topoaa", None, 0.005),
    (1, "rigidbody", 1000, 0.30),
    (2, "caprieval", None, 0.04),
    (3, "seletop", None, 0.005),
    (4, "flexref", 200, 0.30),
    (5, "caprieval", None, 0.02),
    (6, "emref", 200, 0.26),
    (7, "caprieval", None, 0.02),
    (8, "clustfcc", None, 0.035),
    (9, "caprieval", None, 0.015),
]


def estimate_progress(run_dir, step_plan):
    completed_weight = 0.0
    current = None
    for idx, name, expected, weight in step_plan:
        step_dir = os.path.join(run_dir, f"{idx}_{name}")
        if not os.path.isdir(step_dir):
            break
        next_dir = (os.path.join(run_dir, f"{idx + 1}_{step_plan[idx + 1][1]}")
                    if idx + 1 < len(step_plan) else None)
        if next_dir and os.path.isdir(next_dir):
            completed_weight += weight
            continue
        if expected:
            count = len({os.path.basename(p).split(".")[0]
                         for p in glob.glob(os.path.join(step_dir, f"{name}_*.pdb*"))})
            frac = min(count / expected, 1.0)
            current = (idx, name, count, expected)
        else:
            frac = 0.5
            current = (idx, name, None, None)
        completed_weight += weight * frac
        break
    return completed_weight, current


def write_actpass_file(active, passive, out_path):
    with open(out_path, "w") as f:
        f.write(" ".join(str(r) for r in active) + "\n")
        f.write(" ".join(str(r) for r in passive) + "\n")


def run(cmd, **kwargs):
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)


def run_with_heartbeat(cmd, run_dir=None, step_plan=None, interval=20, label=None, **kwargs):
    prefix = f"[{label}] " if label else ""
    start = time.time()
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(interval):
            elapsed = int(time.time() - start)
            if not (run_dir and step_plan):
                print(f"    ... {prefix}still running ({elapsed}s elapsed)", flush=True)
                continue
            pct, current = estimate_progress(run_dir, step_plan)
            if current is None:
                print(f"    ... {prefix}{elapsed}s elapsed, starting up", flush=True)
            else:
                idx, name, count, expected = current
                detail = f"{count}/{expected} models" if expected else "running"
                print(f"    ... {prefix}step {idx + 1}/{len(step_plan)} ({name}): "
                      f"{detail} | overall ~{pct * 100:.0f}%", flush=True)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    log_path = os.path.join(os.path.dirname(run_dir), "haddock3_stdout.log") if run_dir else None
    try:
        if log_path:
            with open(log_path, "w") as log_f:
                subprocess.run(cmd, check=True, stdout=log_f, stderr=subprocess.STDOUT, **kwargs)
        else:
            subprocess.run(cmd, check=True, **kwargs)
    except subprocess.CalledProcessError:
        if log_path:
            print(f"{prefix}HADDOCK3 failed -- last output from {log_path}:")
            with open(log_path) as log_f:
                print("".join(log_f.readlines()[-30:]))
        raise
    finally:
        stop.set()
        thread.join()


def print_table(rows, columns, title=None):
    if not rows:
        print("(no rows)")
        return
    get = lambda r, key: key(r) if callable(key) else str(r[key])
    widths = {label: max(len(label), *(len(get(r, key)) for r in rows)) for label, key in columns}
    if title:
        print(f"\n== {title} ==")
    header_line = "  ".join(label.ljust(widths[label]) for label, _ in columns)
    print(header_line)
    print("-" * len(header_line))
    for r in rows:
        print("  ".join(get(r, key).ljust(widths[label]) for label, key in columns))


def read_final_capri_rows(haddock_run_dir):
    caprieval_dirs = sorted(
        glob.glob(os.path.join(haddock_run_dir, "[0-9]*_caprieval")),
        key=lambda p: int(os.path.basename(p).split("_")[0]),
    )
    if not caprieval_dirs:
        return None, None
    final_dir = caprieval_dirs[-1]
    with open(os.path.join(final_dir, "capri_clt.tsv")) as f:
        lines = [line for line in f if line.strip() and not line.startswith("#")]
    header = lines[0].strip().split("\t")
    rows = [dict(zip(header, line.strip().split("\t"))) for line in lines[1:]]
    rows.sort(key=lambda r: int(r["caprieval_rank"]))
    return final_dir, rows


CAPRI_COLUMNS = [
    ("rank", "caprieval_rank"), ("cluster", "cluster_id"), ("n", "n"),
    ("score", "score"), ("dockq", "dockq"),
    ("irmsd", "irmsd"), ("fnat", "fnat"), ("lrmsd", "lrmsd"),
]


def print_capri_summary(haddock_run_dir):
    final_dir, rows = read_final_capri_rows(haddock_run_dir)
    if rows is None:
        print("No caprieval output found.")
        return
    print_table(rows, CAPRI_COLUMNS, title=f"Final CAPRI cluster results ({os.path.basename(final_dir)})")


# Locked output schema (PROTOCOL_LOCK.md section 6): the one CSV both Ryan and Tyrone emit.
LOCKED_COLUMNS = ["molecule", "mean_best_10", "best_haddock", "sd_best_10",
                  "dockq_vs_9dwv", "cluster_id", "n_models", "protocol_version"]


def locked_stats(haddock_run_dir):
    """Locked-schema per-run stats from the final caprieval's single-model table
    (capri_ss.tsv): mean/best/sd over the 10 best-scoring models, plus the best
    model's cluster id and the number of models scored. dockq_vs_9dwv is left
    blank here (self-referential in 08 -- 09 fills the real value against 9DWV).
    Returns a dict, or None if no models were scored."""
    final_dir, _ = read_final_capri_rows(haddock_run_dir)
    if final_dir is None:
        return None
    ss_path = os.path.join(final_dir, "capri_ss.tsv")
    if not os.path.exists(ss_path):
        return None
    with open(ss_path) as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    if len(lines) < 2:
        return None
    header = lines[0].strip().split("\t")
    recs = [dict(zip(header, l.strip().split("\t"))) for l in lines[1:]]
    recs.sort(key=lambda r: float(r["score"]))
    best10 = [float(r["score"]) for r in recs[:10]]
    return {
        "mean_best_10": round(statistics.fmean(best10), 3),
        "best_haddock": round(best10[0], 3),
        "sd_best_10": round(statistics.stdev(best10), 3) if len(best10) > 1 else 0.0,
        "dockq_vs_9dwv": "",
        "cluster_id": recs[0].get("cluster_id", "-"),
        "n_models": len(recs),
    }


def candidate_already_completed(candidate_name):
    """Detect whether candidate_name has already been FULLY docked in the
    current run tree (RUN_DIR_BASE) -- the whole locked workflow finished and
    the final caprieval (9_caprieval, the same step 07/00_validate_08 read)
    wrote a capri_ss.tsv with at least one scored model.

    This is the real meaning of "already run", and why it beats a ledger-name
    check: a candidate that crashed partway (e.g. the old iCloud-race caprieval
    failure) or was interrupted by hand leaves an earlier caprieval and/or a
    "-" ledger row but NO valid final result. Those must re-run, not be skipped
    forever. So main() skips a candidate iff this returns True and re-runs
    everything else (never-run, failed, or interrupted).
    """
    ss_path = os.path.join(RUN_DIR_BASE, candidate_name, "run1",
                           "9_caprieval", "capri_ss.tsv")  # 9 = final STEP_PLAN step; matches 07/00
    if not os.path.exists(ss_path):
        return False
    try:
        with open(ss_path) as f:
            data_lines = [l for l in f if l.strip() and not l.startswith("#")]
    except OSError:
        return False
    return len(data_lines) >= 2  # header + at least one scored model


def pick_candidate_names():
    """Returns the pool of candidates eligible for this script: whatever's
    currently in 06's progress ledger with a real (non "-") result. Not
    filtered/ranked by 06's dockq -- see module docstring for why."""
    if CANDIDATE_NAME:
        return [CANDIDATE_NAME]
    if not os.path.exists(SIX_PROGRESS_CSV):
        sys.exit(f"{SIX_PROGRESS_CSV} not found -- run 06 first, or set CANDIDATE_NAME directly.")
    with open(SIX_PROGRESS_CSV, newline="") as f:
        names = [r["name"] for r in csv.DictReader(f) if r["dockq"] != "-"]
    if not names:
        sys.exit(f"No candidates with a real result in {SIX_PROGRESS_CSV} yet -- "
                  "run 06 first, or set CANDIDATE_NAME directly.")
    print(f"CANDIDATE_NAME not set -- {len(names)} candidate(s) eligible from {SIX_PROGRESS_CSV} "
          f"(completed 06 without failing): {', '.join(names)}")
    return names


def load_existing_rows(path=RESULTS_CSV):
    """Load whatever's already in `path` (default RESULTS_CSV) so a run only
    adds new molecules -- never re-does a 45 min - 1.5 hr run it already has
    a result for."""
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows or "mean_best_10" not in rows[0]:
        return []  # pre-lock (old-schema) ledger -- ignore it so the locked run starts clean
    loaded = []
    for r in rows:
        if r["mean_best_10"] == "-":
            loaded.append((r["molecule"], None))
        else:
            loaded.append((r["molecule"], {
                "mean_best_10": r["mean_best_10"], "best_haddock": r["best_haddock"],
                "sd_best_10": r["sd_best_10"], "dockq_vs_9dwv": r.get("dockq_vs_9dwv", ""),
                "cluster_id": r["cluster_id"], "n_models": r["n_models"],
            }))
    return loaded


def write_results_csv(all_rows, path=RESULTS_CSV):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(LOCKED_COLUMNS)
        for candidate_name, r in all_rows:
            if r is None:
                writer.writerow([candidate_name, "-", "-", "-", "-", "-", "-", PROTOCOL_VERSION])
            else:
                writer.writerow([candidate_name, r["mean_best_10"], r["best_haddock"], r["sd_best_10"],
                                  r.get("dockq_vs_9dwv", ""), r["cluster_id"], r["n_models"], PROTOCOL_VERSION])
    print(f"Wrote {path}")


def extract_ligand_pdb(candidate_name, out_path):
    """Pull just the ligand (resn LIG, chain C) HETATM lines out of 04's
    CRBN+ligand complex into a standalone PDB."""
    complex_path = os.path.join(VINA_OUT_DIR, candidate_name, "CRBN_candidate_complex.pdb")
    if not os.path.exists(complex_path):
        sys.exit(f"{complex_path} not found -- run 04 for this candidate first "
                  "(MANUAL_CANDIDATE for a one-off compound like a positive control).")
    with open(complex_path) as f:
        lig_lines = [l for l in f if l.startswith("HETATM") and l[17:20] == "LIG"]
    if not lig_lines:
        sys.exit(f"No LIG HETATM records found in {complex_path}.")
    with open(out_path, "w") as f:
        f.writelines(lig_lines)
        f.write("END\n")


def _sanitize_atom_names(content):
    """Mirrors haddock.libs.libligand._sanitize_atom_names: PRODRG can emit
    atom type names containing colons (e.g. "HT:A"), which CNS rejects."""
    lines = []
    for line in content.splitlines(keepends=True):
        if line.lstrip().startswith("!"):
            lines.append(line)
        else:
            lines.append(line.replace(":", ""))
    return "".join(lines)


def _remove_nbonds(par_content):
    """Mirrors haddock.libs.libligand._remove_nbonds: PRODRG's own NBONds
    block can conflict with HADDOCK's internal nonbonded parameters."""
    return re.sub(r"(?s)NBONds.*?END", "", par_content)


def run_prodrg_for_ligand(ligand_pdb, out_dir):
    """Run PRODRG on a ligand PDB and return (top_path, param_path,
    fixed_pdb_path) -- CNS topology, parameters, and a coordinate PDB whose
    atom names actually match that topology (chain C, resnum 900, matching
    this project's convention), all written into out_dir.

    Unlike haddock.libs.libligand.run_prodrg (which discards PRODRG's own
    renamed-atom PDB, DRGFIN.PDB), this keeps it -- see this file's
    docstring for why that renaming matters."""
    prodrg_exec, prodrg_param = get_prodrg_exec()
    if prodrg_exec is None:
        sys.exit("PRODRG binary not found for this platform (haddock.libs.libutil.get_prodrg_exec "
                  "returned None) -- this script requires the bundled PRODRG binary to parametrize "
                  "the ligand.")

    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copy(prodrg_param, os.path.join(tmpdir, os.path.basename(prodrg_param)))
        shutil.copy(ligand_pdb, os.path.join(tmpdir, os.path.basename(ligand_pdb)))

        result = subprocess.run(
            [str(prodrg_exec), os.path.basename(ligand_pdb), os.path.basename(prodrg_param), "PDBELEM"],
            cwd=tmpdir, capture_output=True, text=True,
        )
        if result.returncode != 0:
            sys.exit(f"PRODRG failed:\n{result.stderr}")

        tmp_top = os.path.join(tmpdir, "DRGCNS.TOP")
        tmp_par = os.path.join(tmpdir, "DRGCNS.PAR")
        tmp_pdb = os.path.join(tmpdir, "DRGFIN.PDB")
        if not (os.path.exists(tmp_top) and os.path.exists(tmp_par) and os.path.exists(tmp_pdb)):
            sys.exit(f"PRODRG finished but expected output files are missing in {tmpdir}: "
                      f"{sorted(os.listdir(tmpdir))}")

        stem = os.path.splitext(os.path.basename(ligand_pdb))[0]
        top_path = os.path.join(out_dir, f"{stem}_prodrg.top")
        par_path = os.path.join(out_dir, f"{stem}_prodrg.param")
        fixed_pdb_path = os.path.join(out_dir, f"{stem}_prodrg_fixed.pdb")

        with open(tmp_top) as f:
            top_content = _sanitize_atom_names(f.read())
        with open(tmp_par) as f:
            par_content = _sanitize_atom_names(_remove_nbonds(f.read()))
        with open(tmp_pdb) as f:
            pdb_lines = f.readlines()

        with open(top_path, "w") as f:
            f.write(top_content)
        with open(par_path, "w") as f:
            f.write(par_content)
        with open(fixed_pdb_path, "w") as f:
            for line in pdb_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    line = line[:21] + "C" + " 900" + line[26:]
                f.write(line)
            if not pdb_lines[-1].startswith("END"):
                f.write("END\n")

    return top_path, par_path, fixed_pdb_path


def run_candidate(candidate_name, label):
    """Run the ligand-inclusive 3-body HADDOCK3 routine for one candidate under
    the LOCKED protocol (see PROTOCOL_LOCK.md): fixed restraints (ambig_FIXED.tbl,
    identical for every candidate) and a pinned iniseed. Returns the top
    (rank-1, best-score) cluster row, or None if setup/HADDOCK3 fails or no
    cluster was found."""
    candidate_run_dir = os.path.join(RUN_DIR_BASE, candidate_name)
    os.makedirs(candidate_run_dir, exist_ok=True)
    print(f"[{label}] step 0/{len(STEP_PLAN)}: preparing ligand topology (restraints are locked)")

    # --- Ligand: extract 04's Vina pose (the LOCKED starting-pose recipe), parametrize via PRODRG ---
    # Only the ligand varies between candidates now; the CRBN/PPIL4 restraints are the
    # same committed file (ambig_FIXED.tbl, LOCK #2) for all of them, so there is no
    # per-candidate restraint derivation here anymore -- that is what removes the Vina noise.
    ligand_raw_pdb = os.path.join(candidate_run_dir, "ligand.pdb")
    extract_ligand_pdb(candidate_name, ligand_raw_pdb)
    ligand_top, ligand_param, ligand_pdb = run_prodrg_for_ligand(ligand_raw_pdb, candidate_run_dir)

    ambig_tbl = AMBIG_FIXED  # LOCK #2: one committed restraint file for every candidate/control

    # --- HADDOCK3 config: 3 molecules, ligand_top_fname/ligand_param_fname repeated on
    # every module that runs its own CNS minimization (topoaa, rigidbody, flexref, emref) ---
    haddock_run_dir = os.path.join(candidate_run_dir, "run1")
    cfg_path = os.path.join(candidate_run_dir, "haddock3_ternary_with_ligand.toml")
    cfg = f"""
run_dir = "{haddock_run_dir}"
ncores = {NCORES}

molecules = [
    "{CRBN_RECEPTOR_ONLY_PDB}",
    "{ligand_pdb}",
    "{PPIL4_PDB}"
]

[topoaa]
ligand_top_fname = "{ligand_top}"
ligand_param_fname = "{ligand_param}"

[rigidbody]
ambig_fname = "{ambig_tbl}"
ligand_top_fname = "{ligand_top}"
ligand_param_fname = "{ligand_param}"
sampling = 1000
iniseed = {INISEED}

[caprieval]

[seletop]
select = 200

[flexref]
ambig_fname = "{ambig_tbl}"
ligand_top_fname = "{ligand_top}"
ligand_param_fname = "{ligand_param}"

[caprieval]

[emref]
ambig_fname = "{ambig_tbl}"
ligand_top_fname = "{ligand_top}"
ligand_param_fname = "{ligand_param}"

[caprieval]

[clustfcc]

[caprieval]
"""
    with open(cfg_path, "w") as f:
        f.write(cfg.strip() + "\n")

    if os.path.exists(haddock_run_dir):
        shutil.rmtree(haddock_run_dir)
    print(f"[{label}] starting HADDOCK3 3-body (CRBN+ligand+PPIL4) run "
          f"(rough estimate: 45 min - 1.5 hr, similar to 06)")
    run_with_heartbeat(["haddock3", cfg_path], run_dir=haddock_run_dir, step_plan=STEP_PLAN, label=label)

    print_capri_summary(haddock_run_dir)
    return locked_stats(haddock_run_dir)


def write_finalists(all_rows):
    """Best TOP_N_KEEP of all_rows by HADDOCK score (lower/more negative =
    better) -- not dockq (self-referential, see module docstring). Written
    every time, not gated on the whole pool finishing -- 08 is expensive
    enough that seeing an up-to-date ranking as candidates trickle in is
    more useful than waiting for a batch to fully complete."""
    finished = [(n, r) for n, r in all_rows if r is not None]
    failed = len(all_rows) - len(finished)
    if failed:
        print(f"{failed}/{len(all_rows)} candidate(s) had no HADDOCK3 result (failed/skipped) -- "
              "excluded from the finalist pool.")
    finished.sort(key=lambda pair: float(pair[1]["mean_best_10"]))
    finalists = finished[:TOP_N_KEEP]
    write_results_csv(finalists, path=FINALISTS_CSV)
    print(f"Finalists (best {len(finalists)} of {len(finished)} completed, by mean_best_10): "
          f"{', '.join(n for n, _ in finalists)}")


def run_named_set(names, csv_path, kind):
    """Dock a fixed, named set of molecules (the controls, or the cross-check
    set) through the SAME locked protocol as the candidates, writing them to
    their OWN csv_path -- separate from the candidate finalists. A molecule
    whose 04 Vina pose doesn't exist yet is skipped with a message (08 can't
    make the pose -- no Vina in the haddock venv -- so run 04 for it first).
    Resumes: a molecule that already has a real result is not re-run."""
    result = {name: row for name, row in load_existing_rows(csv_path)}
    to_run = [n for n in names if result.get(n) is None]  # includes prior skips/failures -> retried
    if not to_run:
        print(f"All {len(names)} {kind}s already have results in {csv_path} -- skipping {kind}s.")
        return
    print(f"\n== Docking {len(to_run)} {kind}(s) into {csv_path}: {', '.join(to_run)} ==")
    for name in to_run:
        complex_pdb = os.path.join(VINA_OUT_DIR, name, "CRBN_candidate_complex.pdb")
        if not os.path.exists(complex_pdb):
            print(f"[{kind} {name}] no 04 Vina pose ({os.path.relpath(complex_pdb, SCRIPT_DIR)} missing) -- "
                  f"run 04 for it first (its SMILES is in 04's CONTROLS/CROSSCHECK list). Skipping.")
            result[name] = None
        else:
            try:
                result[name] = run_candidate(name, f"{kind} {name}")
            except subprocess.CalledProcessError:
                print(f"[{kind} {name}] HADDOCK3 failed -- skipping.")
                result[name] = None
        write_results_csv([(n, result.get(n)) for n in names], path=csv_path)
    ok = [n for n in names if result.get(n) is not None]
    print(f"{kind.capitalize()}s: {len(ok)}/{len(names)} docked -> {csv_path} ({', '.join(ok) or 'none yet'})")


def main():
    os.makedirs(RUN_DIR_BASE, exist_ok=True)

    # Locked protocol: restraints come from the committed ambig_FIXED.tbl (identical for
    # every candidate/control) and PPIL4_chainB.pdb is a committed input -- see
    # PROTOCOL_LOCK.md. Verify the locked bundle is present before docking anything.
    for artifact in (AMBIG_FIXED, PPIL4_PDB, CRBN_RECEPTOR_ONLY_PDB):
        if not os.path.exists(artifact):
            sys.exit(f"Locked artifact missing: {artifact} -- rebuild the locked bundle (see PROTOCOL_LOCK.md).")

    # Controls first (funnel mode only): dock the calibration set through the same
    # locked protocol into their OWN CSV (kept out of the candidate finalists).
    if RUN_CONTROLS and not CANDIDATE_NAME:
        run_named_set(CONTROLS, CONTROLS_CSV, "control")
    # Then the cross-check set (Tyrone's molecules through this protocol), also
    # into its own CSV -- same reason: it isn't part of this funnel's ranking.
    if RUN_CROSSCHECK and not CANDIDATE_NAME:
        run_named_set(CROSSCHECK, CROSSCHECK_CSV, "crosscheck")

    candidate_names = pick_candidate_names()

    # --- Resume by DETECTING finished docks, not by trusting ledger names ---
    # results_by_name maps candidate -> its result dict (real) or None (failed/
    # placeholder). A candidate is "already run" only if it has a real result,
    # i.e. its final caprieval actually completed. Interrupted/failed runs (a
    # "-" ledger row and maybe a partial run dir) fall through and re-run.
    results_by_name = {name: r for name, r in load_existing_rows()}

    # Self-heal: a candidate can be fully docked on disk yet missing from the
    # ledger (e.g. the ledger was reset, or the run tree was moved -- as it just
    # was, off the iCloud Desktop). Recompute its result from its run dir so it
    # still ranks and isn't needlessly re-docked.
    for name in candidate_names:
        if results_by_name.get(name) is None and candidate_already_completed(name):
            stats = locked_stats(os.path.join(RUN_DIR_BASE, name, "run1"))
            if stats is not None:
                results_by_name[name] = stats
                print(f"Recovered {name}'s completed result from its run dir (ledger was missing it).")

    def ledger_rows():
        # Every ledger entry (stable order), then any pool candidate not yet present.
        names = list(results_by_name) + [n for n in candidate_names if n not in results_by_name]
        return [(n, results_by_name.get(n)) for n in names]

    done_names = [n for n in candidate_names if results_by_name.get(n) is not None]
    new_candidates = [n for n in candidate_names if results_by_name.get(n) is None]
    if done_names:
        print(f"Skipping {len(done_names)} candidate(s) already fully docked in {RUN_DIR_BASE}: "
              f"{', '.join(done_names)}")
    if not new_candidates:
        print("Nothing new to run this session -- every candidate already has a completed dock.")
        if not CANDIDATE_NAME:
            write_finalists(ledger_rows())
        return
    print(f"Running the ligand-inclusive HADDOCK3 routine on {len(new_candidates)} not-yet-completed "
          f"candidate(s): {', '.join(new_candidates)}")

    candidates_loop_start = time.time()
    for i, candidate_name in enumerate(new_candidates, 1):
        remaining = len(new_candidates) - i
        elapsed_so_far = time.time() - candidates_loop_start
        avg_per_candidate = elapsed_so_far / (i - 1) if i > 1 else None
        eta_str = (f"~{avg_per_candidate * remaining / 60:.0f} min remaining for the run"
                   if avg_per_candidate else "remaining time unknown until candidate 1 finishes")
        label = f"{candidate_name}, {i}/{len(new_candidates)}, {remaining} left"
        print(f"\n[{label}] ({elapsed_so_far / 60:.1f} min elapsed this run, {eta_str} | "
              f"{(time.time() - SCRIPT_START_TIME) / 60:.1f} min total script time)")

        try:
            top_row = run_candidate(candidate_name, label)
        except subprocess.CalledProcessError:
            print(f"[{label}] HADDOCK3 failed for this candidate -- skipping it and continuing with the rest.")
            top_row = None

        results_by_name[candidate_name] = top_row  # update in place -- re-running a failed candidate never duplicates its row
        if not CANDIDATE_NAME:
            write_results_csv(ledger_rows())  # unsorted, for resilience
            write_finalists(ledger_rows())

    if CANDIDATE_NAME:
        print(f"\nCANDIDATE_NAME set -- skipping {RESULTS_CSV}/{FINALISTS_CSV} (those are the funnel's "
              f"shared ranked output; a one-off candidate isn't part of that ranking). Its docking output "
              f"is still under {RUN_DIR_BASE}/{CANDIDATE_NAME}/.")

    total = time.time() - SCRIPT_START_TIME
    print(f"\nTotal script runtime: {total:.0f}s ({total / 60:.1f} min)")


if __name__ == "__main__":
    main()
