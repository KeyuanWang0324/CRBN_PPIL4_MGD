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
SIX_PROGRESS_CSV = os.path.join(SCRIPT_DIR, "06_complete_run_progress_(Ryan).csv")
# This script's own running ledger -- every candidate attempted so far,
# win or lose. Read back in on every run so a later session only adds new
# candidates (06's progress ledger can grow over time as 06 keeps running;
# nothing here requires 06 to be "finished" first).
RESULTS_CSV = os.path.join(SCRIPT_DIR, "08_ternary_with_ligand_progress_(Ryan).csv")
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
RUN_DIR_BASE = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run")
PPIL4_PDB = os.path.join(RUN_DIR_BASE, "PPIL4_chainB.pdb")

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


def load_existing_rows():
    """Load whatever's already in RESULTS_CSV so this run only adds new
    candidates -- never re-does a 45 min - 1.5 hr candidate it already has
    a result for."""
    if not os.path.exists(RESULTS_CSV):
        return []
    with open(RESULTS_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    loaded = []
    for r in rows:
        if r["dockq"] == "-":
            loaded.append((r["name"], None))
        else:
            loaded.append((r["name"], {
                "caprieval_rank": r["cluster_rank"], "cluster_id": r["cluster_id"], "n": r["n"],
                "score": r["score"], "dockq": r["dockq"], "irmsd": r["irmsd"],
                "fnat": r["fnat"], "lrmsd": r["lrmsd"],
            }))
    return loaded


def write_results_csv(all_rows, path=RESULTS_CSV):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "cluster_rank", "cluster_id", "n", "score", "dockq", "irmsd", "fnat", "lrmsd"])
        for candidate_name, r in all_rows:
            if r is None:
                writer.writerow([candidate_name, "-", "-", "-", "-", "-", "-", "-", "-"])
            else:
                writer.writerow([candidate_name, r["caprieval_rank"], r["cluster_id"], r["n"],
                                  r["score"], r["dockq"], r["irmsd"], r["fnat"], r["lrmsd"]])
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


def run_candidate(candidate_name, ppil4_actpass, label):
    """Run the ligand-inclusive 3-body HADDOCK3 routine for one candidate.
    Returns the top (rank-1, best-score) cluster row, or None if setup/
    HADDOCK3 fails or no cluster was found."""
    candidate_run_dir = os.path.join(RUN_DIR_BASE, candidate_name)
    os.makedirs(candidate_run_dir, exist_ok=True)
    print(f"[{label}] step 0/{len(STEP_PLAN)}: preparing ligand topology/restraints")

    # --- CRBN-side setup (from 04's Vina contact residues) ---
    contacts_path = os.path.join(VINA_OUT_DIR, candidate_name, "crbn_contacts.txt")
    if not os.path.exists(contacts_path):
        print(f"[{label}] {contacts_path} not found -- run 04 for this candidate first. Skipping.")
        return None
    with open(contacts_path) as f:
        crbn_active = [int(x) for x in f.readline().split()]
    crbn_active_csv = ",".join(str(r) for r in crbn_active)
    crbn_passive_out = subprocess.run(
        ["haddock3-restraints", "passive_from_active", CRBN_RECEPTOR_ONLY_PDB, crbn_active_csv, "-c", "A"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    crbn_passive = [int(x) for x in crbn_passive_out.split()] if crbn_passive_out else []
    crbn_actpass = os.path.join(candidate_run_dir, "crbn_actpass.txt")
    write_actpass_file(crbn_active, crbn_passive, crbn_actpass)

    # --- Ligand: extract, parametrize via PRODRG, and mark the whole thing "active" ---
    ligand_raw_pdb = os.path.join(candidate_run_dir, "ligand.pdb")
    extract_ligand_pdb(candidate_name, ligand_raw_pdb)
    ligand_top, ligand_param, ligand_pdb = run_prodrg_for_ligand(ligand_raw_pdb, candidate_run_dir)
    ligand_actpass = os.path.join(candidate_run_dir, "ligand_actpass.txt")
    write_actpass_file([900], [], ligand_actpass)  # the whole ligand is one CNS residue (resnum 900)

    # --- Restraints: ligand<->CRBN and ligand<->PPIL4, concatenated into one ambig.tbl ---
    crbn_ligand_tbl = subprocess.run(
        ["haddock3-restraints", "active_passive_to_ambig", crbn_actpass, ligand_actpass,
         "--segid-one", "A", "--segid-two", "C"],
        check=True, capture_output=True, text=True,
    ).stdout
    ppil4_ligand_tbl = subprocess.run(
        ["haddock3-restraints", "active_passive_to_ambig", ppil4_actpass, ligand_actpass,
         "--segid-one", "B", "--segid-two", "C"],
        check=True, capture_output=True, text=True,
    ).stdout
    ambig_tbl = os.path.join(candidate_run_dir, "ambig.tbl")
    with open(ambig_tbl, "w") as f:
        f.write(crbn_ligand_tbl)
        f.write(ppil4_ligand_tbl)

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
    _, rows = read_final_capri_rows(haddock_run_dir)
    return rows[0] if rows else None


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
    finished.sort(key=lambda pair: float(pair[1]["score"]))
    finalists = finished[:TOP_N_KEEP]
    write_results_csv(finalists, path=FINALISTS_CSV)
    print(f"Finalists (best {len(finalists)} of {len(finished)} completed, by HADDOCK score): "
          f"{', '.join(n for n, _ in finalists)}")


def main():
    os.makedirs(RUN_DIR_BASE, exist_ok=True)
    candidate_names = pick_candidate_names()

    all_rows = load_existing_rows()
    existing_names = {name for name, _ in all_rows}
    if existing_names:
        print(f"{RESULTS_CSV} already has {len(existing_names)} candidate(s) -- keeping them, only adding new ones.")

    new_candidates = [n for n in candidate_names if n not in existing_names]
    already_done = [n for n in candidate_names if n in existing_names]
    if already_done:
        print(f"Skipping {len(already_done)} candidate(s) already in {RESULTS_CSV}: {', '.join(already_done)}")
    if not new_candidates:
        print("Nothing new to run this session.")
        if not CANDIDATE_NAME:
            write_finalists(all_rows)
        return
    print(f"Running the ligand-inclusive HADDOCK3 routine on {len(new_candidates)} candidate(s): "
          f"{', '.join(new_candidates)}")

    # PPIL4-side restraints don't depend on the candidate -- compute once,
    # not once per candidate. Real RRM-domain interface residues from
    # FPFT-2216 (PDB 9DWV) -- same as 05/06's ppil4_pocket_residues().
    with open(PPIL4_PDB, "w") as out:
        run(["pdb_chain", "-B", PPIL4_SOURCE_PDB], stdout=out)
    ppil4_active = [249, 250, 273, 275, 276, 277, 278, 279]
    ppil4_active_csv = ",".join(str(r) for r in ppil4_active)
    ppil4_passive_out = subprocess.run(
        ["haddock3-restraints", "passive_from_active", PPIL4_PDB, ppil4_active_csv, "-c", "B"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    ppil4_passive = [int(x) for x in ppil4_passive_out.split()] if ppil4_passive_out else []
    ppil4_actpass = os.path.join(RUN_DIR_BASE, "ppil4_actpass.txt")
    write_actpass_file(ppil4_active, ppil4_passive, ppil4_actpass)

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
            top_row = run_candidate(candidate_name, ppil4_actpass, label)
        except subprocess.CalledProcessError:
            print(f"[{label}] HADDOCK3 failed for this candidate -- skipping it and continuing with the rest.")
            top_row = None

        all_rows.append((candidate_name, top_row))
        if not CANDIDATE_NAME:
            write_results_csv(all_rows)  # unsorted, for resilience
            write_finalists(all_rows)

    if CANDIDATE_NAME:
        print(f"\nCANDIDATE_NAME set -- skipping {RESULTS_CSV}/{FINALISTS_CSV} (those are the funnel's "
              f"shared ranked output; a one-off candidate isn't part of that ranking). Its docking output "
              f"is still under {RUN_DIR_BASE}/{CANDIDATE_NAME}/.")

    total = time.time() - SCRIPT_START_TIME
    print(f"\nTotal script runtime: {total:.0f}s ({total / 60:.1f} min)")


if __name__ == "__main__":
    main()
