"""
Extract the representative 3D structure (the rank-1, best-scoring model)
for every candidate in 08's final ranked output, so they can be opened
directly in PyMOL for visual inspection/comparison.

08_haddock3_ternary_with_ligand_(Ryan).py is the pipeline's actual final
ranking stage now (see chat discussion -- 06's dockq turned out to be
structurally insensitive to candidate chemistry, since 05/06 never put
the ligand in the CNS topology). Its models already contain the real
ligand natively as chain C, docked as part of the actual 3-body CNS force
field computation -- not bolted on afterward. That makes this script much
simpler than its old (06-based) version: no PyMOL alignment/merge hack is
needed anymore, since there's nothing missing to merge in. This just
resolves each candidate's rank-1 model (from docking_tmp/
haddock3_ternary_with_ligand_run/<candidate>/run1/9_caprieval/capri_ss.tsv,
same TSV lookup logic as 08/09) and copies it to ternary_structures/ with
a clear name.

Run with any Python that can read 08_final_ternary_with_ligand_results_(Ryan).csv
(no haddock3/rdkit/vina/pymol dependency needed):
    python3 "07_extract_top_structures_(Ryan).py"
"""
import csv
import gzip
import os
import shutil
import sys
import time

SCRIPT_START_TIME = time.time()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

RESULTS_CSV = os.path.join(SCRIPT_DIR, "08_final_ternary_with_ligand_results_(Ryan).csv")
RUN_DIR_BASE = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "ternary_structures")

# How many of 08's top-ranked (by HADDOCK score) candidates to extract a
# structure for. None = every candidate in RESULTS_CSV. RESULTS_CSV is
# already sorted best-score-first.
TOP_N = None


def find_top_model_path(candidate_name):
    """Return the absolute path to candidate_name's rank-1 model file
    (still gzipped if HADDOCK3 left it that way), or None if not found."""
    caprieval_dir = os.path.join(RUN_DIR_BASE, candidate_name, "run1", "9_caprieval")
    tsv_path = os.path.join(caprieval_dir, "capri_ss.tsv")
    if not os.path.exists(tsv_path):
        return None

    with open(tsv_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    top = next((r for r in rows if r["caprieval_rank"] == "1"), None)
    if top is None:
        return None

    rel_path = top["model"]  # e.g. "../6_emref/emref_3.pdb"
    resolved = os.path.normpath(os.path.join(caprieval_dir, rel_path))
    if os.path.exists(resolved):
        return resolved
    if os.path.exists(resolved + ".gz"):
        return resolved + ".gz"
    return None


def extract_structure(candidate_name):
    """Decompress (if needed) candidate_name's rank-1 model into OUTPUT_DIR.
    Returns the written path, or None if the model couldn't be found."""
    src_path = find_top_model_path(candidate_name)
    if src_path is None:
        print(f"[{candidate_name}] no rank-1 model found (has 08 been run for this candidate?) -- skipping.")
        return None

    out_path = os.path.join(OUTPUT_DIR, f"07_best_model_{candidate_name}_(Ryan).pdb")
    if src_path.endswith(".gz"):
        with gzip.open(src_path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copyfile(src_path, out_path)
    print(f"[{candidate_name}] wrote {out_path} (from {os.path.relpath(src_path, SCRIPT_DIR)})")
    return out_path


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(RESULTS_CSV):
        sys.exit(f"{RESULTS_CSV} not found -- run 08 first.")
    with open(RESULTS_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    candidates = [r["name"] for r in (rows[:TOP_N] if TOP_N else rows)]
    print(f"Extracting structures for {len(candidates)} candidate(s) from {RESULTS_CSV}: "
          f"{', '.join(candidates)}")

    written = []
    for i, candidate_name in enumerate(candidates, 1):
        print(f"[{i}/{len(candidates)}] {candidate_name}")
        out_path = extract_structure(candidate_name)
        if out_path:
            written.append(out_path)

    print(f"\nWrote {len(written)} file(s) to {OUTPUT_DIR}")

    total = time.time() - SCRIPT_START_TIME
    print(f"Total script runtime: {total:.0f}s ({total / 60:.1f} min)")


if __name__ == "__main__":
    main()
