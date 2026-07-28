"""
Population-scale version of 00_validate_docking_interface_(Ryan).py's
agree() check, run against 05's actual HADDOCK3 ternary-docking models
instead of 04's Vina poses (see 00_validate_docking_interface_04_(Ryan).py
for that version) -- this is the more direct comparison, since 05's
models are real protein-protein (CRBN chain A vs PPIL4 chain B) docking
results, the same kind of structure REAL_CRBN/REAL_PPIL4 and the original
cand_5 example (9/10 CRBN, 0/8 PPIL4) were measured from -- unlike 04,
which only checks ligand-to-protein contacts.

REAL_CRBN/REAL_PPIL4 are the same ground-truth residue sets from
session4_handson.md (Part 2a) -- pasted verbatim, not recomputed, not
invented.

YOUR_CRBN/YOUR_PPIL4 per candidate are computed the same way
session4_handson.md says the original REAL_CRBN/REAL_PPIL4/cand_5 numbers
were: heavy-atom (no hydrogens) contacts within 5 A across chain A
(CRBN) and chain B (PPIL4), taken from that candidate's own rank-1 model
(same TSV lookup logic as 07_extract_top_structures/09) out of 05's
run1/5_caprieval/capri_ss.tsv -- 05's actual HADDOCK3 output, not
invented, not re-docked.

Run with the SYSTEM python:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "00_validate_docking_interface_05_(Ryan).py"
"""
import csv
import gzip
import os
import statistics

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

RESULTS_CSV = os.path.join(SCRIPT_DIR, "05_ternary_docking_scores_for_06_(Ryan).csv")
RUN_DIR_BASE = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_novel_candidate_run")

# Real interface residues from the FPFT-2216 / PDB 9DWV structure (see
# session4_handson.md's REAL_CRBN/REAL_PPIL4).
REAL_CRBN = {351, 353, 355, 357, 372, 373, 386, 388, 397, 400}
REAL_PPIL4 = {249, 250, 273, 275, 276, 277, 278, 279}

CUTOFF = 5.0  # heavy-atom contact distance, matches how REAL_CRBN/REAL_PPIL4 were measured


def agree(your_set, real_set):
    found = your_set & real_set
    verdict = "FOUND THE SPOT" if len(found) >= len(real_set) / 2 else "WRONG SPOT"
    return found, verdict


def find_top_model_path(candidate_name):
    """Same resolution logic as 07_extract_top_structures/09: the rank-1
    model from 05's own (lite) HADDOCK3 run for this candidate."""
    caprieval_dir = os.path.join(RUN_DIR_BASE, candidate_name, "run1", "5_caprieval")
    tsv_path = os.path.join(caprieval_dir, "capri_ss.tsv")
    if not os.path.exists(tsv_path):
        return None
    with open(tsv_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    top = next((r for r in rows if r["caprieval_rank"] == "1"), None)
    if top is None:
        return None
    resolved = os.path.normpath(os.path.join(caprieval_dir, top["model"]))
    if os.path.exists(resolved):
        return resolved
    if os.path.exists(resolved + ".gz"):
        return resolved + ".gz"
    return None


def read_pdb_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return f.readlines()


def chain_heavy_atoms(lines, chain):
    """(resnum, x, y, z) for heavy (non-hydrogen) ATOM records on one chain."""
    atoms = []
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        element = line[76:78].strip()
        if element == "H":
            continue
        resnum = int(line[22:26])
        atoms.append((resnum, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def chain_contacts(atoms_a, atoms_b, cutoff=CUTOFF):
    """Residue numbers on side A within cutoff of any heavy atom on side B."""
    contacts = set()
    for resnum, x, y, z in atoms_a:
        for _, bx, by, bz in atoms_b:
            if (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2 <= cutoff ** 2:
                contacts.add(resnum)
                break
    return contacts


def main():
    with open(RESULTS_CSV, newline="") as f:
        rows = list(csv.DictReader(f))

    crbn_overlaps, ppil4_overlaps = [], []
    crbn_found, ppil4_found = 0, 0
    missing = 0

    for row in rows:
        name = row["name"]
        model_path = find_top_model_path(name)
        if model_path is None:
            missing += 1
            continue

        lines = read_pdb_lines(model_path)
        crbn_atoms = chain_heavy_atoms(lines, "A")
        ppil4_atoms = chain_heavy_atoms(lines, "B")

        your_crbn = chain_contacts(crbn_atoms, ppil4_atoms)
        your_ppil4 = chain_contacts(ppil4_atoms, crbn_atoms)

        found, verdict = agree(your_crbn, REAL_CRBN)
        crbn_overlaps.append(len(found))
        if verdict == "FOUND THE SPOT":
            crbn_found += 1

        found, verdict = agree(your_ppil4, REAL_PPIL4)
        ppil4_overlaps.append(len(found))
        if verdict == "FOUND THE SPOT":
            ppil4_found += 1

    n = len(rows)
    print(f"05's ternary docking CSV: {n} candidates (missing rank-1 model: {missing})\n")

    print(f"CRBN  (real interface has {len(REAL_CRBN)} residues, threshold {len(REAL_CRBN) / 2:.0f}+ to FOUND THE SPOT)")
    print(f"  FOUND THE SPOT: {crbn_found}/{len(crbn_overlaps)} ({crbn_found / len(crbn_overlaps):.0%})")
    print(f"  overlap count -- mean: {statistics.mean(crbn_overlaps):.2f}, "
          f"median: {statistics.median(crbn_overlaps):.1f}, "
          f"min: {min(crbn_overlaps)}, max: {max(crbn_overlaps)}")

    print(f"\nPPIL4 (real interface has {len(REAL_PPIL4)} residues, threshold {len(REAL_PPIL4) / 2:.0f}+ to FOUND THE SPOT)")
    print(f"  FOUND THE SPOT: {ppil4_found}/{len(ppil4_overlaps)} ({ppil4_found / len(ppil4_overlaps):.0%})")
    print(f"  overlap count -- mean: {statistics.mean(ppil4_overlaps):.2f}, "
          f"median: {statistics.median(ppil4_overlaps):.1f}, "
          f"min: {min(ppil4_overlaps)}, max: {max(ppil4_overlaps)}")


if __name__ == "__main__":
    main()
