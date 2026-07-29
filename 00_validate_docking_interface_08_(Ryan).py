"""
00_validate_docking_interface, run against 08's actual FINAL ranked
output (08_final_ternary_with_ligand_results_(Ryan).csv) -- the top 5
candidates by 08's own HADDOCK score, the pipeline's real final ranking
(see chat discussion for why 06's dockq was replaced by 08's
ligand-inclusive score). Unlike 00_validate_docking_interface_04/05
(which report aggregate stats over the whole screening pool), this
reports each of the top 5 individually, matching the original single-
candidate 00_validate_docking_interface_(Ryan).py's per-candidate agree()
output -- 5 candidates is small enough to just look at directly.

REAL_CRBN/REAL_PPIL4 are the same ground-truth residue sets from
session4_handson.md (Part 2a) -- pasted verbatim, not recomputed, not
invented.

YOUR_CRBN/YOUR_PPIL4 per candidate are heavy-atom (no hydrogens) contacts
within 5 A across chain A (CRBN) and chain B (PPIL4), from that
candidate's own rank-1 model out of 08's run1/9_caprieval/capri_ss.tsv.
08's models have a 3rd chain (C = the ligand); it's ignored here since
REAL_CRBN/REAL_PPIL4 are about the CRBN-PPIL4 protein-protein interface
specifically, same as 00_validate_docking_interface_05's methodology.

Run with the SYSTEM python:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "00_validate_docking_interface_08_(Ryan).py"
"""
import csv
import gzip
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

RESULTS_CSV = os.path.join(SCRIPT_DIR, "08_final_ternary_with_ligand_results_(Ryan).csv")
RUN_DIR_BASE = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run")

TOP_N = 5

# Real interface residues from the FPFT-2216 / PDB 9DWV structure (see
# session4_handson.md's REAL_CRBN/REAL_PPIL4).
REAL_CRBN = {351, 353, 355, 357, 372, 373, 386, 388, 397, 400}
REAL_PPIL4 = {249, 250, 273, 275, 276, 277, 278, 279}

CUTOFF = 5.0  # heavy-atom contact distance, matches how REAL_CRBN/REAL_PPIL4 were measured


def agree(your_set, real_set, label):
    found = your_set & real_set
    print(f"{label}: found {len(found)} of {len(real_set)} real contact residues")
    print(f"  matched: {sorted(found)}")
    if len(found) >= len(real_set) / 2:
        print("FOUND THE SPOT")
    else:
        print("WRONG SPOT")


def find_top_model_path(candidate_name):
    """The rank-1 model from 08's ligand-inclusive run for this candidate
    (same TSV lookup logic as 07/09 -- 08's STEP_PLAN has its final
    caprieval at step index 9)."""
    caprieval_dir = os.path.join(RUN_DIR_BASE, candidate_name, "run1", "9_caprieval")
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
    """(resnum, x, y, z) for heavy (non-hydrogen) ATOM records on one
    chain -- HETATM (e.g. the structural Zn) excluded too, same as
    00_validate_docking_interface_05's methodology."""
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
    top_rows = rows[:TOP_N]
    print(f"Checking 08's top {len(top_rows)} (of {len(rows)}) finalists by HADDOCK score "
          f"against the real FPFT-2216/9DWV interface:\n")

    for rank, row in enumerate(top_rows, 1):
        name = row["name"]
        print(f"=== #{rank}: {name} (score={row['score']}) ===")
        model_path = find_top_model_path(name)
        if model_path is None:
            print("  no rank-1 model found -- skipping\n")
            continue

        lines = read_pdb_lines(model_path)
        crbn_atoms = chain_heavy_atoms(lines, "A")
        ppil4_atoms = chain_heavy_atoms(lines, "B")

        your_crbn = chain_contacts(crbn_atoms, ppil4_atoms)
        your_ppil4 = chain_contacts(ppil4_atoms, crbn_atoms)

        agree(your_crbn, REAL_CRBN, "  CRBN ")
        agree(your_ppil4, REAL_PPIL4, "  PPIL4")
        print()


if __name__ == "__main__":
    main()
