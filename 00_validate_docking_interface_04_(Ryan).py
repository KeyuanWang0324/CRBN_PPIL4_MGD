"""
Population-scale version of 00_validate_docking_interface_(Ryan).py's
agree() check -- that script hardcoded a single example (cand_5, docked
BEFORE the PPIL4 restraint fix in 04/05/06/08). This runs the same
real-interface comparison over every candidate in 04's current Vina
screening output, now that PPIL4 docking is restrained to its real
RRM-domain interface (249-279, from FPFT-2216 / PDB 9DWV) instead of the
old CypA-homology guess.

REAL_CRBN/REAL_PPIL4 are the same ground-truth residue sets from
session4_handson.md (Part 2a) -- pasted verbatim, not recomputed, not
invented.

YOUR_CRBN per candidate comes straight from 04's own saved
crbn_contacts.txt (the CRBN residues that candidate's selected Vina pose
actually touched). YOUR_PPIL4 isn't saved by 04 directly, so it's
recomputed here the exact same way 04 computes it internally
(heavy-atom contacts within 4.5 A, see find_ligand_contacts() in
04_vina_dock_candidates_(Ryan).py) -- from that SAME already-selected
PPIL4 pose, found by matching each model's "REMARK VINA RESULT" affinity
in <candidate>_ppil4_poses.pdbqt against the ppil4_affinity 04 already
recorded in its CSV. No re-docking, nothing invented -- just re-deriving
a value from the real Vina output that 04 didn't happen to persist.

Run with the SYSTEM python:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "00_validate_docking_interface_04_(Ryan).py"
"""
import csv
import os
import statistics

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SCREENING_CSV = os.path.join(SCRIPT_DIR, "04_vina_screening_scores_for_05_(Ryan).csv")
VINA_OUT_DIR = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_novel_candidate")
PPIL4_SOURCE_PDB = os.path.join(SCRIPT_DIR, "PPIL4_alphafold_(Ryan).pdb")

# Real interface residues from the FPFT-2216 / PDB 9DWV structure (see
# session4_handson.md's REAL_CRBN/REAL_PPIL4).
REAL_CRBN = {351, 353, 355, 357, 372, 373, 386, 388, 397, 400}
REAL_PPIL4 = {249, 250, 273, 275, 276, 277, 278, 279}


def agree(your_set, real_set):
    found = your_set & real_set
    verdict = "FOUND THE SPOT" if len(found) >= len(real_set) / 2 else "WRONG SPOT"
    return found, verdict


def read_crbn_contacts(name):
    path = os.path.join(VINA_OUT_DIR, name, "crbn_contacts.txt")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return set(int(x) for x in f.readline().split())


def parse_pdbqt_models(path):
    """Yield (affinity, [(x, y, z), ...]) for each MODEL in a Vina poses pdbqt."""
    models = []
    affinity, atoms = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                affinity, atoms = None, []
            elif line.startswith("REMARK VINA RESULT:"):
                affinity = float(line.split()[3])
            elif line.startswith("ENDMDL"):
                models.append((affinity, atoms))
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return models


def find_contacts(receptor_pdb, atoms, cutoff=4.5):
    """Same logic as 04's find_ligand_contacts()."""
    protein_atoms = []
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                resnum = int(line[22:26])
                protein_atoms.append((resnum, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    contacts = set()
    for resnum, x, y, z in protein_atoms:
        for lx, ly, lz in atoms:
            if (x - lx) ** 2 + (y - ly) ** 2 + (z - lz) ** 2 <= cutoff ** 2:
                contacts.add(resnum)
                break
    return contacts


def read_your_ppil4(name, target_affinity, tolerance=0.01):
    path = os.path.join(VINA_OUT_DIR, name, f"{name}_ppil4_poses.pdbqt")
    if not os.path.exists(path):
        return None
    for affinity, atoms in parse_pdbqt_models(path):
        if affinity is not None and abs(affinity - target_affinity) <= tolerance:
            return find_contacts(PPIL4_SOURCE_PDB, atoms)
    return None


def main():
    with open(SCREENING_CSV, newline="") as f:
        rows = list(csv.DictReader(f))

    crbn_overlaps, ppil4_overlaps = [], []
    crbn_found, ppil4_found = 0, 0
    crbn_missing, ppil4_missing = 0, 0

    for row in rows:
        name = row["name"]

        your_crbn = read_crbn_contacts(name)
        if your_crbn is None:
            crbn_missing += 1
        else:
            found, verdict = agree(your_crbn, REAL_CRBN)
            crbn_overlaps.append(len(found))
            if verdict == "FOUND THE SPOT":
                crbn_found += 1

        your_ppil4 = read_your_ppil4(name, float(row["ppil4_affinity"]))
        if your_ppil4 is None:
            ppil4_missing += 1
        else:
            found, verdict = agree(your_ppil4, REAL_PPIL4)
            ppil4_overlaps.append(len(found))
            if verdict == "FOUND THE SPOT":
                ppil4_found += 1

    n = len(rows)
    print(f"04's screening CSV: {n} candidates\n")

    print(f"CRBN  (real interface has {len(REAL_CRBN)} residues, threshold {len(REAL_CRBN) / 2:.0f}+ to FOUND THE SPOT)")
    print(f"  checked: {len(crbn_overlaps)}  (missing crbn_contacts.txt: {crbn_missing})")
    if crbn_overlaps:
        print(f"  FOUND THE SPOT: {crbn_found}/{len(crbn_overlaps)} "
              f"({crbn_found / len(crbn_overlaps):.0%})")
        print(f"  overlap count -- mean: {statistics.mean(crbn_overlaps):.2f}, "
              f"median: {statistics.median(crbn_overlaps):.1f}, "
              f"min: {min(crbn_overlaps)}, max: {max(crbn_overlaps)}")

    print(f"\nPPIL4 (real interface has {len(REAL_PPIL4)} residues, threshold {len(REAL_PPIL4) / 2:.0f}+ to FOUND THE SPOT)")
    print(f"  checked: {len(ppil4_overlaps)}  (missing ppil4 pose match: {ppil4_missing})")
    if ppil4_overlaps:
        print(f"  FOUND THE SPOT: {ppil4_found}/{len(ppil4_overlaps)} "
              f"({ppil4_found / len(ppil4_overlaps):.0%})")
        print(f"  overlap count -- mean: {statistics.mean(ppil4_overlaps):.2f}, "
              f"median: {statistics.median(ppil4_overlaps):.1f}, "
              f"min: {min(ppil4_overlaps)}, max: {max(ppil4_overlaps)}")


if __name__ == "__main__":
    main()
