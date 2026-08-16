"""
Dock multiple candidate CRBN-glue candidates into both CRBN's thalidomide
pocket and PPIL4's real RRM-domain interface via Vina, as a fast
pre-filter for the HADDOCK3 ternary-docking step (see
05_haddock3_ternary_novel_candidate_(Ryan).py).

Run with the SYSTEM python (has vina/meeko/rdkit installed), not the
haddock3 venv:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "04_vina_dock_candidates_(Ryan).py"

CRBN receptor: CRBN_receptor_thalidomide_Ryan.pdb (CRBN chain A, apo of
ligand). Box is centered on where thalidomide sits in the reference crystal structure
(CRBN-Thalidomide-SALL4_(Ryan).pdb), same pocket, since each candidate
shares the identical glutarimide-isoindolinone CRBN-binding degron.

PPIL4 receptor: PPIL4_alphafold_(Ryan).pdb, box centered on the same
real RRM-domain interface residues (249-279, from FPFT-2216, PDB 9DWV) used for
restraint generation in 05/06 -- see ppil4_pocket_residues() below.
NOTE: this is a simplification, not a mechanistic model -- in a real
molecular-glue ternary complex the small molecule typically stays bound to
CRBN and presents a new protein-protein interface to the neo-substrate,
rather than independently occupying a pocket on it. Docking each candidate
into PPIL4's (unrelated, enzymatic) active site is a cheap second signal
for ranking, not a claim about the true ternary mechanism.

Since the CRBN and PPIL4 dockings are independent, nothing stops Vina's
best pose for each from using the *same* substructure of the ligand to
make its key contacts -- which is geometrically impossible in reality (the
molecule can't bury the same atoms in two separate protein pockets at
once). To guard against that, both dockings keep their top N_POSES modes
(Vina computes these internally either way), and for each candidate we
search all CRBN-pose x PPIL4-pose combinations for the best-combined-
affinity pair whose ligand-contact atoms don't substantially overlap
(POSE_OVERLAP_THRESHOLD). If no combination clears that bar, the candidate
is flagged rather than silently reporting an inconsistent pair. This only
affects screening/ranking quality -- 05's actual restraints don't depend
on this script's PPIL4 pose at all, only the fixed RRM-domain residue set.

This is a fast screening pass (Vina only, no HADDOCK3) meant to rank many
candidates by combined (CRBN + PPIL4) predicted affinity. 05 then runs the
slow full ternary HADDOCK3 docking on only the top-ranked candidates --
see TOP_N there.
"""
import csv
import math
import os
import subprocess
import sys
import time

SCRIPT_START_TIME = time.time()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SYSTEM_PYTHON = "/Library/Frameworks/Python.framework/Versions/3.12/bin/python3"


def _has_required_packages():
    try:
        import vina, meeko, rdkit  # noqa: F401
        return True
    except ImportError:
        return False


# This script needs vina/meeko/rdkit, which live in the SYSTEM python, not
# the .venv-haddock3 venv used by 05/06/08. If they're missing -- e.g. the
# wrong venv is active, or the IDE's Run button used a different interpreter
# -- relaunch under the system python automatically instead of failing deep
# inside prepare_receptor(). The sys.executable check guards against looping
# if the system python itself is ever missing these packages.
if not _has_required_packages() and sys.executable != SYSTEM_PYTHON:
    os.execv(SYSTEM_PYTHON, [SYSTEM_PYTHON] + sys.argv)

REFERENCE_PDB = os.path.join(SCRIPT_DIR, "CRBN-Thalidomide-SALL4_(Ryan).pdb")
CRBN_RECEPTOR_PDB = os.path.join(SCRIPT_DIR, "CRBN_receptor_thalidomide_Ryan.pdb")
PPIL4_SOURCE_PDB = os.path.join(SCRIPT_DIR, "PPIL4_alphafold_(Ryan).pdb")
OUT_DIR = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_novel_candidate")
SCREENING_SUMMARY_CSV = os.path.join(OUT_DIR, "screening_summary.csv")
# Root-level copy, named after this script per project convention (consumed
# by 05) -- the docking_tmp copy above stays put too, next to the raw
# per-candidate poses/contacts it's derived from.
SCREENING_SUMMARY_CSV_ROOT = os.path.join(SCRIPT_DIR, "04_vina_screening_scores_for_05_(Ryan).csv")

# Vina poses considered per protein when picking a mutually-compatible
# CRBN/PPIL4 pair, and the max fraction of shared ligand-contact atoms
# between them still considered geometrically plausible.
N_POSES = 10
POSE_OVERLAP_THRESHOLD = 0.3

# Same CRBN-glutarimide degron chemotype as thalidomide/lenalidomide/pomalidomide
# (scored P(CRBN-glue)=1.000 by the Step-3 RF classifier), each with a
# different candidate extension -- picked from crbn_glue_compounds_(Ryan).txt.
_FALLBACK_CANDIDATES = [
    ("novel_candidate_1", "O=C1CCC(N2Cc3cc(NC(=O)c4cn5cc(Cl)ccc5n4)ccc3C2=O)C(=O)N1"),
]

# Top-ranked subset from 01's 2500 analogs, filtered by 02+03's
# combined rank-aggregated scoring -- see build_active_candidates() in
# 03_crbn_binder_scaffold_model_(Ryan).py for how this file gets built.
ACTIVE_CANDIDATES_CSV = os.path.join(SCRIPT_DIR, "03_active_candidates_for_04_(Ryan).csv")

# Cap how many of 03_active_candidates_for_04_(Ryan).csv's rows get docked, since each
# one costs two Vina dockings (CRBN + PPIL4). The file is already sorted by
# combined_rank (best first, see build_active_candidates() in 03), so this
# keeps the top TOP_FRACTION of it -- e.g. 0.5 keeps the best 50%. Set to
# None (or 1.0) to run all of them.
#
# Set to None (dock everything): 03's combined_rank is a classifier-based
# pre-filter, not a validated ranking (the binder-classifier half is
# saturated -- 0.82-1.0 across the whole library, see chat discussion --
# so it barely discriminates candidates). Vina docking is cheap (~10s per
# candidate observed), so there's little to gain and real signal to lose
# by pre-cutting with that score before the first real (docking-based)
# discrimination happens here.
TOP_FRACTION = None

# Set to a single (name, smiles) pair to dock ONLY that one candidate --
# e.g. a known positive-control compound -- instead of the active-
# candidates funnel. Bypasses ACTIVE_CANDIDATES_CSV/TOP_FRACTION entirely.
# Its results are NOT written to SCREENING_SUMMARY_CSV/_ROOT (those are the
# funnel's shared ranked output that 05 reads; a one-off manual candidate
# isn't part of that ranking) -- see the MANUAL_CANDIDATE branch in main().
# Leave as None for normal funnel behavior.
MANUAL_CANDIDATE = None

# Skip any molecule already docked (its crbn_contacts.txt exists), so re-running
# 04 only docks what's MISSING instead of re-churning the whole 625-candidate pool.
RESUME = True

# Control molecules docked alongside the candidates (kept OUT of the funnel
# summary CSV). FPFT + isoindolinone already have poses (RESUME skips them); the 3
# negatives get their poses here so 08 can dock them. Same set/SMILES as 08's CONTROLS.
CONTROLS = [
    ("FPFT_2216_positive_control", "COc1cscc1c2cn(nn2)[C@H]3CCC(=O)NC3=O"),
    ("Z6466608628_positive_control", "O=C1CCC(C(N1)=O)N2CC3=C(C2=O)C=CC=C3C4=CC=C(C=C4)N5CCNCC5"),  # positive (isoindolinone-piperazine, Enamine Z6466608628); pose already exists -> RESUME skips
    ("TM-04-064-02_positive_control", "O=C(NCCCCCCCCNC1=CC=CC2=C1C(N(C2=O)C3CCC(NC3=O)=O)=O)CCN4C=NC5=C(N=CN=C54)N"),  # positive (adenine-linked thalidomide); NEW -> 04 docks its pose so 08 can use it
    ("negctrl_no_glutarimide", "CN1Cc2c(cccc2-c2ccc(N3CCNCC3)cc2)C1=O"),
    ("negctrl_no_crbn_ppil4_arm", "COc1cscc1-c1c[nH]nn1"),
    ("negctrl_no_ppil4_arm", "O=C1CCC(N2Cc3ccccc3C2=O)C(=O)N1"),
]

# Cross-check molecules: candidates picked by Tyrone's independent guided
# screen (10_*), run through Ryan's pipeline so both rankings can be compared
# on the SAME molecules under the same locked protocol (PROTOCOL_LOCK.md §8).
# Named tyrone_<his number> so they can never be confused with Ryan's own
# cand_<n> (which are 01's generated analogs, numbered independently).
#
# All three happen to also exist in 01's library under a different name --
# tyrone_118 = cand_1981, tyrone_45 = cand_467, tyrone_6 = cand_85 (identical
# SMILES) -- but they are re-docked here under their own names rather than
# reusing those poses: Vina's search is stochastic, and the existing
# cand_* poses on disk no longer match the affinities recorded for them in
# SCREENING_SUMMARY_CSV_ROOT (that file predates RESUME). Re-docking gives one
# self-consistent pose + affinity + overlap set, which is what 08 then uses.
#
# Kept OUT of the funnel summary (same as CONTROLS -- they aren't part of 03's
# candidate pool) and written to their own CSV instead.
CROSSCHECK = [
    ("tyrone_118", "COc1ccc2c(c1C#N)C(=O)N(C1CCC(=O)NC1=O)C2=O"),
    ("tyrone_45", "Cc1ccc2c(c1Br)C(=O)N(C1CCC(=O)NC1=O)C2=O"),
    ("tyrone_6", "CCc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O"),
]
CROSSCHECK_CSV = os.path.join(SCRIPT_DIR, "04_crosscheck_vina_scores_(Ryan).csv")


def load_candidates():
    if MANUAL_CANDIDATE is not None:
        print(f"MANUAL_CANDIDATE set -- docking only {MANUAL_CANDIDATE[0]}, "
              f"bypassing {ACTIVE_CANDIDATES_CSV}.")
        return [MANUAL_CANDIDATE]
    if os.path.exists(ACTIVE_CANDIDATES_CSV):
        import csv
        with open(ACTIVE_CANDIDATES_CSV) as f:
            candidates = [(row["name"], row["smiles"]) for row in csv.DictReader(f)]
        print(f"Loaded {len(candidates)} candidates from {ACTIVE_CANDIDATES_CSV}")
        if TOP_FRACTION is not None and TOP_FRACTION < 1.0:
            n_keep = max(1, round(len(candidates) * TOP_FRACTION))
            candidates = candidates[:n_keep]
            print(f"Capped to top {n_keep} ({TOP_FRACTION:.0%}) by combined_rank (TOP_FRACTION)")
        return candidates
    print(f"{ACTIVE_CANDIDATES_CSV} not found -- using the single fallback candidate. "
          "Run 01, then 02 and 03, to generate a real screened candidate list.")
    return _FALLBACK_CANDIDATES


CANDIDATES = load_candidates()


def thalidomide_box(reference_pdb, padding=14, min_size=20):
    coords = []
    with open(reference_pdb) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == "EF2" and line[21] == "A":
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    xs, ys, zs = zip(*coords)
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    size_x = max(max(xs) - min(xs) + padding, min_size)
    size_y = max(max(ys) - min(ys) + padding, min_size)
    size_z = max(max(zs) - min(zs) + padding, min_size)
    return (cx, cy, cz), (size_x, size_y, size_z)


def ppil4_pocket_residues():
    """PPIL4's real CRBN-facing interface, taken directly from the actual
    experimental FPFT-2216 ternary complex (PDB 9DWV chain C) instead of a
    homology guess. This replaces the old CypA-active-site-homology mapping
    (which landed in PPIL4's N-terminal cyclophilin-like domain, ~residues
    44-180) -- 00_validate_docking_interface_(Ryan).py showed that mapping
    gets 0/8 overlap with FPFT-2216/9DWV's real contact residues, because
    the real glue-mediated interface is entirely in PPIL4's RRM domain
    (~240-318) instead. These 8 residues (249-279) are FPFT-2216/9DWV's
    real CRBN-contact set (see session4_handson.md's REAL_PPIL4), used
    verbatim -- not widened or re-derived -- since we now have the actual
    structure instead of needing to infer the pocket by homology."""
    return sorted({249, 250, 273, 275, 276, 277, 278, 279})


def residue_box(pdb_path, chain, residues, padding=10, min_size=20):
    residues = set(residues)
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain and int(line[22:26]) in residues:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    xs, ys, zs = zip(*coords)
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    size_x = max(max(xs) - min(xs) + padding, min_size)
    size_y = max(max(ys) - min(ys) + padding, min_size)
    size_z = max(max(zs) - min(zs) + padding, min_size)
    return (cx, cy, cz), (size_x, size_y, size_z)


def prepare_receptor(pdb_path, out_base, center, size):
    # Invoke via `-m` (module path) rather than relying on the
    # mk_prepare_receptor.py console script being on PATH -- meeko installs
    # that script under ~/Library/Python/3.12/bin, which isn't always on
    # PATH depending on how this script gets launched (IDE run button vs.
    # a shell that's sourced ~/.zprofile).
    subprocess.run(
        [sys.executable, "-m", "meeko.cli.mk_prepare_receptor",
         "--read_pdb", pdb_path, "-o", out_base, "-p", "-v",
         "--box_center", str(center[0]), str(center[1]), str(center[2]),
         "--box_size", str(size[0]), str(size[1]), str(size[2])],
        check=True,
    )


def smiles_to_ligand_pdbqt(smiles, out_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid SMILES")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise ValueError("3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]
    with open(out_path, "w") as f:
        f.write(pdbqt_string)


def parse_poses(pdbqt_path):
    """Return each MODEL's ligand atoms as (name, x, y, z, atom_type) tuples,
    in file atom order -- consistent across independent dockings of the same
    ligand pdbqt, so atom index doubles as a stable atom identity."""
    poses = []
    current = None
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current = []
            elif line.startswith("ENDMDL"):
                if current is not None:
                    poses.append(current)
                current = None
            elif (line.startswith("ATOM") or line.startswith("HETATM")) and current is not None:
                name = line[12:16].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atype = line[76:78].strip() if len(line) >= 78 else name[0]
                current.append((name, x, y, z, atype))
    return poses


def dock_ligand(receptor_pdbqt, ligand_pdbqt, center, size, out_poses_pdbqt, exhaustiveness=16, n_poses=N_POSES):
    """Dock and return a list of (affinity, pose_atoms) tuples, best affinity first."""
    from vina import Vina

    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(out_poses_pdbqt, n_poses=n_poses, overwrite=True)
    affinities = [float(e[0]) for e in v.energies(n_poses=n_poses)]
    poses = parse_poses(out_poses_pdbqt)
    return list(zip(affinities, poses))


def pose_atoms_to_pdb_lines(pose_atoms, chain="C", resname="LIG", resnum=900):
    """PDB column spec: atom name cols 13-16, altLoc col 17, resName cols
    18-20, chainID col 22, resSeq cols 23-26. The altLoc column (a blank)
    must be written explicitly between name and resname, or every field
    after it (resName/chain/resSeq) lands one column early -- which a
    strict column-based reader (e.g. PyMOL) will misparse."""
    lines = []
    serial = 9000
    for name, x, y, z, atype in pose_atoms:
        serial += 1
        lines.append(
            f"HETATM{serial:5d} {name:<4s} {resname:>3s} {chain}{resnum:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atype[:2].strip():>2s}\n"
        )
    return lines


def merge_receptor_and_ligand(receptor_pdb, ligand_lines, out_pdb):
    with open(receptor_pdb) as f:
        receptor_lines = [l for l in f if l.startswith("ATOM") or l.startswith("HETATM")]
    with open(out_pdb, "w") as f:
        f.writelines(receptor_lines)
        f.writelines(ligand_lines)
        f.write("END\n")


def find_ligand_contacts(receptor_pdb, pose_atoms, cutoff=4.5):
    """Receptor residue numbers within cutoff of any ligand atom in this pose."""
    protein_atoms = []
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                resnum = int(line[22:26])
                protein_atoms.append((resnum, float(line[30:38]), float(line[38:46]), float(line[46:54])))

    contacts = set()
    for resnum, x, y, z in protein_atoms:
        for _, lx, ly, lz, _ in pose_atoms:
            if (x - lx) ** 2 + (y - ly) ** 2 + (z - lz) ** 2 <= cutoff ** 2:
                contacts.add(resnum)
                break
    return sorted(contacts)


def contact_atom_indices(receptor_pdb, pose_atoms, cutoff=4.5):
    """Ligand atom indices (into pose_atoms) within cutoff of any receptor atom."""
    protein_coords = []
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                protein_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

    contacts = set()
    for i, (_, lx, ly, lz, _) in enumerate(pose_atoms):
        for x, y, z in protein_coords:
            if (x - lx) ** 2 + (y - ly) ** 2 + (z - lz) ** 2 <= cutoff ** 2:
                contacts.add(i)
                break
    return contacts


def best_compatible_pair(crbn_results, ppil4_results, crbn_receptor_pdb, ppil4_receptor_pdb,
                          overlap_threshold=POSE_OVERLAP_THRESHOLD):
    """Search all (CRBN pose, PPIL4 pose) combinations for the best-combined-affinity
    pair whose ligand-contact atoms don't substantially overlap. Neither protein's
    pose is privileged -- both pose sets are searched jointly. Falls back to the
    best combined affinity overall (flagged as inconsistent) if no pair clears
    the threshold, rather than dropping the candidate."""
    crbn_contacts = [contact_atom_indices(crbn_receptor_pdb, pose) for _, pose in crbn_results]
    ppil4_contacts = [contact_atom_indices(ppil4_receptor_pdb, pose) for _, pose in ppil4_results]

    pairs = []
    for i, (crbn_aff, _) in enumerate(crbn_results):
        for j, (ppil4_aff, _) in enumerate(ppil4_results):
            shared = crbn_contacts[i] & ppil4_contacts[j]
            union = crbn_contacts[i] | ppil4_contacts[j]
            overlap = len(shared) / len(union) if union else 0.0
            pairs.append((overlap <= overlap_threshold, crbn_aff + ppil4_aff, i, j, crbn_aff, ppil4_aff, overlap, shared))

    consistent_pairs = [p for p in pairs if p[0]]
    pool = consistent_pairs if consistent_pairs else pairs
    _, combined, i, j, crbn_aff, ppil4_aff, overlap, shared = min(pool, key=lambda p: p[1])
    shared_atom_names = [crbn_results[i][1][idx][0] for idx in sorted(shared)]
    return {
        "crbn_pose": i, "ppil4_pose": j, "crbn_affinity": crbn_aff, "ppil4_affinity": ppil4_aff,
        "combined_affinity": combined, "overlap": overlap, "consistent": bool(consistent_pairs),
        "overlap_atoms": shared_atom_names,
    }


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


def check_environment():
    """Fail fast with a clear message if run under the wrong interpreter.

    This script needs the SYSTEM python (vina/meeko/rdkit installed there),
    not the .venv-haddock3 venv used by 05/06/08 -- that venv is for the
    haddock3 CLI only and doesn't have these packages.
    """
    missing = []
    for module in ("vina", "meeko", "rdkit"):
        try:
            __import__(module)
        except ImportError:
            missing.append(module)
    if missing:
        sys.exit(
            f"Missing packages: {', '.join(missing)}. You're running this with:\n"
            f"  {sys.executable}\n"
            "This script needs the SYSTEM python, not the .venv-haddock3 venv. Run:\n"
            '  /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 '
            '"04_vina_dock_candidates_(Ryan).py"'
        )


def main():
    check_environment()
    os.makedirs(OUT_DIR, exist_ok=True)

    print("== Computing CRBN docking box from reference thalidomide position ==")
    crbn_center, crbn_size = thalidomide_box(REFERENCE_PDB)
    print(f"Box center: {crbn_center}, size: {crbn_size}")

    crbn_receptor_base = os.path.join(OUT_DIR, "CRBN_receptor")
    print("== Preparing CRBN Vina receptor ==")
    prepare_receptor(CRBN_RECEPTOR_PDB, crbn_receptor_base, crbn_center, crbn_size)
    crbn_receptor_pdbqt = crbn_receptor_base + ".pdbqt"

    print("\n== Computing PPIL4 docking box from real RRM-domain interface residues (FPFT-2216) ==")
    ppil4_active = ppil4_pocket_residues()
    print("PPIL4 pocket residues:", ppil4_active)
    ppil4_center, ppil4_size = residue_box(PPIL4_SOURCE_PDB, "A", ppil4_active)
    print(f"Box center: {ppil4_center}, size: {ppil4_size}")

    ppil4_receptor_base = os.path.join(OUT_DIR, "PPIL4_receptor")
    print("== Preparing PPIL4 Vina receptor ==")
    prepare_receptor(PPIL4_SOURCE_PDB, ppil4_receptor_base, ppil4_center, ppil4_size)
    ppil4_receptor_pdbqt = ppil4_receptor_base + ".pdbqt"

    # Run list: candidates + controls (controls excluded from the funnel summary).
    # RESUME skips anything already docked, so this only docks what's missing.
    control_names, crosscheck_names = set(), set()
    run_list = list(CANDIDATES)
    if MANUAL_CANDIDATE is None:
        control_names = {n for n, _ in CONTROLS}
        crosscheck_names = {n for n, _ in CROSSCHECK}
        have = {n for n, _ in run_list}
        run_list += [(n, s) for n, s in CONTROLS + CROSSCHECK if n not in have]

    prior_summary = {}
    if MANUAL_CANDIDATE is None and os.path.exists(SCREENING_SUMMARY_CSV_ROOT):
        with open(SCREENING_SUMMARY_CSV_ROOT) as f:
            prior_summary = {row["name"]: row for row in csv.DictReader(f)}

    results = []
    skipped = 0
    loop_start = time.time()
    for i, (candidate_name, candidate_smiles) in enumerate(run_list, 1):
        candidate_dir = os.path.join(OUT_DIR, candidate_name)
        contacts_path = os.path.join(candidate_dir, "crbn_contacts.txt")
        if RESUME and os.path.exists(contacts_path):
            skipped += 1
            if skipped <= 5 or skipped % 100 == 0:
                print(f"[{i}/{len(run_list)}] {candidate_name} already docked -- skipping "
                      f"({skipped} skipped so far).")
            continue
        done = i - 1 - skipped
        elapsed = time.time() - loop_start
        eta = (elapsed / done) * (len(run_list) - i + 1) if done > 0 else 0
        tag = ("control" if candidate_name in control_names else
               "crosscheck" if candidate_name in crosscheck_names else "candidate")
        print(f"\n== [{i}/{len(run_list)}] {tag}: {candidate_name} "
              f"({elapsed:.0f}s in this step, ~{eta:.0f}s remaining | "
              f"{time.time() - SCRIPT_START_TIME:.0f}s total script time) ==")
        os.makedirs(candidate_dir, exist_ok=True)

        print(f"== Preparing ligand: {candidate_name} ==")
        ligand_pdbqt = os.path.join(candidate_dir, f"{candidate_name}.pdbqt")
        smiles_to_ligand_pdbqt(candidate_smiles, ligand_pdbqt)

        print(f"== Docking against CRBN with Vina ({N_POSES} poses) ==")
        crbn_poses_pdbqt = os.path.join(candidate_dir, f"{candidate_name}_crbn_poses.pdbqt")
        crbn_results = dock_ligand(crbn_receptor_pdbqt, ligand_pdbqt, crbn_center, crbn_size, crbn_poses_pdbqt)

        print(f"== Docking against PPIL4 with Vina ({N_POSES} poses) ==")
        ppil4_poses_pdbqt = os.path.join(candidate_dir, f"{candidate_name}_ppil4_poses.pdbqt")
        ppil4_results = dock_ligand(ppil4_receptor_pdbqt, ligand_pdbqt, ppil4_center, ppil4_size, ppil4_poses_pdbqt)

        print("== Selecting best geometrically-compatible CRBN/PPIL4 pose pair ==")
        selection = best_compatible_pair(crbn_results, ppil4_results, CRBN_RECEPTOR_PDB, PPIL4_SOURCE_PDB)
        status = "consistent" if selection["consistent"] else "CONFLICT -- no compatible pair found"
        print(f"CRBN pose #{selection['crbn_pose']} ({selection['crbn_affinity']:.2f} kcal/mol), "
              f"PPIL4 pose #{selection['ppil4_pose']} ({selection['ppil4_affinity']:.2f} kcal/mol), "
              f"ligand-atom overlap {selection['overlap']:.0%} [{status}]")

        crbn_pose_atoms = crbn_results[selection["crbn_pose"]][1]

        print("== Writing selected CRBN pose and merging with CRBN receptor ==")
        ligand_lines = pose_atoms_to_pdb_lines(crbn_pose_atoms)
        merged_pdb = os.path.join(candidate_dir, "CRBN_candidate_complex.pdb")
        merge_receptor_and_ligand(CRBN_RECEPTOR_PDB, ligand_lines, merged_pdb)
        print(f"Wrote {merged_pdb}")

        print("== Finding CRBN contact residues ==")
        contacts = find_ligand_contacts(CRBN_RECEPTOR_PDB, crbn_pose_atoms)
        print("CRBN active (candidate-contact) residues:", contacts)

        contacts_path = os.path.join(candidate_dir, "crbn_contacts.txt")
        with open(contacts_path, "w") as f:
            f.write(" ".join(str(r) for r in contacts) + "\n")
            f.write(f"{selection['crbn_affinity']}\n")
        print(f"Wrote {contacts_path}")

        results.append({
            "name": candidate_name, "crbn_affinity": selection["crbn_affinity"],
            "ppil4_affinity": selection["ppil4_affinity"], "combined_affinity": selection["combined_affinity"],
            "overlap": selection["overlap"], "consistent": selection["consistent"],
            "overlap_atoms": selection["overlap_atoms"],
            "n_contacts": len(contacts),
        })

    print(f"\n== Docked {len(results)} molecule(s) this run; skipped {skipped} already-docked ==")
    if results:
        results.sort(key=lambda r: r["combined_affinity"])
        print_table(
            results,
            [("name", "name"),
             ("crbn (kcal/mol)", lambda r: f"{r['crbn_affinity']:.2f}"),
             ("ppil4 (kcal/mol)", lambda r: f"{r['ppil4_affinity']:.2f}"),
             ("combined", lambda r: f"{r['combined_affinity']:.2f}"),
             ("overlap", lambda r: f"{r['overlap']:.0%}"),
             ("consistent", lambda r: "yes" if r["consistent"] else "NO"),
             ("n_contacts", lambda r: str(r["n_contacts"])),
             ("overlap_atoms", lambda r: ",".join(r["overlap_atoms"]) if r["overlap_atoms"] else "-")],
            title="Vina results docked this run (best combined affinity first)",
        )
        flagged = [r["name"] for r in results if not r["consistent"]]
        if flagged:
            print(f"\nWARNING: no geometrically-compatible CRBN/PPIL4 pose pair found for: {', '.join(flagged)} "
                  "-- reported affinities use the best available (overlapping) pair.")

    if MANUAL_CANDIDATE is not None:
        print(f"\nMANUAL_CANDIDATE set -- not writing the funnel summary (one-offs aren't part of that "
              f"ranking). Docking output under {os.path.join(OUT_DIR, MANUAL_CANDIDATE[0])}/.")
    else:
        # Merge newly-docked candidate rows into the preserved prior summary, dropping
        # controls -- so a RESUME run never clobbers the candidates it skipped.
        fields = ["name", "crbn_affinity", "ppil4_affinity", "combined_affinity",
                  "overlap", "consistent", "n_contacts", "overlap_atoms"]

        def as_row(r):
            return {"name": r["name"], "crbn_affinity": r["crbn_affinity"],
                    "ppil4_affinity": r["ppil4_affinity"], "combined_affinity": r["combined_affinity"],
                    "overlap": r["overlap"], "consistent": r["consistent"],
                    "n_contacts": r["n_contacts"], "overlap_atoms": ",".join(r["overlap_atoms"])}

        def write_summary(out_paths, rows_by_name, label):
            rows_out = sorted(rows_by_name.values(), key=lambda x: float(x["combined_affinity"]))
            for out_path in out_paths:
                with open(out_path, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=fields)
                    writer.writeheader()
                    for row in rows_out:
                        writer.writerow(row)
                print(f"Wrote {out_path} ({len(rows_out)} {label})")

        # Merge newly-docked candidate rows into the preserved prior summary, dropping
        # controls and cross-check molecules -- so a RESUME run never clobbers the
        # candidates it skipped, and neither set leaks into the funnel's ranked pool.
        excluded = control_names | crosscheck_names
        merged = {n: row for n, row in prior_summary.items() if n not in excluded}
        for r in results:
            if r["name"] not in excluded:
                merged[r["name"]] = as_row(r)
        write_summary((SCREENING_SUMMARY_CSV, SCREENING_SUMMARY_CSV_ROOT), merged, "candidates")

        # Cross-check molecules get their own CSV (same columns). Merged the same
        # way, so a RESUME run that skips an already-docked one keeps its numbers.
        if crosscheck_names:
            prior_crosscheck = {}
            if os.path.exists(CROSSCHECK_CSV):
                with open(CROSSCHECK_CSV) as f:
                    prior_crosscheck = {row["name"]: row for row in csv.DictReader(f)}
            cross = {n: row for n, row in prior_crosscheck.items() if n in crosscheck_names}
            for r in results:
                if r["name"] in crosscheck_names:
                    cross[r["name"]] = as_row(r)
            if cross:
                write_summary((CROSSCHECK_CSV,), cross, "cross-check molecules")

    total = time.time() - SCRIPT_START_TIME
    print(f"\nTotal script runtime: {total:.0f}s ({total / 60:.1f} min)")


if __name__ == "__main__":
    main()
