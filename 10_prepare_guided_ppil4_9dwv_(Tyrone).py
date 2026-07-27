#!/usr/bin/env python3
"""Prepare a separate, 9DWV-guided PPIL4-RRM HADDOCK screen.

This protocol is for exploratory FPFT-2216/9DWV-like molecular-glue hypotheses.
It uses broad ambiguous interface restraints derived from 9DWV residue sets,
not the exact native contact-pair restraints used by the positive control.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

from validate_9dwv_blind_vina import transform_4ci1_to_9dwv


ROOT = Path("haddock3_ternary/guided_ppil4_9DWV")
SOURCE_MANIFEST = Path("haddock3_ternary/ppil4_top250/candidate_manifest.csv")
VINA_POSES = Path("crbn_vina_rf_top500")
VALIDATION = Path("haddock3_ternary/validation_9DWV")
CONTROL = VALIDATION / "bound_control"
REFERENCE = VALIDATION / "blind_redocking/reference_9DWV_ternary.pdb"
DDB1 = CONTROL / "ddb1_9DWV_A.pdb"
CRBN = CONTROL / "crbn_9DWV_B.pdb"
PPIL4 = CONTROL / "ppil4_rrm_9DWV_C.pdb"
EXPERIMENTAL_GLUE = CONTROL / "fpft2216_9DWV_A1BC8_E.pdb"


def protein_residues(chain):
    return [residue for residue in chain if residue.id[0] == " " and "CA" in residue]


def heavy_coordinates(residue) -> np.ndarray:
    return np.asarray([atom.coord for atom in residue.get_atoms() if atom.element.upper() != "H"])


def interface_residue_sets() -> tuple[list[int], list[int]]:
    """Return compact, high-contact 9DWV interface residue sets for AIRs."""
    reference = PDBParser(QUIET=True).get_structure("9dwv", REFERENCE)[0]
    crbn, ppil4 = reference["B"], reference["C"]
    crbn_counts, ppil4_counts = {}, {}
    for residue_c in protein_residues(ppil4):
        atoms_c = heavy_coordinates(residue_c)
        for residue_b in protein_residues(crbn):
            atoms_b = heavy_coordinates(residue_b)
            if len(atoms_c) and len(atoms_b) and np.min(np.linalg.norm(atoms_c[:, None, :] - atoms_b[None, :, :], axis=2)) <= 5.0:
                ppil4_counts[residue_c.id[1]] = ppil4_counts.get(residue_c.id[1], 0) + 1
                crbn_counts[residue_b.id[1]] = crbn_counts.get(residue_b.id[1], 0) + 1
    if not crbn_counts or not ppil4_counts:
        raise RuntimeError("No 9DWV CRBN--PPIL4 interface residues were found.")
    # Keep lines safely below CNS's parser limit.  These are residue *sets*,
    # not pairwise restraints: every selected PPIL4 residue may contact any
    # selected CRBN residue.
    crbn_active = [residue for residue, _ in sorted(crbn_counts.items(), key=lambda item: (-item[1], item[0]))[:8]]
    ppil4_active = [residue for residue, _ in sorted(ppil4_counts.items(), key=lambda item: (-item[1], item[0]))[:5]]
    return sorted(crbn_active), sorted(ppil4_active)


def write_ambiguous_restraints(crbn_residues: list[int], ppil4_residues: list[int]) -> Path:
    """Each PPIL4 interface residue may contact any CRBN interface residue."""
    path = ROOT / "ambiguous_9DWV_interface.tbl"
    crbn_selection = " or ".join(f"resid {residue}" for residue in crbn_residues)
    with path.open("w") as handle:
        handle.write("! Broad ambiguous restraints derived from 9DWV interface residue sets.\n")
        handle.write("! Not exact native contact pairs; use only for 9DWV-like hypotheses.\n")
        for residue in ppil4_residues:
            handle.write(
                f"assign (segid C and resid {residue} and name CA) "
                f"(segid B and ({crbn_selection}) and name CA) 2.0 5.0 0.0\n"
            )
    return path


def pdbqt_atoms(path: Path):
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line[77:].strip().split()[-1] if line[77:].strip() else "C"
        # Vina uses two-character atom types for halogens.  Falling back to
        # atom_type[0] changed Br into B (boron), which made HADDOCK/PRODRG
        # reject every brominated candidate at topology preparation.
        element = {
            "A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H",
            "CL": "Cl", "BR": "Br", "I": "I", "F": "F",
        }.get(atom_type.upper(), atom_type[0].upper())
        coordinate = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        yield element, coordinate


def experimental_centroid() -> np.ndarray:
    points = []
    for line in EXPERIMENTAL_GLUE.read_text().splitlines():
        if line.startswith("HETATM") and line[76:78].strip().upper() != "H":
            points.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.mean(points, axis=0)


def write_transformed_candidate(pose: Path, destination: Path, rotation, translation) -> tuple[float, int]:
    lines, heavy_points = [], []
    for serial, (element, coordinate) in enumerate(pdbqt_atoms(pose), start=1):
        x, y, z = np.dot(coordinate, rotation) + translation
        # Keep both characters of two-letter elements in the PDB atom name.
        # The ligands here contain fewer than 100 atoms.
        atom_name = f"{element.upper()}{serial:02d}"[-4:]
        lines.append(
            f"HETATM{serial:5d} {atom_name:>4} LIG E 502    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2}\n"
        )
        if element != "H":
            heavy_points.append([x, y, z])
    destination.write_text("".join(lines) + "TER\nEND\n")
    return np.asarray(heavy_points), len(lines)


def write_config(destination: Path, ligand: Path, restraint_file: Path, sampling: int) -> Path:
    config = destination / "guided_ppil4_haddock.cfg"
    config.write_text(
        "# Exploratory 9DWV-guided PPIL4-RRM protocol; not a general PPIL4 predictor.\n"
        f"run_dir = \"{destination / 'run'}\"\n"
        "molecules = [\n"
        f"  \"{DDB1}\",\n  \"{CRBN}\",\n  \"{ligand}\",\n  \"{PPIL4}\",\n]\n"
        "ncores = 2\nmode = \"local\"\npostprocess = false\nclean = false\n\n"
        "[topoaa]\nautotoppar = true\nhydrogen_build = \"all\"\n\n"
        "[rigidbody]\n"
        f"ambig_fname = \"{restraint_file}\"\nambig_scale = 20\n"
        "mol_fix_origin_1 = true\nmol_fix_origin_2 = true\nmol_fix_origin_3 = true\n"
        f"cmrest = true\ncmtight = true\nsampling = {sampling}\niniseed = 917\n\n"
        "[emref]\n\n[emscoring]\n"
    )
    return config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-rank", type=int, default=1)
    parser.add_argument("--end-rank", type=int, default=250)
    parser.add_argument("--sampling", type=int, default=200, help="Rigid-body models per eligible candidate.")
    parser.add_argument("--max-pocket-centroid-distance", type=float, default=5.0)
    args = parser.parse_args()
    for required in (SOURCE_MANIFEST, REFERENCE, DDB1, CRBN, PPIL4, EXPERIMENTAL_GLUE):
        if not required.is_file():
            raise FileNotFoundError(required)
    with SOURCE_MANIFEST.open(newline="") as handle:
        candidates = list(csv.DictReader(handle))
    if not 1 <= args.start_rank <= args.end_rank <= len(candidates):
        raise ValueError(f"Ranks must be within 1..{len(candidates)}.")
    ROOT.mkdir(parents=True, exist_ok=True)
    crbn_residues, ppil4_residues = interface_residue_sets()
    restraint_file = write_ambiguous_restraints(crbn_residues, ppil4_residues)
    rotation, translation, fit_rmsd, offset, n_ca = transform_4ci1_to_9dwv()
    reference_centroid = experimental_centroid()
    manifest = []
    for rank in range(args.start_rank, args.end_rank + 1):
        row = candidates[rank - 1]
        molecule, input_rank = row["molecule"], int(row["vina_pose_pdbqt"].split("/")[-1].split("_")[0])
        pose = VINA_POSES / f"{input_rank:02d}_{molecule}_pose.pdbqt"
        run_dir = ROOT / f"rank{rank:03d}_molecule_{molecule}"
        result = {"vina_rank": rank, "molecule": molecule, "canonical_smiles": row["canonical_smiles"], "vina_score_kcal_mol": row["vina_score_kcal_mol"], "pocket_centroid_distance_angstrom": "", "eligible": False, "config": "", "status": ""}
        if not pose.is_file():
            result["status"] = "missing_vina_pose"
            manifest.append(result)
            continue
        run_dir.mkdir(parents=True, exist_ok=True)
        ligand = run_dir / f"rank{rank:03d}_molecule_{molecule}_aligned_to_9DWV.pdb"
        points, atom_count = write_transformed_candidate(pose, ligand, rotation, translation)
        distance = float(np.linalg.norm(np.mean(points, axis=0) - reference_centroid))
        result["pocket_centroid_distance_angstrom"] = f"{distance:.3f}"
        if distance > args.max_pocket_centroid_distance:
            result["status"] = "pocket_incompatible"
            manifest.append(result)
            continue
        config = write_config(run_dir, ligand, restraint_file, args.sampling)
        result.update(eligible=True, config=str(config), status="prepared")
        manifest.append(result)
    output = ROOT / "guided_candidate_manifest.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Prepared {sum(bool(row['eligible']) for row in manifest)} eligible guided configurations in {ROOT}.")
    print(f"Broad interface restraints: {len(crbn_residues)} CRBN residues x {len(ppil4_residues)} PPIL4 residues.")
    print(f"4CI1-to-9DWV CRBN alignment: {fit_rmsd:.3f} A over {n_ca} CA atoms (residue offset {offset}).")
    print(f"Manifest: {output}")


if __name__ == "__main__":
    main()
