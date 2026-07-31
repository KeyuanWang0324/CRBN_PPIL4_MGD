#!/usr/bin/env python3
"""Run Ryan's 18 finalists plus one supplied molecule through Vina and HADDOCK.

All 19 ligands are newly embedded, converted to PDBQT, docked into the fixed
4CI1 CRBN box, transformed to the 9DWV CRBN frame with corrected halogen
mapping, and run through the same 200-model guided PPIL4 HADDOCK workflow.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

from dock_crbn_vina import DEFAULT_BOX_SIZE, dock_ligand, fixed_box_from_cocrystal, make_3d_sdf, prepare_ligand
from prepare_guided_ppil4_9dwv import experimental_centroid, interface_residue_sets, transform_4ci1_to_9dwv, write_config, write_transformed_candidate
from validate_9dwv_blind_vina import ensure_tools_on_path

RYAN_FINAL = Path("/Users/tyronezong/Downloads/08_final_ternary_with_ligand_results_(Ryan).csv")
RYAN_SCORES = Path("/Users/tyronezong/Downloads/02_crbn_glue_scores_for_03_(Ryan).csv")
ROOT = Path("haddock3_ternary/ryan18_friend19_corrected")
VINA_RECEPTOR = Path("crbn_4ci1_receptor.pdbqt")
VINA_REFERENCE = Path("4CI1.pdb")
FRIEND_SMILES = "O=C1CCC(C(N1)=O)N2CC3=C(C2=O)C=CC=C3C4=CC=C(C=C4)N5CCNCC5"


def write_restraints() -> Path:
    crbn, ppil4 = interface_residue_sets()
    path = ROOT / "ambiguous_9DWV_interface.tbl"
    selection = " or ".join(f"resid {residue}" for residue in crbn)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("! Broad ambiguous CRBN--PPIL4 interface restraints derived from 9DWV.\n" + "".join(
        f"assign (segid C and resid {residue} and name CA) (segid B and ({selection}) and name CA) 2.0 5.0 0.0\n" for residue in ppil4))
    return path


def set_ncores(config: Path, cores: int) -> None:
    text, count = re.subn(r"(?m)^ncores\s*=\s*\d+\s*$", f"ncores = {cores}", config.read_text())
    if count != 1:
        raise RuntimeError(f"Expected one ncores entry in {config}")
    # Housekeeping only: preserve final scores/models while removing completed
    # intermediate modules, so this 19-job campaign fits on the available disk.
    text, clean_count = re.subn(r"(?m)^clean\s*=\s*false\s*$", "clean = true", text)
    if clean_count != 1:
        raise RuntimeError(f"Expected one clean entry in {config}")
    config.write_text(text)


def is_complete(config: Path) -> bool:
    return (config.parent / "run/3_emscoring/emscoring.tsv").is_file()


def input_rows() -> list[dict[str, str]]:
    with RYAN_FINAL.open(newline="") as handle:
        final = list(csv.DictReader(handle))
    with RYAN_SCORES.open(newline="") as handle:
        scores = {row["name"]: row for row in csv.DictReader(handle)}
    if len(final) != 18 or len({row["name"] for row in final}) != 18:
        raise RuntimeError("Ryan final-results file must contain exactly 18 unique candidates.")
    rows = []
    for row in final:
        matched = scores.get(row["name"])
        if not matched or not matched.get("smiles"):
            raise RuntimeError(f"No SMILES match for {row['name']} in Ryan's CRBN-score file.")
        rows.append({
            "kind": "ryan_candidate", "molecule": row["name"], "label": row["name"], "smiles": matched["smiles"],
            "rf_active_probability": matched["p_crbn_glue"], "ryan_cluster_rank": row["cluster_rank"],
            "ryan_cluster_id": row["cluster_id"], "ryan_cluster_size": row["n"],
            "ryan_haddock_score": row["score"], "ryan_dockq": row["dockq"], "ryan_irmsd": row["irmsd"],
            "ryan_fnat": row["fnat"], "ryan_lrmsd": row["lrmsd"],
        })
    rows.append({
        "kind": "positive_control", "molecule": "positive_control", "label": "positive_control", "smiles": FRIEND_SMILES,
        "rf_active_probability": "not supplied", "ryan_cluster_rank": "", "ryan_cluster_id": "", "ryan_cluster_size": "",
        "ryan_haddock_score": "", "ryan_dockq": "", "ryan_irmsd": "", "ryan_fnat": "", "ryan_lrmsd": "",
    })
    return rows


def prepare(sampling: int, cores: int) -> list[dict[str, str]]:
    for required in (RYAN_FINAL, RYAN_SCORES, VINA_RECEPTOR, VINA_REFERENCE):
        if not required.is_file():
            raise FileNotFoundError(required)
    ensure_tools_on_path()
    rotation, translation, fit_rmsd, _, n_ca = transform_4ci1_to_9dwv()
    reference_centroid, restraints = experimental_centroid(), write_restraints()
    centre, box_description = fixed_box_from_cocrystal(VINA_REFERENCE, None)
    vina_dir = ROOT / "vina"
    vina_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for index, row in enumerate(input_rows(), start=1):
        safe = row["molecule"].replace("/", "_")
        sdf, ligand, pose, log = (vina_dir / f"{safe}.sdf", vina_dir / f"{safe}.pdbqt", vina_dir / f"{safe}_pose_4CI1.pdbqt", vina_dir / f"{safe}.vina.log")
        canonical = make_3d_sdf(row["smiles"], sdf)
        prepare_ligand(sdf, ligand)
        vina_score = dock_ligand(ligand, VINA_RECEPTOR, centre, DEFAULT_BOX_SIZE, pose, log)
        folder = ROOT / f"{index:02d}_{safe}"
        folder.mkdir(parents=True, exist_ok=True)
        transformed = folder / f"{safe}_vina_pose_corrected_aligned_to_9DWV.pdb"
        points, _ = write_transformed_candidate(pose, transformed, rotation, translation)
        config = write_config(folder, transformed, restraints, sampling)
        set_ncores(config, cores)
        jobs.append({**row, "canonical_smiles": canonical, "vina_score_kcal_mol": f"{vina_score:.3f}",
                     "pocket_centroid_distance_angstrom": f"{np.linalg.norm(np.mean(points, axis=0) - reference_centroid):.3f}",
                     "vina_box_reference": box_description, "vina_pose_file": str(pose), "guided_config": str(config)})
    print(f"4CI1-to-9DWV CRBN alignment: {fit_rmsd:.3f} A over {n_ca} CA atoms")
    return jobs


def write_manifest(jobs: list[dict[str, str]]) -> Path:
    path = ROOT / "ryan18_friend19_manifest.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(jobs[0]))
        writer.writeheader(); writer.writerows(jobs)
    return path


def execute(jobs: list[dict[str, str]]) -> int:
    executable = Path(os.environ.get("CONDA_PREFIX", "")) / "bin/haddock3"
    if not executable.is_file():
        found = shutil.which("haddock3")
        if not found: raise RuntimeError("Activate my-rdkit-env so haddock3 is available.")
        executable = Path(found)
    failures = 0
    for n, job in enumerate(jobs, 1):
        config = Path(job["guided_config"])
        if is_complete(config):
            print(f"[{n}/19] {job['label']}: already complete; skipped."); continue
        print(f"[{n}/19] starting {job['label']}")
        with (config.parent / "launcher.log").open("a") as log:
            result = subprocess.run([str(executable), str(config), "--restart", "0"], stdout=log, stderr=subprocess.STDOUT, check=False)
        if result.returncode or not is_complete(config):
            failures += 1; print(f"[{n}/19] {job['label']}: FAILED (exit {result.returncode}); continuing.")
        else: print(f"[{n}/19] {job['label']}: complete.")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampling", type=int, default=200); parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.sampling < 1 or args.cores < 1: raise SystemExit("--sampling and --cores must be positive.")
    jobs = prepare(args.sampling, args.cores)
    print(f"Prepared 19/19 matched Vina-to-HADDOCK jobs. Manifest: {write_manifest(jobs)}")
    if args.execute:
        failures = execute(jobs)
        if failures: raise SystemExit(f"{failures} job(s) failed; inspect each launcher.log")


if __name__ == "__main__": main()

