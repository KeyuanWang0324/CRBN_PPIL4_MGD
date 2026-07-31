#!/usr/bin/env python3
"""Calibrate all original successful candidates with corrected ligand generation.

The 18 originally completed candidates are rebuilt from their existing Vina
poses using the corrected PDBQT element/atom-name conversion. FPFT-2216 is
run anew from SMILES through 3D preparation, ligand PDBQT generation, fixed-box
Vina docking, 4CI1-to-9DWV transformation, and the same guided HADDOCK run.
Nothing in prior campaign directories is overwritten.
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
from prepare_guided_ppil4_9dwv import (
    EXPERIMENTAL_GLUE,
    SOURCE_MANIFEST,
    VINA_POSES,
    experimental_centroid,
    interface_residue_sets,
    transform_4ci1_to_9dwv,
    write_config,
    write_transformed_candidate,
)
from validate_9dwv_blind_vina import SMILES as FPFT2216_SMILES, ensure_tools_on_path


ORIGINAL_ROOT = Path("haddock3_ternary/guided_ppil4_9DWV")
ROOT = Path("haddock3_ternary/guided_ppil4_9DWV_corrected_full_calibration")
SHORTLIST = ORIGINAL_ROOT / "pre_haddock_triage/guided_haddock_shortlist_30.csv"
STATUS = ORIGINAL_ROOT / "pre_haddock_triage/guided_haddock_batch_status.csv"
VINA_RECEPTOR = Path("crbn_4ci1_receptor.pdbqt")
VINA_REFERENCE = Path("4CI1.pdb")


def write_restraints() -> Path:
    crbn, ppil4 = interface_residue_sets()
    path = ROOT / "ambiguous_9DWV_interface.tbl"
    selection = " or ".join(f"resid {residue}" for residue in crbn)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "! Broad ambiguous CRBN--PPIL4 interface restraints derived from 9DWV.\n"
        + "".join(
            f"assign (segid C and resid {residue} and name CA) "
            f"(segid B and ({selection}) and name CA) 2.0 5.0 0.0\n"
            for residue in ppil4
        )
    )
    return path


def set_ncores(config: Path, cores: int) -> None:
    text, count = re.subn(r"(?m)^ncores\s*=\s*\d+\s*$", f"ncores = {cores}", config.read_text())
    if count != 1:
        raise RuntimeError(f"Expected exactly one ncores entry in {config}")
    config.write_text(text)


def is_complete(config: Path) -> bool:
    return (config.parent / "run/3_emscoring/emscoring.tsv").is_file()


def fpft_vina_pose() -> tuple[Path, str, float]:
    """Run the FPFT-2216 control from Vina, not from deposited coordinates."""
    vina_dir = ROOT / "fpft2216_from_vina"
    vina_dir.mkdir(parents=True, exist_ok=True)
    sdf = vina_dir / "fpft2216.sdf"
    ligand_pdbqt = vina_dir / "fpft2216.pdbqt"
    pose = vina_dir / "fpft2216_vina_pose_4CI1.pdbqt"
    log = vina_dir / "fpft2216_vina.log"
    canonical = make_3d_sdf(FPFT2216_SMILES, sdf)
    prepare_ligand(sdf, ligand_pdbqt)
    centre, _ = fixed_box_from_cocrystal(VINA_REFERENCE, None)
    score = dock_ligand(ligand_pdbqt, VINA_RECEPTOR, centre, DEFAULT_BOX_SIZE, pose, log)
    return pose, canonical, score


def prepare_jobs(sampling: int, cores: int) -> list[dict[str, str]]:
    for path in (SHORTLIST, STATUS, SOURCE_MANIFEST, VINA_RECEPTOR, VINA_REFERENCE, EXPERIMENTAL_GLUE):
        if not path.is_file():
            raise FileNotFoundError(path)
    ensure_tools_on_path()
    with SHORTLIST.open(newline="") as handle:
        shortlist = {row["guided_selection_rank"]: row for row in csv.DictReader(handle)}
    with STATUS.open(newline="") as handle:
        original_completed = [row for row in csv.DictReader(handle) if row["status"] == "complete"]
    with SOURCE_MANIFEST.open(newline="") as handle:
        source = list(csv.DictReader(handle))
    rotation, translation, fit_rmsd, _, n_ca = transform_4ci1_to_9dwv()
    reference_centroid = experimental_centroid()
    restraints = write_restraints()
    jobs = []

    for status in original_completed:
        selected = shortlist[status["guided_selection_rank"]]
        vina_rank, molecule = int(selected["vina_rank"]), selected["molecule"]
        input_rank = int(Path(source[vina_rank - 1]["vina_pose_pdbqt"]).name.split("_")[0])
        vina_pose = VINA_POSES / f"{input_rank:02d}_{molecule}_pose.pdbqt"
        if not vina_pose.is_file():
            raise FileNotFoundError(vina_pose)
        folder = ROOT / f"rank{vina_rank:03d}_molecule_{molecule}"
        folder.mkdir(parents=True, exist_ok=True)
        ligand = folder / f"rank{vina_rank:03d}_molecule_{molecule}_corrected.pdb"
        points, _ = write_transformed_candidate(vina_pose, ligand, rotation, translation)
        config = write_config(folder, ligand, restraints, sampling)
        set_ncores(config, cores)
        jobs.append({
            "kind": "candidate", "molecule": molecule, "label": f"candidate_{molecule}",
            "guided_selection_rank": status["guided_selection_rank"], "vina_rank": str(vina_rank),
            "canonical_smiles": selected["canonical_smiles"], "rf_active_probability": selected["rf_active_probability"],
            "vina_score_kcal_mol": selected["vina_score_kcal_mol"],
            "pocket_centroid_distance_angstrom": f"{np.linalg.norm(np.mean(points, axis=0) - reference_centroid):.3f}",
            "vina_source": str(vina_pose), "guided_config": str(config),
        })

    vina_pose, canonical, vina_score = fpft_vina_pose()
    folder = ROOT / "baseline_fpft2216_from_vina"
    folder.mkdir(parents=True, exist_ok=True)
    ligand = folder / "fpft2216_vina_pose_corrected_aligned_to_9DWV.pdb"
    points, _ = write_transformed_candidate(vina_pose, ligand, rotation, translation)
    config = write_config(folder, ligand, restraints, sampling)
    set_ncores(config, cores)
    jobs.append({
        "kind": "vina_baseline", "molecule": "FPFT-2216", "label": "FPFT-2216_Vina_to_HADDOCK_baseline",
        "guided_selection_rank": "", "vina_rank": "", "canonical_smiles": canonical,
        "rf_active_probability": "not applicable: reference control", "vina_score_kcal_mol": f"{vina_score:.3f}",
        "pocket_centroid_distance_angstrom": f"{np.linalg.norm(np.mean(points, axis=0) - reference_centroid):.3f}",
        "vina_source": str(vina_pose), "guided_config": str(config),
    })
    print(f"4CI1-to-9DWV CRBN alignment: {fit_rmsd:.3f} A over {n_ca} CA atoms")
    return jobs


def write_manifest(jobs: list[dict[str, str]]) -> Path:
    manifest = ROOT / "corrected_calibration_manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(jobs[0]))
        writer.writeheader()
        writer.writerows(jobs)
    return manifest


def execute(jobs: list[dict[str, str]]) -> int:
    executable = Path(os.environ.get("CONDA_PREFIX", "")) / "bin/haddock3"
    if not executable.is_file():
        candidate = shutil.which("haddock3")
        if not candidate:
            raise RuntimeError("Activate my-rdkit-env so haddock3 is on PATH.")
        executable = Path(candidate)
    failures = 0
    for index, job in enumerate(jobs, start=1):
        config = Path(job["guided_config"])
        if is_complete(config):
            print(f"[{index}/{len(jobs)}] {job['label']}: already complete; skipped.")
            continue
        print(f"[{index}/{len(jobs)}] starting {job['label']}")
        with (config.parent / "calibration_launcher.log").open("a") as log:
            result = subprocess.run([str(executable), str(config), "--restart", "0"], stdout=log, stderr=subprocess.STDOUT, check=False)
        if result.returncode or not is_complete(config):
            failures += 1
            print(f"[{index}/{len(jobs)}] {job['label']}: FAILED (exit {result.returncode}); continuing.")
        else:
            print(f"[{index}/{len(jobs)}] {job['label']}: complete.")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampling", type=int, default=200)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.sampling < 1 or args.cores < 1:
        raise SystemExit("--sampling and --cores must be positive.")
    jobs = prepare_jobs(args.sampling, args.cores)
    manifest = write_manifest(jobs)
    print(f"Prepared {len(jobs) - 1} original candidates plus FPFT-2216 from Vina.")
    print(f"Manifest: {manifest}")
    if args.execute:
        failures = execute(jobs)
        if failures:
            raise SystemExit(f"{failures} run(s) failed; inspect each calibration_launcher.log")


if __name__ == "__main__":
    main()

