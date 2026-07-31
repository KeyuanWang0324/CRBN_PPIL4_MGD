#!/usr/bin/env python3
"""Rerun failed brominated guided-HADDOCK candidates plus an FPFT-2216 control.

The original 30-candidate campaign wrote Vina's ``Br`` atom type as boron.
This script uses the corrected converter in prepare_guided_ppil4_9dwv.py,
creates a separate retry directory, and adds the 9DWV FPFT-2216 ligand as a
protocol baseline.  It never overwrites the original runs.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from pathlib import Path

from prepare_guided_ppil4_9dwv import (
    EXPERIMENTAL_GLUE,
    SOURCE_MANIFEST,
    VINA_POSES,
    interface_residue_sets,
    transform_4ci1_to_9dwv,
    write_config,
    write_transformed_candidate,
)


ORIGINAL_ROOT = Path("haddock3_ternary/guided_ppil4_9DWV")
ROOT = Path("haddock3_ternary/guided_ppil4_9DWV_halogen_retry_baseline")
SHORTLIST = ORIGINAL_ROOT / "pre_haddock_triage/guided_haddock_shortlist_30.csv"
STATUS = ORIGINAL_ROOT / "pre_haddock_triage/guided_haddock_batch_status.csv"


def write_restraints() -> Path:
    crbn, ppil4 = interface_residue_sets()
    path = ROOT / "ambiguous_9DWV_interface.tbl"
    crbn_selection = " or ".join(f"resid {residue}" for residue in crbn)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "! Broad ambiguous CRBN--PPIL4 interface restraints derived from 9DWV.\n"
        + "".join(
            f"assign (segid C and resid {residue} and name CA) "
            f"(segid B and ({crbn_selection}) and name CA) 2.0 5.0 0.0\n"
            for residue in ppil4
        )
    )
    return path


def set_ncores(config: Path, cores: int) -> None:
    text, count = re.subn(r"(?m)^ncores\s*=\s*\d+\s*$", f"ncores = {cores}", config.read_text())
    if count != 1:
        raise RuntimeError(f"Expected one ncores entry in {config}")
    config.write_text(text)


def is_complete(config: Path) -> bool:
    return (config.parent / "run/3_emscoring/emscoring.tsv").is_file()


def make_jobs(sampling: int, cores: int, limit_failed: int | None = None) -> list[dict[str, str]]:
    for required in (SHORTLIST, STATUS, SOURCE_MANIFEST, EXPERIMENTAL_GLUE):
        if not required.is_file():
            raise FileNotFoundError(required)
    with SHORTLIST.open(newline="") as handle:
        shortlist = {row["guided_selection_rank"]: row for row in csv.DictReader(handle)}
    with STATUS.open(newline="") as handle:
        failed = [row for row in csv.DictReader(handle) if row["status"] == "failed"]
    if limit_failed is not None:
        failed = failed[:limit_failed]
    with SOURCE_MANIFEST.open(newline="") as handle:
        source = list(csv.DictReader(handle))
    rotation, translation, _, _, _ = transform_4ci1_to_9dwv()
    restraints = write_restraints()
    jobs = []
    for failed_row in failed:
        selection = failed_row["guided_selection_rank"]
        selected = shortlist[selection]
        vina_rank, molecule = int(selected["vina_rank"]), selected["molecule"]
        original = source[vina_rank - 1]
        input_rank = int(Path(original["vina_pose_pdbqt"]).name.split("_")[0])
        pose = VINA_POSES / f"{input_rank:02d}_{molecule}_pose.pdbqt"
        if not pose.is_file():
            raise FileNotFoundError(pose)
        folder = ROOT / f"rank{vina_rank:03d}_molecule_{molecule}"
        folder.mkdir(parents=True, exist_ok=True)
        ligand = folder / f"rank{vina_rank:03d}_molecule_{molecule}_corrected_halogen.pdb"
        write_transformed_candidate(pose, ligand, rotation, translation)
        config = write_config(folder, ligand, restraints, sampling)
        set_ncores(config, cores)
        jobs.append({
            "kind": "candidate", "guided_selection_rank": selection, "vina_rank": str(vina_rank),
            "molecule": molecule, "label": f"candidate_{molecule}", "config": str(config),
            "canonical_smiles": selected["canonical_smiles"], "rf_active_probability": selected["rf_active_probability"],
            "vina_score_kcal_mol": selected["vina_score_kcal_mol"],
        })

    baseline_folder = ROOT / "baseline_fpft2216_9DWV"
    baseline_folder.mkdir(parents=True, exist_ok=True)
    baseline_ligand = baseline_folder / "fpft2216_9DWV_native_pose.pdb"
    shutil.copyfile(EXPERIMENTAL_GLUE, baseline_ligand)
    baseline_config = write_config(baseline_folder, baseline_ligand, restraints, sampling)
    set_ncores(baseline_config, cores)
    jobs.append({
        "kind": "baseline", "guided_selection_rank": "", "vina_rank": "", "molecule": "FPFT-2216",
        "label": "FPFT-2216_9DWV_native_pose_baseline", "config": str(baseline_config),
        "canonical_smiles": "", "rf_active_probability": "", "vina_score_kcal_mol": "",
    })
    return jobs


def write_manifest(jobs: list[dict[str, str]]) -> Path:
    path = ROOT / "retry_and_baseline_manifest.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(jobs[0]))
        writer.writeheader()
        writer.writerows(jobs)
    return path


def main() -> None:
    global ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampling", type=int, default=200, help="Models for each candidate/control (default: 200).")
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--limit-failed", type=int, help="Prepare only the first N failed candidates (for a short preflight).")
    parser.add_argument("--output-root", type=Path, help="Separate output folder; useful for a non-destructive preflight.")
    parser.add_argument("--execute", action="store_true", help="Run HADDOCK after preparing all jobs.")
    args = parser.parse_args()
    if args.sampling < 1 or args.cores < 1 or (args.limit_failed is not None and args.limit_failed < 1):
        raise SystemExit("--sampling, --cores, and --limit-failed must be positive.")
    if args.output_root:
        ROOT = args.output_root
    jobs = make_jobs(args.sampling, args.cores, args.limit_failed)
    manifest = write_manifest(jobs)
    print(f"Prepared {len(jobs) - 1} corrected failed candidates plus one FPFT-2216 baseline.")
    print(f"Manifest: {manifest}")
    if not args.execute:
        return
    executable = Path(os.environ.get("CONDA_PREFIX", "")) / "bin/haddock3"
    if not executable.is_file():
        found = shutil.which("haddock3")
        if not found:
            raise RuntimeError("Activate my-rdkit-env so haddock3 is available.")
        executable = Path(found)
    failures = 0
    for index, job in enumerate(jobs, start=1):
        config = Path(job["config"])
        if is_complete(config):
            print(f"[{index}/{len(jobs)}] {job['label']}: already complete; skipped.")
            continue
        log = config.parent / "retry_launcher.log"
        print(f"[{index}/{len(jobs)}] starting {job['label']}")
        with log.open("a") as handle:
            result = subprocess.run([str(executable), str(config), "--restart", "0"], stdout=handle, stderr=subprocess.STDOUT, check=False)
        if result.returncode or not is_complete(config):
            failures += 1
            print(f"[{index}/{len(jobs)}] {job['label']}: FAILED (exit {result.returncode}); continuing.")
        else:
            print(f"[{index}/{len(jobs)}] {job['label']}: complete.")
    if failures:
        raise SystemExit(f"{failures} retry/baseline job(s) failed; inspect retry_launcher.log files.")


if __name__ == "__main__":
    main()

