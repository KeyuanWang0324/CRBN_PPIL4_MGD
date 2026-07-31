#!/usr/bin/env python3
"""Create one 31-entry table from the fully corrected guided-HADDOCK campaign.

Includes 18 originally successful candidates rerun with corrected ligand
generation, 12 corrected halogen retries, and FPFT-2216 starting from Vina.
The earlier FPFT-2216 *native-pose* control is deliberately excluded.
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path


SHORTLIST = Path("haddock3_ternary/guided_ppil4_9DWV/pre_haddock_triage/guided_haddock_shortlist_30.csv")
CAL_ROOT = Path("haddock3_ternary/guided_ppil4_9DWV_corrected_full_calibration")
RETRY_ROOT = Path("haddock3_ternary/guided_ppil4_9DWV_halogen_retry_baseline")
CAL_MANIFEST = CAL_ROOT / "corrected_calibration_manifest.csv"
RETRY_MANIFEST = RETRY_ROOT / "retry_and_baseline_manifest.csv"
OUTPUT = Path("haddock3_ternary/final_corrected_31_results")


def model_rows(config: Path) -> list[dict[str, str]]:
    table = config.parent / "run/3_emscoring/emscoring.tsv"
    if not table.is_file():
        return []
    with table.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return sorted(rows, key=lambda row: float(row["score"]))


def score_summary(rows: list[dict[str, str]], config: Path) -> dict[str, str]:
    values = [float(row["score"]) for row in rows]
    top = values[: min(10, len(values))]
    return {
        "models_scored": str(len(values)),
        "best_haddock_score": f"{values[0]:.3f}",
        "mean_best_10_haddock_score": f"{statistics.mean(top):.3f}",
        "sd_best_10_haddock_score": f"{statistics.stdev(top):.3f}" if len(top) > 1 else "0.000",
        "median_haddock_score": f"{statistics.median(values):.3f}",
        "worst_haddock_score": f"{values[-1]:.3f}",
        "best_model_structure": rows[0]["structure"],
        "best_model_path": str(config.parent / "run/3_emscoring" / rows[0]["structure"]),
    }


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    for required in (SHORTLIST, CAL_MANIFEST, RETRY_MANIFEST):
        if not required.is_file():
            raise FileNotFoundError(required)
    with SHORTLIST.open(newline="") as handle:
        shortlist = {row["molecule"]: row for row in csv.DictReader(handle)}
    with CAL_MANIFEST.open(newline="") as handle:
        calibration = list(csv.DictReader(handle))
    with RETRY_MANIFEST.open(newline="") as handle:
        retries = [row for row in csv.DictReader(handle) if row["kind"] == "candidate"]

    all_jobs = [("corrected rerun", row) for row in calibration] + [("corrected halogen retry", row) for row in retries]
    complete, incomplete = [], []
    for source, job in all_jobs:
        config = Path(job["guided_config"] if "guided_config" in job else job["config"])
        models = model_rows(config)
        if not models:
            incomplete.append({"source": source, "molecule": job["molecule"], "config": str(config)})
            continue
        is_baseline = job["molecule"] == "FPFT-2216"
        selected = shortlist.get(job["molecule"], {})
        complete.append({
            "source": "FPFT-2216 Vina-to-HADDOCK baseline" if is_baseline else source,
            "molecule": job["molecule"], "label": job["label"],
            "guided_selection_rank": job.get("guided_selection_rank", ""), "vina_rank": job.get("vina_rank", ""),
            "canonical_smiles": job.get("canonical_smiles", ""),
            "rf_active_probability": job.get("rf_active_probability", "not applicable: reference control"),
            "vina_score_kcal_mol": job.get("vina_score_kcal_mol", ""),
            # The calibration manifest stores this directly; the earlier
            # halogen-retry manifest does not, so inherit the unchanged
            # Vina-pose geometry from the original shortlist.
            "pocket_centroid_distance_angstrom": job.get("pocket_centroid_distance_angstrom") or selected.get("pocket_centroid_distance_angstrom", "not applicable: reference control"),
            "ppil4_min_distance_angstrom": selected.get("ppil4_min_distance_angstrom", "not applicable: reference control"),
            "geometry_status": selected.get("geometry_status", "not applicable: reference control"),
            "rf_component": selected.get("rf_component", ""), "vina_component": selected.get("vina_component", ""),
            "pocket_component": selected.get("pocket_component", ""),
            "ppil4_geometry_component": selected.get("ppil4_geometry_component", ""),
            "composite_pre_haddock_score": selected.get("composite_pre_haddock_score", ""),
            "max_similarity_to_other_selected": selected.get("max_similarity_to_other_selected", ""),
            "initial_ligand_pose": "Vina pose transformed from 4CI1 into 9DWV CRBN frame",
            "vina_pose_file": job.get("vina_source", ""), "guided_config": str(config),
            **score_summary(models, config),
        })

    OUTPUT.mkdir(parents=True, exist_ok=True)
    fields = [
        "final_corrected_haddock_rank", "source", "molecule", "label", "guided_selection_rank", "vina_rank",
        "canonical_smiles", "rf_active_probability", "vina_score_kcal_mol", "pocket_centroid_distance_angstrom",
        "ppil4_min_distance_angstrom", "geometry_status", "rf_component", "vina_component", "pocket_component",
        "ppil4_geometry_component", "composite_pre_haddock_score", "max_similarity_to_other_selected",
        "initial_ligand_pose", "models_scored", "best_haddock_score", "mean_best_10_haddock_score",
        "sd_best_10_haddock_score", "median_haddock_score", "worst_haddock_score", "best_model_structure",
        "best_model_path", "vina_pose_file", "guided_config",
    ]
    complete.sort(key=lambda row: float(row["mean_best_10_haddock_score"]))
    for rank, row in enumerate(complete, start=1):
        row["final_corrected_haddock_rank"] = str(rank)
    write_csv(OUTPUT / "final_corrected_31_ranked.csv", complete, fields)
    write_csv(OUTPUT / "incomplete_runs.csv", incomplete, ["source", "molecule", "config"])
    print(f"Completed entries: {len(complete)} of 31 expected.")
    print(f"Ranking: {OUTPUT / 'final_corrected_31_ranked.csv'}")
    if incomplete:
        print(f"Incomplete runs: {OUTPUT / 'incomplete_runs.csv'}")


if __name__ == "__main__":
    main()

