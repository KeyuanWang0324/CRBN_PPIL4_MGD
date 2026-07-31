#!/usr/bin/env python3
"""Merge the original screen, corrected halogen retries, and FPFT-2216 control.

Rank is the mean final HADDOCK score of the best ten models. Lower is better
only within this identical, exploratory guided-HADDOCK protocol. FPFT-2216 is
labelled as a *native-pose protocol baseline*: it starts from its experimental
9DWV coordinates, whereas candidates start from transformed Vina poses.
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path


ORIGINAL = Path("haddock3_ternary/guided_ppil4_9DWV/pre_haddock_triage/haddock_results/completed_candidates_ranked.csv")
SHORTLIST = Path("haddock3_ternary/guided_ppil4_9DWV/pre_haddock_triage/guided_haddock_shortlist_30.csv")
ORIGINAL_STATUS = Path("haddock3_ternary/guided_ppil4_9DWV/pre_haddock_triage/guided_haddock_batch_status.csv")
RETRY_ROOT = Path("haddock3_ternary/guided_ppil4_9DWV_halogen_retry_baseline")
MANIFEST = RETRY_ROOT / "retry_and_baseline_manifest.csv"
OUTPUT = RETRY_ROOT / "haddock_results"


def score_rows(config: Path) -> list[dict[str, str]]:
    path = config.parent / "run/3_emscoring/emscoring.tsv"
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return sorted(rows, key=lambda row: float(row["score"]))


def stats(rows: list[dict[str, str]], config: Path) -> dict[str, str]:
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


def write(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    for required in (ORIGINAL, SHORTLIST, ORIGINAL_STATUS, MANIFEST):
        if not required.is_file():
            raise FileNotFoundError(required)
    with ORIGINAL.open(newline="") as handle:
        original = list(csv.DictReader(handle))
    with MANIFEST.open(newline="") as handle:
        jobs = list(csv.DictReader(handle))
    with SHORTLIST.open(newline="") as handle:
        shortlist = {row["guided_selection_rank"]: row for row in csv.DictReader(handle)}
    with ORIGINAL_STATUS.open(newline="") as handle:
        original_status = {row["guided_selection_rank"]: row for row in csv.DictReader(handle)}

    retry_rows, status_rows = [], []
    for job in jobs:
        config = Path(job["config"])
        models = score_rows(config)
        state = "complete" if models else "not_complete"
        status_rows.append({**job, "retry_status": state, "models_scored": str(len(models))})
        if not models:
            continue
        source_row = shortlist.get(job["guided_selection_rank"], {})
        first_attempt = original_status.get(job["guided_selection_rank"], {})
        row = {
            "source": "FPFT-2216 native-pose protocol baseline" if job["kind"] == "baseline" else "corrected halogen retry",
            "molecule": job["molecule"], "label": job["label"],
            "guided_selection_rank": job["guided_selection_rank"], "vina_rank": job["vina_rank"],
            "canonical_smiles": job["canonical_smiles"],
            "rf_active_probability": job["rf_active_probability"], "vina_score_kcal_mol": job["vina_score_kcal_mol"],
            "pocket_centroid_distance_angstrom": source_row.get("pocket_centroid_distance_angstrom", "not applicable: experimental native pose"),
            "ppil4_min_distance_angstrom": source_row.get("ppil4_min_distance_angstrom", "not applicable: experimental native pose"),
            "geometry_status": source_row.get("geometry_status", "not applicable: experimental native pose"),
            "rf_component": source_row.get("rf_component", ""), "vina_component": source_row.get("vina_component", ""),
            "pocket_component": source_row.get("pocket_component", ""), "ppil4_geometry_component": source_row.get("ppil4_geometry_component", ""),
            "composite_pre_haddock_score": source_row.get("composite_pre_haddock_score", ""),
            "max_similarity_to_other_selected": source_row.get("max_similarity_to_other_selected", ""),
            "initial_ligand_pose": "9DWV experimental FPFT-2216 pose" if job["kind"] == "baseline" else "Vina pose transformed from 4CI1 into 9DWV CRBN frame",
            "original_batch_status": first_attempt.get("status", "not applicable"),
            "original_started_at": first_attempt.get("started_at", ""), "original_finished_at": first_attempt.get("finished_at", ""),
            "retry_status": state, "guided_config": job["config"],
            **stats(models, config),
        }
        retry_rows.append(row)

    combined = []
    for row in original:
        config = Path(row["guided_config"])
        models = score_rows(config)
        if not models:
            raise RuntimeError(f"Missing final score table for original completed candidate {row['molecule']}: {config}")
        selected = shortlist[row["guided_selection_rank"]]
        combined.append({
            "source": "original candidate run", "molecule": row["molecule"], "label": f"candidate_{row['molecule']}",
            "guided_selection_rank": row["guided_selection_rank"], "vina_rank": row["vina_rank"],
            "canonical_smiles": row["canonical_smiles"], "rf_active_probability": row["rf_active_probability"],
            "vina_score_kcal_mol": row["vina_score_kcal_mol"],
            "pocket_centroid_distance_angstrom": selected["pocket_centroid_distance_angstrom"],
            "ppil4_min_distance_angstrom": selected["ppil4_min_distance_angstrom"], "geometry_status": selected["geometry_status"],
            "rf_component": selected["rf_component"], "vina_component": selected["vina_component"],
            "pocket_component": selected["pocket_component"], "ppil4_geometry_component": selected["ppil4_geometry_component"],
            "composite_pre_haddock_score": selected["composite_pre_haddock_score"],
            "max_similarity_to_other_selected": selected["max_similarity_to_other_selected"],
            "initial_ligand_pose": "Vina pose transformed from 4CI1 into 9DWV CRBN frame",
            "original_batch_status": row["batch_status"], "original_started_at": row["started_at"],
            "original_finished_at": row["finished_at"], "retry_status": "not rerun", "guided_config": row["guided_config"],
            **stats(models, config),
        })
    combined.extend(retry_rows)
    combined.sort(key=lambda row: float(row["mean_best_10_haddock_score"]))
    for rank, row in enumerate(combined, start=1):
        row["comparative_haddock_rank"] = str(rank)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    fields = [
        "comparative_haddock_rank", "source", "molecule", "label", "guided_selection_rank", "vina_rank",
        "canonical_smiles", "rf_active_probability", "vina_score_kcal_mol",
        "pocket_centroid_distance_angstrom", "ppil4_min_distance_angstrom", "geometry_status",
        "rf_component", "vina_component", "pocket_component", "ppil4_geometry_component",
        "composite_pre_haddock_score", "max_similarity_to_other_selected", "initial_ligand_pose",
        "original_batch_status", "original_started_at", "original_finished_at", "retry_status", "models_scored",
        "best_haddock_score", "mean_best_10_haddock_score", "sd_best_10_haddock_score",
        "median_haddock_score", "worst_haddock_score", "best_model_structure", "best_model_path", "guided_config",
    ]
    write(OUTPUT / "comparative_ranking_with_fpft2216_baseline.csv", combined, fields)
    write(OUTPUT / "retry_and_baseline_status.csv", status_rows, list(status_rows[0]))
    print(f"Original completed candidates: {len(original)}")
    print(f"Completed corrected retries/baseline: {len(retry_rows)}")
    print(f"Combined table: {OUTPUT / 'comparative_ranking_with_fpft2216_baseline.csv'}")


if __name__ == "__main__":
    main()

