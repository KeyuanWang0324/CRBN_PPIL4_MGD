#!/usr/bin/env python3
"""Collect and rank results from the 30-candidate guided PPIL4 HADDOCK screen.

Completed candidates are ranked by the mean HADDOCK emscoring score of their
best ten rigid-body/refined models (lower/more negative is better *within this
single, identical exploratory protocol*).  This is not a binding-affinity,
selectivity, or experimentally validated PPIL4-recruitment prediction.
"""

from __future__ import annotations

import csv
import statistics
from datetime import datetime
from pathlib import Path


ROOT = Path("haddock3_ternary/guided_ppil4_9DWV")
TRIAGE = ROOT / "pre_haddock_triage"
SHORTLIST = TRIAGE / "guided_haddock_shortlist_30.csv"
BATCH_STATUS = TRIAGE / "guided_haddock_batch_status.csv"
OUTPUT = TRIAGE / "haddock_results"


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def elapsed_minutes(started: str, finished: str) -> str:
    if not started or not finished:
        return ""
    try:
        delta = datetime.fromisoformat(finished) - datetime.fromisoformat(started)
    except ValueError:
        return ""
    return f"{delta.total_seconds() / 60:.2f}"


def load_scores(config: Path) -> list[dict[str, object]]:
    table = config.parent / "run/3_emscoring/emscoring.tsv"
    if not table.is_file():
        return []
    with table.open(newline="") as handle:
        return [{**row, "score": float(row["score"])} for row in csv.DictReader(handle, delimiter="\t")]


def main() -> None:
    for required in (SHORTLIST, BATCH_STATUS):
        if not required.is_file():
            raise FileNotFoundError(required)
    with SHORTLIST.open(newline="") as handle:
        shortlist = list(csv.DictReader(handle))
    with BATCH_STATUS.open(newline="") as handle:
        status_by_selection = {row["guided_selection_rank"]: row for row in csv.DictReader(handle)}

    OUTPUT.mkdir(parents=True, exist_ok=True)
    completed, campaign, raw_models = [], [], []
    for candidate in shortlist:
        selection_rank = candidate["guided_selection_rank"]
        status = status_by_selection.get(selection_rank, {})
        config = Path(candidate["guided_config"])
        models = load_scores(config)
        base = {
            "guided_selection_rank": selection_rank,
            "vina_rank": candidate["vina_rank"],
            "molecule": candidate["molecule"],
            "canonical_smiles": candidate["canonical_smiles"],
            "rf_active_probability": candidate["rf_active_probability"],
            "vina_score_kcal_mol": candidate["vina_score_kcal_mol"],
            "pocket_centroid_distance_angstrom": candidate["pocket_centroid_distance_angstrom"],
            "ppil4_min_distance_angstrom": candidate["ppil4_min_distance_angstrom"],
            "composite_pre_haddock_score": candidate["composite_pre_haddock_score"],
            "max_similarity_to_other_selected": candidate["max_similarity_to_other_selected"],
            "batch_status": status.get("status", "not_started"),
            "started_at": status.get("started_at", ""),
            "finished_at": status.get("finished_at", ""),
            "elapsed_minutes": elapsed_minutes(status.get("started_at", ""), status.get("finished_at", "")),
            "return_code": status.get("return_code", ""),
            "guided_config": str(config),
            "run_directory": str(config.parent / "run"),
            "launcher_log": status.get("log_file", str(config.parent / "guided_haddock_launcher.log")),
        }
        if models:
            models.sort(key=lambda item: float(item["score"]))
            top = models[: min(10, len(models))]
            summary = {
                "models_scored": len(models),
                "best_haddock_score": f"{models[0]['score']:.3f}",
                "mean_best_10_haddock_score": f"{statistics.mean(float(row['score']) for row in top):.3f}",
                "sd_best_10_haddock_score": f"{statistics.stdev(float(row['score']) for row in top):.3f}" if len(top) > 1 else "0.000",
                "median_haddock_score": f"{statistics.median(float(row['score']) for row in models):.3f}",
                "worst_haddock_score": f"{models[-1]['score']:.3f}",
            }
            completed.append({**base, **summary})
            for model_rank, model in enumerate(models, start=1):
                raw_models.append({
                    "guided_selection_rank": selection_rank,
                    "molecule": candidate["molecule"], "vina_rank": candidate["vina_rank"],
                    "model_score_rank": model_rank, "structure": model["structure"],
                    "original_name": model["original_name"], "haddock_score": f"{model['score']:.3f}",
                    "model_path": str(config.parent / "run/3_emscoring" / model["structure"]),
                })
            campaign.append({**base, **summary})
        else:
            campaign.append({**base, "models_scored": 0, "best_haddock_score": "", "mean_best_10_haddock_score": "", "sd_best_10_haddock_score": "", "median_haddock_score": "", "worst_haddock_score": ""})

    completed.sort(key=lambda row: float(row["mean_best_10_haddock_score"]))
    for rank, row in enumerate(completed, start=1):
        row["haddock_rank"] = rank
    campaign_by_key = {row["guided_selection_rank"]: row for row in completed}
    for row in campaign:
        row["haddock_rank"] = campaign_by_key.get(row["guided_selection_rank"], {}).get("haddock_rank", "")

    shared = [
        "haddock_rank", "guided_selection_rank", "vina_rank", "molecule", "canonical_smiles",
        "rf_active_probability", "vina_score_kcal_mol", "pocket_centroid_distance_angstrom",
        "ppil4_min_distance_angstrom", "composite_pre_haddock_score",
        "max_similarity_to_other_selected", "models_scored", "best_haddock_score",
        "mean_best_10_haddock_score", "sd_best_10_haddock_score", "median_haddock_score",
        "worst_haddock_score", "batch_status", "started_at", "finished_at", "elapsed_minutes",
        "return_code", "guided_config", "run_directory", "launcher_log",
    ]
    write_csv(OUTPUT / "completed_candidates_ranked.csv", completed, shared)
    write_csv(OUTPUT / "all_30_campaign_status.csv", campaign, shared)
    write_csv(OUTPUT / "all_completed_model_scores.csv", raw_models, [
        "guided_selection_rank", "molecule", "vina_rank", "model_score_rank", "structure",
        "original_name", "haddock_score", "model_path",
    ])
    print(f"Collected {len(completed)} completed candidates and {len(raw_models)} scored HADDOCK models.")
    print(f"Ranked results: {OUTPUT / 'completed_candidates_ranked.csv'}")
    print(f"All 30 statuses: {OUTPUT / 'all_30_campaign_status.csv'}")
    print(f"All raw model scores: {OUTPUT / 'all_completed_model_scores.csv'}")


if __name__ == "__main__":
    main()

