#!/usr/bin/env python3
"""Aggregate one molecule from the final post-emref caprieval table."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


SCHEMA = (
    "molecule",
    "mean_best_10",
    "best_haddock",
    "sd_best_10",
    "dockq_vs_9dwv",
    "cluster_id",
    "n_models",
    "protocol_version",
)


def read_capri_rows(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    if not lines:
        raise ValueError(f"No tabular rows found in {path}")
    rows = list(csv.DictReader(lines, delimiter="\t"))
    required = {"model", "score"}
    missing = required.difference(rows[0] if rows else {})
    if missing:
        raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")
    return rows


def finite_float(value: str | None) -> float | None:
    try:
        number = float(value) if value not in (None, "", "-") else None
    except (TypeError, ValueError):
        return None
    return number if number is not None and math.isfinite(number) else None


def aggregate(molecule: str, rows: list[dict[str, str]], protocol_version: str) -> dict[str, object]:
    successful: list[tuple[float, str, dict[str, str]]] = []
    for row in rows:
        score = finite_float(row.get("score"))
        model = row.get("model", "")
        if score is not None and model:
            successful.append((score, model, row))

    successful.sort(key=lambda item: (item[0], item[1]))
    n_models = len(successful)
    if n_models < 10:
        return {
            "molecule": molecule,
            "mean_best_10": "",
            "best_haddock": "",
            "sd_best_10": "",
            "dockq_vs_9dwv": "",
            "cluster_id": "",
            "n_models": n_models,
            "protocol_version": protocol_version,
        }

    top = successful[:10]
    scores = [item[0] for item in top]
    best_row = top[0][2]
    dockq = finite_float(best_row.get("dockq"))
    cluster_id = best_row.get("cluster_id", "")
    if cluster_id == "-":
        cluster_id = "unclustered"
    return {
        "molecule": molecule,
        "mean_best_10": f"{math.fsum(scores) / 10:.6f}",
        "best_haddock": f"{scores[0]:.6f}",
        "sd_best_10": f"{statistics.pstdev(scores):.6f}",
        "dockq_vs_9dwv": "" if dockq is None else f"{dockq:.6f}",
        "cluster_id": cluster_id,
        "n_models": n_models,
        "protocol_version": protocol_version,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("molecule")
    parser.add_argument("capri_ss_tsv", type=Path, help="Final post-emref caprieval capri_ss.tsv")
    parser.add_argument("--output", type=Path, help="Write the one-row CSV here; defaults to stdout")
    parser.add_argument(
        "--protocol-version-file",
        type=Path,
        default=Path(__file__).with_name("PROTOCOL_VERSION.txt"),
    )
    args = parser.parse_args()

    protocol_version = args.protocol_version_file.read_text().strip()
    result = aggregate(args.molecule, read_capri_rows(args.capri_ss_tsv), protocol_version)

    if args.output:
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SCHEMA)
            writer.writeheader()
            writer.writerow(result)
    else:
        import sys

        writer = csv.DictWriter(sys.stdout, fieldnames=SCHEMA)
        writer.writeheader()
        writer.writerow(result)


if __name__ == "__main__":
    main()
