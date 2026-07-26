#!/usr/bin/env python3
"""Reduce the 250 guided PPIL4 candidates before expensive HADDOCK runs.

This is a reproducible, pre-HADDOCK triage—not a binding-affinity predictor.
It combines RF CRBN probability, Vina score, 9DWV pocket proximity, PPIL4-RRM
geometry (clash/proximity), and RDKit fingerprint diversity.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


GUIDED_ROOT = Path("haddock3_ternary/guided_ppil4_9DWV")
GUIDED_MANIFEST = GUIDED_ROOT / "guided_candidate_manifest.csv"
VINA_RESULTS = Path("crbn_vina_rf_top500/top_250_vina_docked.csv")
PPIL4 = Path("haddock3_ternary/validation_9DWV/bound_control/ppil4_rrm_9DWV_C.pdb")
OUTPUT_DIR = GUIDED_ROOT / "pre_haddock_triage"


def pdb_heavy_coordinates(path: Path) -> np.ndarray:
    points = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[76:78].strip().upper() == "H":
            continue
        points.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not points:
        raise RuntimeError(f"No heavy atoms found in {path}")
    return np.asarray(points)


def min_distance(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.min(np.linalg.norm(first[:, None, :] - second[None, :, :], axis=2)))


def scaled(values: list[float], reverse: bool = False) -> list[float]:
    low, high = min(values), max(values)
    if math.isclose(low, high):
        return [1.0] * len(values)
    result = [(value - low) / (high - low) for value in values]
    return [1.0 - value for value in result] if reverse else result


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path.name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefilter-n", type=int, default=75)
    parser.add_argument("--select-n", type=int, default=30)
    parser.add_argument("--max-pocket-distance", type=float, default=5.0)
    parser.add_argument("--min-ppil4-distance", type=float, default=2.0, help="Reject steric clashes below this distance (Å).")
    parser.add_argument("--max-ppil4-distance", type=float, default=6.0, help="Require a ligand-facing PPIL4 surface within this distance (Å).")
    args = parser.parse_args()
    if args.prefilter_n < args.select_n or args.select_n < 1:
        raise SystemExit("Require --prefilter-n >= --select-n >= 1.")
    for required in (GUIDED_MANIFEST, VINA_RESULTS, PPIL4):
        if not required.is_file():
            raise FileNotFoundError(required)

    with GUIDED_MANIFEST.open(newline="") as handle:
        guided = list(csv.DictReader(handle))
    with VINA_RESULTS.open(newline="") as handle:
        vina = {row["molecule"]: row for row in csv.DictReader(handle)}
    ppil4_points = pdb_heavy_coordinates(PPIL4)
    rows = []
    for row in guided:
        molecule = row["molecule"]
        if molecule not in vina or row["status"] != "prepared":
            continue
        ligand_file = Path(row["config"]).parent / f"rank{int(row['vina_rank']):03d}_molecule_{molecule}_aligned_to_9DWV.pdb"
        ligand_points = pdb_heavy_coordinates(ligand_file)
        ppil4_distance = min_distance(ligand_points, ppil4_points)
        pocket_distance = float(row["pocket_centroid_distance_angstrom"])
        geometry_status = "eligible"
        if pocket_distance > args.max_pocket_distance:
            geometry_status = "pocket_incompatible"
        elif ppil4_distance < args.min_ppil4_distance:
            geometry_status = "ppil4_clash"
        elif ppil4_distance > args.max_ppil4_distance:
            geometry_status = "ppil4_too_distant"
        rows.append({
            "vina_rank": int(row["vina_rank"]), "molecule": molecule,
            "canonical_smiles": row["canonical_smiles"],
            "rf_active_probability": float(vina[molecule]["rf_active_probability"]),
            "vina_score_kcal_mol": float(row["vina_score_kcal_mol"]),
            "pocket_centroid_distance_angstrom": pocket_distance,
            "ppil4_min_distance_angstrom": ppil4_distance,
            "geometry_status": geometry_status, "guided_config": row["config"],
        })

    eligible = [row for row in rows if row["geometry_status"] == "eligible"]
    if len(eligible) < args.select_n:
        raise RuntimeError(f"Only {len(eligible)} candidates pass geometry filters; lower thresholds or select fewer.")
    rf_scaled = scaled([row["rf_active_probability"] for row in eligible])
    vina_scaled = scaled([row["vina_score_kcal_mol"] for row in eligible], reverse=True)
    pocket_scaled = scaled([row["pocket_centroid_distance_angstrom"] for row in eligible], reverse=True)
    # Ideal ligand-to-PPIL4 clearance is near 3.5 Å: close enough to mediate
    # contacts, but not an atomic clash.
    interface_scaled = [math.exp(-((row["ppil4_min_distance_angstrom"] - 3.5) ** 2) / (2 * 1.25**2)) for row in eligible]
    for row, rf, vina_score, pocket, interface in zip(eligible, rf_scaled, vina_scaled, pocket_scaled, interface_scaled):
        row["rf_component"] = rf
        row["vina_component"] = vina_score
        row["pocket_component"] = pocket
        row["ppil4_geometry_component"] = interface
        row["composite_pre_haddock_score"] = 0.35 * rf + 0.30 * vina_score + 0.20 * pocket + 0.15 * interface
    for row in rows:
        if row["geometry_status"] != "eligible":
            row.update(rf_component="", vina_component="", pocket_component="", ppil4_geometry_component="", composite_pre_haddock_score="")

    rows.sort(key=lambda row: float(row["composite_pre_haddock_score"]) if row["geometry_status"] == "eligible" else -1.0, reverse=True)
    # A “top 75” file means *up to* 75 viable candidates.  Do not pad it with
    # geometry-rejected molecules: that would make it unsafe to use directly
    # as the input set for a costly HADDOCK campaign.
    prefilter = eligible[:args.prefilter_n]
    candidates = prefilter
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = []
    for row in candidates:
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None:
            raise RuntimeError(f"Invalid SMILES for molecule {row['molecule']}")
        fingerprints.append(generator.GetFingerprint(molecule))

    # Greedy score/diversity selection: retain high composite-score candidates
    # while avoiding near-duplicate halogen/methyl analogue series.
    selected_indices = [0]
    while len(selected_indices) < min(args.select_n, len(candidates)):
        best_index, best_value = None, -float("inf")
        for index, row in enumerate(candidates):
            if index in selected_indices:
                continue
            max_similarity = max(DataStructs.TanimotoSimilarity(fingerprints[index], fingerprints[chosen]) for chosen in selected_indices)
            value = 0.70 * float(row["composite_pre_haddock_score"]) + 0.30 * (1.0 - max_similarity)
            if value > best_value:
                best_index, best_value = index, value
        selected_indices.append(best_index)
    shortlist = []
    for selection_rank, index in enumerate(selected_indices, start=1):
        row = dict(candidates[index])
        max_similarity = max((DataStructs.TanimotoSimilarity(fingerprints[index], fingerprints[chosen]) for chosen in selected_indices if chosen != index), default=0.0)
        row["guided_selection_rank"] = selection_rank
        row["max_similarity_to_other_selected"] = max_similarity
        shortlist.append(row)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "all_250_scored.csv", rows)
    write_csv(OUTPUT_DIR / f"top_{args.prefilter_n}_pre_haddock.csv", prefilter)
    write_csv(OUTPUT_DIR / f"guided_haddock_shortlist_{args.select_n}.csv", shortlist)
    print(f"Scored {len(rows)} candidates; {len(eligible)} pass 9DWV pocket/PPIL4 geometry filters.")
    print(f"Wrote {len(prefilter)} eligible prefilter candidates and {len(shortlist)} diverse guided-HADDOCK candidates to {OUTPUT_DIR}.")


if __name__ == "__main__":
    main()

