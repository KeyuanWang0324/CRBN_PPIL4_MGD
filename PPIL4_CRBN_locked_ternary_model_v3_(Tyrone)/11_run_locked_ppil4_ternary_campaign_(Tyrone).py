#!/usr/bin/env python3
"""Prepare and resume the locked 30-candidate PPIL4 HADDOCK campaign."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import tempfile


LOCK = Path(__file__).resolve().parent
ROOT = LOCK.parent
QUEUE = LOCK / "candidates_30_final_corrected_ranked.csv"
HADDOCK_BIN = Path("/opt/anaconda3/envs/my-rdkit-env/bin/haddock3")
PROTOCOL_VERSION = (LOCK / "PROTOCOL_VERSION.txt").read_text().strip()
MIN_FREE_BYTES = 3 * 1024**3
RESULT_FIELDS = [
    "locked_candidate_rank", "molecule", "mean_best_10", "best_haddock",
    "sd_best_10", "dockq_vs_9dwv", "cluster_id", "n_models",
    "protocol_version", "status",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_queue() -> list[dict[str, str]]:
    with QUEUE.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 30:
        raise RuntimeError(f"Expected 30 locked candidates, found {len(rows)}")
    return rows


def find_pose(row: dict[str, str]) -> Path:
    molecule = row["molecule"]
    packaged = LOCK / "candidate_ligand_prep" / f"candidate_{molecule}.pdb"
    if packaged.exists():
        return packaged
    if "halogen" in row["source"]:
        pattern = (
            f"haddock3_ternary/guided_ppil4_9DWV_halogen_retry_baseline/"
            f"*molecule_{molecule}/*molecule_{molecule}_corrected_halogen.pdb"
        )
    else:
        pattern = (
            f"haddock3_ternary/guided_ppil4_9DWV_corrected_full_calibration/"
            f"*molecule_{molecule}/*molecule_{molecule}_corrected.pdb"
        )
    hits = sorted(ROOT.glob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"Expected one corrected pose for molecule {molecule}, found {hits}")
    return hits[0]


def prepare_candidate(row: dict[str, str], output_root: Path) -> Path:
    rank = int(row["locked_candidate_rank"])
    molecule = row["molecule"]
    candidate = f"candidate_{molecule}"
    work = output_root / f"{rank:02d}_{candidate}"
    inputs = work / "inputs"
    ligand_prep = work / "ligand_prep"
    inputs.mkdir(parents=True, exist_ok=True)
    ligand_prep.mkdir(parents=True, exist_ok=True)

    for name in (
        "CRBN_receptor_thalidomide_Ryan.pdb",
        "PPIL4_chainB.pdb",
        "reference_9DWV_CRBN_A_PPIL4_B_FPFT_C.pdb",
    ):
        shutil.copy2(LOCK / "inputs" / name, inputs / name)
    shutil.copy2(LOCK / "ambig_FIXED.tbl", work / "ambig_FIXED.tbl")

    source_pose = find_pose(row)
    packaged_files = [
        LOCK / "candidate_ligand_prep" / f"{candidate}.pdb",
        LOCK / "candidate_ligand_prep" / f"{candidate}_prodrg_fixed.pdb",
        LOCK / "candidate_ligand_prep" / f"{candidate}_prodrg.top",
        LOCK / "candidate_ligand_prep" / f"{candidate}_prodrg.param",
    ]
    missing = [path.name for path in packaged_files if not path.exists()]
    if missing:
        raise RuntimeError(f"Missing packaged ligand files for {candidate}: {missing}")
    for source in packaged_files:
        shutil.copy2(source, ligand_prep / source.name)

    template = (LOCK / "guided_ppil4_haddock.cfg").read_text()
    config = template.replace("<candidate>", candidate)
    cfg_path = work / f"{candidate}.cfg"
    cfg_path.write_text(config)
    metadata = {
        "locked_candidate_rank": rank,
        "molecule": molecule,
        "candidate": candidate,
        "smiles": row["canonical_smiles"],
        "source_pose": str(source_pose.relative_to(LOCK)),
        "source_pose_sha256": sha256(source_pose),
        "protocol_version": PROTOCOL_VERSION,
        "status": "prepared",
    }
    (work / "candidate_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return work


def free_bytes(path: Path) -> int:
    return shutil.disk_usage(path).free


def read_results(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {row["molecule"]: row for row in csv.DictReader(handle)}


def write_results(path: Path, rows: dict[str, dict[str, object]]) -> None:
    ordered = sorted(rows.values(), key=lambda row: int(row["locked_candidate_rank"]))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(ordered)


def aggregate_capri(molecule: str, capri_path: Path) -> dict[str, object]:
    lines = [line for line in capri_path.read_text().splitlines() if line and not line.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    valid = []
    for row in rows:
        try:
            score = float(row["score"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(score) and row.get("model"):
            valid.append((score, row["model"], row))
    valid.sort(key=lambda item: (item[0], item[1]))
    if len(valid) < 10:
        return {
            "molecule": molecule, "mean_best_10": "", "best_haddock": "",
            "sd_best_10": "", "dockq_vs_9dwv": "", "cluster_id": "",
            "n_models": len(valid), "protocol_version": PROTOCOL_VERSION,
            "status": "insufficient_models",
        }
    top = valid[:10]
    scores = [item[0] for item in top]
    best = top[0][2]
    try:
        dockq = float(best.get("dockq", ""))
        dockq_text = f"{dockq:.6f}" if math.isfinite(dockq) else ""
    except ValueError:
        dockq_text = ""
    cluster = best.get("cluster_id", "")
    if cluster == "-":
        cluster = "unclustered"
    return {
        "molecule": molecule,
        "mean_best_10": f"{math.fsum(scores) / 10:.6f}",
        "best_haddock": f"{scores[0]:.6f}",
        "sd_best_10": f"{statistics.pstdev(scores):.6f}",
        "dockq_vs_9dwv": dockq_text,
        "cluster_id": cluster,
        "n_models": len(valid),
        "protocol_version": PROTOCOL_VERSION,
        "status": "complete",
    }


def fpft_gate_open() -> bool:
    gate = LOCK / "REFERENCE_FPFT.json"
    if not gate.exists():
        return False
    data = json.loads(gate.read_text())
    return (
        data.get("protocol_version") == PROTOCOL_VERSION
        and data.get("independently_reproduced") is True
    )


def run_queue(output_root: Path, max_candidates: int | None) -> None:
    if not fpft_gate_open():
        raise RuntimeError(
            "FPFT gate is closed. Create REFERENCE_FPFT.json in the model folder "
            "with this protocol version and independently_reproduced=true only after Ryan reproduces it."
        )
    results_path = output_root / "campaign_results.csv"
    results = read_results(results_path)
    launched = 0
    for row in load_queue():
        molecule = row["molecule"]
        if results.get(molecule, {}).get("status") == "complete":
            continue
        if max_candidates is not None and launched >= max_candidates:
            break
        if free_bytes(output_root) < MIN_FREE_BYTES:
            print("Stopping safely: less than 3 GiB free.", file=sys.stderr)
            break
        rank = int(row["locked_candidate_rank"])
        candidate = f"candidate_{molecule}"
        work = output_root / f"{rank:02d}_{candidate}"
        cfg = work / f"{candidate}.cfg"
        if not cfg.exists():
            raise RuntimeError(f"Candidate is not prepared: {cfg}")
        log = work / "haddock3.log"
        with log.open("a") as handle:
            subprocess.run([str(HADDOCK_BIN), cfg.name], cwd=work, stdout=handle, stderr=subprocess.STDOUT, check=True)
        capri = work / "run_ppil4_ternary" / "9_caprieval" / "capri_ss.tsv"
        result = aggregate_capri(molecule, capri)
        result["locked_candidate_rank"] = rank
        results[molecule] = result
        write_results(results_path, results)
        subprocess.run([str(HADDOCK_BIN) + "-clean", "run_ppil4_ternary"], cwd=work, check=False)
        launched += 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "locked_campaign_runs" / "ppil4_crbn_ternary_v3_Tyrone",
    )
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--max-candidates", type=int)
    args = parser.parse_args()
    if not args.prepare and not args.run:
        parser.error("Select --prepare and/or --run")
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    if args.prepare:
        for row in load_queue():
            work = prepare_candidate(row, output_root)
            print(f"prepared {work.name}")
    if args.run:
        run_queue(output_root, args.max_candidates)


if __name__ == "__main__":
    main()
