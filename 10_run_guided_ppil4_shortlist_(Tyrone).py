#!/usr/bin/env python3
"""Run the 30-candidate 9DWV-guided PPIL4 HADDOCK3 shortlist safely.

Runs one candidate at a time, so --cores 8 uses the Mac's eight cores for the
current candidate only.  Completed runs are skipped on later invocations.
This is an exploratory, 9DWV-guided protocol and not evidence of PPIL4
recruitment or experimental binding affinity.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


ROOT = Path("haddock3_ternary/guided_ppil4_9DWV")
SHORTLIST = ROOT / "pre_haddock_triage/guided_haddock_shortlist_30.csv"
STATUS_FILE = ROOT / "pre_haddock_triage/guided_haddock_batch_status.csv"


def is_complete(config: Path) -> bool:
    """A normal guided protocol finishes by writing the emscoring module IO."""
    run_dir = config.parent / "run"
    return (run_dir / "3_emscoring/io.json").is_file()


def set_ncores(config: Path, cores: int) -> None:
    text = config.read_text()
    updated, count = re.subn(r"(?m)^ncores\s*=\s*\d+\s*$", f"ncores = {cores}", text)
    if count != 1:
        raise RuntimeError(f"Expected exactly one ncores line in {config}, found {count}.")
    if updated != text:
        config.write_text(updated)


def write_status(rows: list[dict[str, str]]) -> None:
    STATUS_FILE.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "guided_selection_rank", "vina_rank", "molecule", "status",
        "started_at", "finished_at", "return_code", "guided_config", "log_file",
    ]
    with STATUS_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def haddock_executable() -> str:
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate = Path(conda_prefix) / "bin/haddock3"
        if candidate.is_file():
            return str(candidate)
    fallback = shutil.which("haddock3")
    if fallback:
        return fallback
    raise RuntimeError(
        "HADDOCK3 was not found. Activate my-rdkit-env first: conda activate my-rdkit-env"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cores", type=int, default=8, help="Cores for one sequential HADDOCK job (default: 8).")
    parser.add_argument("--start", type=int, default=1, help="First shortlist selection rank to run.")
    parser.add_argument("--end", type=int, default=30, help="Last shortlist selection rank to run.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned jobs without starting HADDOCK3.")
    args = parser.parse_args()
    if args.cores < 1 or args.start < 1 or args.end < args.start:
        raise SystemExit("--cores must be positive and require 1 <= --start <= --end.")
    if not SHORTLIST.is_file():
        raise FileNotFoundError(f"Shortlist not found: {SHORTLIST}. Run reduce_guided_ppil4_candidates.py first.")

    with SHORTLIST.open(newline="") as handle:
        shortlist = list(csv.DictReader(handle))
    selected = [row for row in shortlist if args.start <= int(row["guided_selection_rank"]) <= args.end]
    if not selected:
        raise SystemExit("No shortlist candidates matched that selection-rank range.")
    executable = haddock_executable()
    statuses: list[dict[str, str]] = []
    print(f"Processing {len(selected)} candidates sequentially with {args.cores} core(s) each.")

    for number, row in enumerate(selected, start=1):
        config = Path(row["guided_config"])
        if not config.is_file():
            status = "missing_config"
            print(f"[{number}/{len(selected)}] molecule {row['molecule']}: {status}")
            statuses.append({**row, "status": status, "started_at": "", "finished_at": "", "return_code": "", "guided_config": str(config), "log_file": ""})
            write_status(statuses)
            continue
        log_file = config.parent / "guided_haddock_launcher.log"
        run_dir = config.parent / "run"
        if is_complete(config):
            status = "skipped_complete"
            print(f"[{number}/{len(selected)}] molecule {row['molecule']}: already complete; skipped.")
            statuses.append({**row, "status": status, "started_at": "", "finished_at": "", "return_code": "0", "guided_config": str(config), "log_file": str(log_file)})
            write_status(statuses)
            continue
        if run_dir.exists():
            # Restarting from a partially written HADDOCK directory can mix
            # old and new module output. Preserve it for diagnosis and let
            # the batch continue with clean candidates.
            status = "incomplete_needs_review"
            print(f"[{number}/{len(selected)}] molecule {row['molecule']}: existing incomplete run; skipped safely.")
            statuses.append({**row, "status": status, "started_at": "", "finished_at": "", "return_code": "", "guided_config": str(config), "log_file": str(log_file)})
            write_status(statuses)
            continue

        set_ncores(config, args.cores)
        command = [executable, str(config), "--restart", "0"]
        if args.dry_run:
            print(f"[{number}/{len(selected)}] would run molecule {row['molecule']}: {' '.join(command)}")
            statuses.append({**row, "status": "dry_run", "started_at": "", "finished_at": "", "return_code": "", "guided_config": str(config), "log_file": str(log_file)})
            write_status(statuses)
            continue

        started = datetime.now().isoformat(timespec="seconds")
        print(f"[{number}/{len(selected)}] starting molecule {row['molecule']} (shortlist rank {row['guided_selection_rank']}).")
        with log_file.open("a") as log:
            log.write(f"\n\n=== {started}: {' '.join(command)} ===\n")
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        finished = datetime.now().isoformat(timespec="seconds")
        status = "complete" if result.returncode == 0 and is_complete(config) else "failed"
        print(f"[{number}/{len(selected)}] molecule {row['molecule']}: {status} (exit {result.returncode}).")
        statuses.append({**row, "status": status, "started_at": started, "finished_at": finished, "return_code": str(result.returncode), "guided_config": str(config), "log_file": str(log_file)})
        write_status(statuses)

    complete = sum(row["status"] in {"complete", "skipped_complete"} for row in statuses)
    failed = sum(row["status"] in {"failed", "missing_config", "incomplete_needs_review"} for row in statuses)
    print(f"Finished launcher: {complete} complete/skipped, {failed} failed. Status table: {STATUS_FILE}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()

