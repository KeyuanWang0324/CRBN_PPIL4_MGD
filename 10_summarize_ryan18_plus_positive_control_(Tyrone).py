#!/usr/bin/env python3
"""Write a detailed 19-entry ranked table for the Ryan/friend corrected campaign."""
from __future__ import annotations
import csv
import statistics
from pathlib import Path

ROOT = Path("haddock3_ternary/ryan18_friend19_corrected")
MANIFEST = ROOT / "ryan18_friend19_manifest.csv"
OUTPUT = ROOT / "results"

def models(config: Path):
    table = config.parent / "run/3_emscoring/emscoring.tsv"
    if not table.is_file(): return []
    with table.open(newline="") as f: rows = list(csv.DictReader(f, delimiter="\t"))
    return sorted(rows, key=lambda row: float(row["score"]))

def stats(rows, config: Path):
    values=[float(row["score"]) for row in rows]; top=values[:10]
    return {"models_scored":str(len(values)), "best_haddock_score":f"{values[0]:.3f}",
            "mean_best_10_haddock_score":f"{statistics.mean(top):.3f}",
            "sd_best_10_haddock_score":f"{statistics.stdev(top):.3f}" if len(top)>1 else "0.000",
            "median_haddock_score":f"{statistics.median(values):.3f}", "worst_haddock_score":f"{values[-1]:.3f}",
            "best_model_structure":rows[0]["structure"], "best_model_path":str(config.parent / "run/3_emscoring" / rows[0]["structure"])}

def write(path, rows, fields):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def main():
    if not MANIFEST.is_file(): raise FileNotFoundError(f"Run run_ryan18_friend19_corrected.py first: {MANIFEST}")
    with MANIFEST.open(newline="") as f: jobs=list(csv.DictReader(f))
    complete=[]; incomplete=[]
    for job in jobs:
        config=Path(job["guided_config"]); rows=models(config)
        if not rows:
            incomplete.append({"molecule":job["molecule"],"label":job["label"],"guided_config":str(config)}); continue
        complete.append({**job, "source":"Ryan matched candidate" if job["kind"]=="ryan_candidate" else "positive control", **stats(rows,config)})
    complete.sort(key=lambda row:float(row["mean_best_10_haddock_score"]))
    for rank,row in enumerate(complete,1): row["final_corrected_haddock_rank"]=str(rank)
    OUTPUT.mkdir(parents=True,exist_ok=True)
    fields=["final_corrected_haddock_rank","source","kind","molecule","label","smiles","canonical_smiles","rf_active_probability",
            "vina_score_kcal_mol","pocket_centroid_distance_angstrom","vina_box_reference","vina_pose_file","models_scored",
            "best_haddock_score","mean_best_10_haddock_score","sd_best_10_haddock_score","median_haddock_score","worst_haddock_score",
            "best_model_structure","best_model_path","guided_config"]
    write(OUTPUT / "ryan18_friend19_final_corrected_ranked.csv",complete,fields)
    write(OUTPUT / "incomplete_runs.csv",incomplete,["molecule","label","guided_config"])
    print(f"Completed entries: {len(complete)} of 19 expected.")
    print(f"Results: {OUTPUT / 'ryan18_friend19_final_corrected_ranked.csv'}")

if __name__=="__main__": main()

