"""
Positive control: does this project's docking METHODOLOGY (04's Vina pose
-> 08's ligand-inclusive 3-body HADDOCK3 protocol) actually reproduce a
KNOWN real CRBN-PPIL4 ternary complex, when given the real drug that
complex was solved with?

This is NOT a per-candidate scoring loop anymore (see chat discussion for
why the earlier version of this script -- which scored every one of 06's
candidates against the real 9DWV structure -- was methodologically
invalid): 05/06/08's PPIL4 restraint (ppil4_pocket_residues() /
[249, 250, 273, 275, 276, 277, 278, 279]) is ITSELF derived from 9DWV's
own real contact residues. Scoring a candidate's dockq against 9DWV after
handing its restraints the answer key isn't independent validation --
PPIL4 lands close to the real interface almost by construction (see
00_validate_docking_interface_05_(Ryan).py: 100% "found the spot" on
PPIL4 across the whole 05 population, mean 6.3/8 residues). A high
dockq-vs-9DWV per candidate was never evidence that candidate is a good
binder -- it mostly just confirmed the restraint was applied. Candidate
ranking now comes from 08's own HADDOCK score (see that script).

What DOES make sense: dock the real drug -- FPFT-2216, PDB ligand code
A1BC8, SMILES verified against RCSB's official ligand page AND
cross-checked against the actual atom composition in the cached 9DWV
structure (12 C / 4 N / 3 O / 1 S in the structure, matching C12H12N4O3S
exactly -- not invented, not assumed) -- through this project's OWN
pipeline, the same 04 Vina step and 08 ligand-inclusive HADDOCK3 protocol
every library candidate goes through, and see whether THAT reproduces
the real 9DWV structure. This tests the methodology itself, using a case
where the right answer is already known, instead of testing a candidate
whose right answer is exactly what the pipeline is trying to find out.

Setup (do this BEFORE running this script):
  1. In 04_vina_dock_candidates_(Ryan).py, set:
         MANUAL_CANDIDATE = (POSITIVE_CONTROL_NAME, POSITIVE_CONTROL_SMILES)
     (same constants as below) and run it -- docks FPFT-2216 into CRBN's
     thalidomide pocket and PPIL4's RRM domain via Vina, same as any
     library candidate, producing crbn_contacts.txt and
     CRBN_candidate_complex.pdb under docking_tmp/haddock3_novel_candidate/.
  2. In 08_haddock3_ternary_with_ligand_(Ryan).py, set:
         CANDIDATE_NAME = POSITIVE_CONTROL_NAME
     and run it -- the full ligand-inclusive 3-body HADDOCK3 docking
     (~45 min - 1.5 hr), producing a rank-1 ternary model.
  3. Run this script to score that model against the real 9DWV structure.

No re-docking is done here -- like the old version of this script, this
only re-evaluates the rank-1 model 08 already produced, scoring it with
HADDOCK3's own CAPRI class (haddock.libs.libcapri.CAPRI) against the real
FPFT-2216/9DWV reference. 08's model has 3 chains (A=CRBN, C=ligand,
B=PPIL4); the ligand chain is stripped before scoring, since the
reference (and dockq/irmsd/fnat/lrmsd) are about the CRBN-PPIL4
protein-protein interface specifically -- matching how the reference
itself was built (DDB1/ligand/Zn dropped, see build_reference_pdb()).

Outputs (this script used to only print, which is why 08's dockq_vs_9dwv
column stayed empty even after 09 had been run):
  - fills dockq_vs_9dwv on the FPFT_2216_positive_control row of
    08_controls_results_(Ryan).csv -- the one molecule 09 scores. The 18
    candidates' dockq_vs_9dwv stays blank on purpose (see above: circular).
  - writes REFERENCE_FPFT.json, PROTOCOL_LOCK.md §5's golden reference.

Run with the same environment 08 uses (needs haddock3 + biopython installed):
    python3 "09_score_vs_fpft2216_reference_(Ryan).py"
"""
import csv
import gzip
import json
import os
import shutil
import statistics
import sys
import time
import urllib.request

SCRIPT_START_TIME = time.time()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Must match 04's MANUAL_CANDIDATE / 08's CANDIDATE_NAME exactly.
POSITIVE_CONTROL_NAME = "FPFT_2216_positive_control"
# RCSB ligand A1BC8's official canonical SMILES -- cross-checked against
# the actual atoms in the cached 9DWV structure (see module docstring).
POSITIVE_CONTROL_SMILES = "COc1cscc1c2cn(nn2)[C@H]3CCC(=O)NC3=O"

RUN_DIR_BASE = os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run")
SCRATCH_DIR = os.path.join(SCRIPT_DIR, "docking_tmp", "09_capri_vs_fpft2216_scratch")

REFERENCE_DIR = os.path.join(SCRIPT_DIR, "reference_structures")
CIF_CACHE_PATH = os.path.join(REFERENCE_DIR, "FPFT-2216_9DWV.cif")
REFERENCE_PDB_PATH = os.path.join(REFERENCE_DIR, "FPFT-2216_9DWV_reference_(Ryan).pdb")

# Where this script's result goes. 08 writes the positive control's row with
# dockq_vs_9dwv deliberately blank (its own caprieval dockq is self-referential
# -- reference = the run's own best model), and its locked_stats() docstring
# says "09 fills the real value against 9DWV". This is that fill: the ONLY
# molecule 09 scores is POSITIVE_CONTROL_NAME, which lives in the controls CSV.
# The 18 candidates' dockq_vs_9dwv stays blank on purpose -- see this module's
# docstring on why a per-candidate dockq vs 9DWV is circular (their PPIL4
# restraints are derived from 9DWV's own contact residues).
CONTROLS_CSV = os.path.join(SCRIPT_DIR, "08_controls_results_(Ryan).csv")
# PROTOCOL_LOCK.md §5: the golden reference both wrappers must reproduce before
# candidate runs are trusted. Committed alongside the locked protocol bundle.
REFERENCE_JSON_PATH = os.path.join(SCRIPT_DIR, "REFERENCE_FPFT.json")
PROTOCOL_VERSION = "ppil4_lock_v1"

# FPFT-2216's structure is deposited under PDB accession 9DWV -- that's
# what RCSB's download URL and the mmCIF's own internal ID use, even
# though we refer to it by compound name everywhere else in this script.
FPFT2216_PDB_ID = "9DWV"
FPFT2216_CIF_URL = f"https://files.rcsb.org/download/{FPFT2216_PDB_ID}.cif"

# Confirmed via RCSB's mmCIF _atom_site.auth_asym_id: chain A = DDB1,
# chain B = CRBN, chain C = PPIL4 (RRM domain, residues ~240-318).
CRBN_SOURCE_CHAIN = "B"
PPIL4_SOURCE_CHAIN = "C"


def download_cif():
    if os.path.exists(CIF_CACHE_PATH):
        return
    os.makedirs(REFERENCE_DIR, exist_ok=True)
    print(f"Downloading FPFT-2216 structure (PDB {FPFT2216_PDB_ID}, mmCIF) from RCSB -> {CIF_CACHE_PATH} ...")
    import ssl
    import certifi
    ctx = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(FPFT2216_CIF_URL, context=ctx) as resp, open(CIF_CACHE_PATH, "wb") as out:
        shutil.copyfileobj(resp, out)
    print("Download complete.")


def build_reference_pdb():
    """Extract the FPFT-2216 structure's CRBN and PPIL4 protein chains
    only (no DDB1, no ligand, no Zn) and relabel them chain A (CRBN) /
    chain B (PPIL4) to match our own models' chain layout."""
    if os.path.exists(REFERENCE_PDB_PATH):
        return
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("FPFT-2216", CIF_CACHE_PATH)
    model = structure[0]

    crbn = model[CRBN_SOURCE_CHAIN]
    ppil4 = model[PPIL4_SOURCE_CHAIN]
    crbn.detach_parent()
    ppil4.detach_parent()
    crbn.id = "A"
    ppil4.id = "B"

    new_structure = Structure("FPFT-2216_reference")
    new_model = Model(0)
    new_model.add(crbn)
    new_model.add(ppil4)
    new_structure.add(new_model)

    class ProteinOnly(Select):
        def accept_residue(self, residue):
            return residue.id[0] == " "  # drop HETATM (ligand, Zn, waters)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(REFERENCE_PDB_PATH, select=ProteinOnly())
    print(f"Wrote reference structure: {REFERENCE_PDB_PATH} "
          f"(chain A = CRBN, chain B = PPIL4, extracted from the FPFT-2216 "
          f"structure's chains {CRBN_SOURCE_CHAIN}/{PPIL4_SOURCE_CHAIN})")


def find_top_model_path():
    """The rank-1 model from 08's ligand-inclusive run for the positive
    control candidate (same TSV lookup logic as 07/08 -- 08's STEP_PLAN
    has its final caprieval at step index 9)."""
    caprieval_dir = os.path.join(RUN_DIR_BASE, POSITIVE_CONTROL_NAME, "run1", "9_caprieval")
    tsv_path = os.path.join(caprieval_dir, "capri_ss.tsv")
    if not os.path.exists(tsv_path):
        return None
    with open(tsv_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    top = next((r for r in rows if r["caprieval_rank"] == "1"), None)
    if top is None:
        return None
    resolved = os.path.normpath(os.path.join(caprieval_dir, top["model"]))
    if os.path.exists(resolved):
        return resolved
    if os.path.exists(resolved + ".gz"):
        return resolved + ".gz"
    return None


def locked_stats():
    """The positive control's own locked-schema stats, recomputed here from the
    same capri_ss.tsv 08 reads (mean/best/sd over the 10 best-scoring models,
    plus the best model's cluster and the number of models scored). Duplicated
    rather than imported because 08's module name isn't a valid identifier and
    importing it would fire its whole docking main(). Returns None if no models
    were scored -- callers just skip the JSON in that case."""
    caprieval_dir = os.path.join(RUN_DIR_BASE, POSITIVE_CONTROL_NAME, "run1", "9_caprieval")
    ss_path = os.path.join(caprieval_dir, "capri_ss.tsv")
    if not os.path.exists(ss_path):
        return None
    with open(ss_path) as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    if len(lines) < 2:
        return None
    header = lines[0].strip().split("\t")
    recs = [dict(zip(header, l.strip().split("\t"))) for l in lines[1:]]
    recs.sort(key=lambda r: float(r["score"]))
    best10 = [float(r["score"]) for r in recs[:10]]
    return {
        "mean_best_10": round(statistics.fmean(best10), 3),
        "best_haddock": round(best10[0], 3),
        "sd_best_10": round(statistics.stdev(best10), 3) if len(best10) > 1 else 0.0,
        "cluster_id": recs[0].get("cluster_id", "-"),
        "n_models": len(recs),
    }


def write_dockq_to_controls_csv(dockq):
    """Fill dockq_vs_9dwv on the positive control's row in 08's controls CSV,
    leaving every other row and column byte-identical. 08 preserves whatever it
    finds in this column on re-runs (it only ever writes ""), so this value
    survives a later 08 invocation."""
    if not os.path.exists(CONTROLS_CSV):
        print(f"NOTE: {os.path.basename(CONTROLS_CSV)} not found -- skipping CSV update. "
              f"Run 08 with RUN_CONTROLS = True first.")
        return
    with open(CONTROLS_CSV, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if "dockq_vs_9dwv" not in (fieldnames or []):
        print(f"NOTE: {os.path.basename(CONTROLS_CSV)} has no dockq_vs_9dwv column "
              f"(pre-lock schema?) -- skipping CSV update.")
        return
    target = next((r for r in rows if r["molecule"] == POSITIVE_CONTROL_NAME), None)
    if target is None:
        print(f"NOTE: no {POSITIVE_CONTROL_NAME} row in {os.path.basename(CONTROLS_CSV)} "
              f"-- skipping CSV update.")
        return
    previous = target["dockq_vs_9dwv"]
    target["dockq_vs_9dwv"] = f"{dockq:.3f}"
    with open(CONTROLS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    changed = f" (was {previous!r})" if previous else ""
    print(f"Wrote dockq_vs_9dwv={dockq:.3f} to {POSITIVE_CONTROL_NAME} in "
          f"{os.path.basename(CONTROLS_CSV)}{changed}")


def write_reference_json(result):
    """PROTOCOL_LOCK.md §5's golden reference: the numbers a second wrapper must
    reproduce (same cluster, DockQ within ~0.05) before its candidate runs are
    trusted. Commit this file."""
    stats = locked_stats()
    if stats is None:
        print("NOTE: no capri_ss.tsv for the positive control -- skipping "
              f"{os.path.basename(REFERENCE_JSON_PATH)}.")
        return
    payload = {
        "molecule": POSITIVE_CONTROL_NAME,
        "smiles": POSITIVE_CONTROL_SMILES,
        "protocol_version": PROTOCOL_VERSION,
        "reference_pdb_id": FPFT2216_PDB_ID,
        **stats,
        "dockq_vs_9dwv": round(result.dockq, 3),
        "irmsd": round(result.irmsd, 3),
        "fnat": round(result.fnat, 3),
        "lrmsd": round(result.lrmsd, 3),
        "ilrmsd": round(result.ilrmsd, 3),
        "rmsd": round(result.rmsd, 3),
    }
    with open(REFERENCE_JSON_PATH, "w") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")
    print(f"Wrote golden reference (PROTOCOL_LOCK §5): {os.path.basename(REFERENCE_JSON_PATH)}")


def strip_ligand_chain(model_path):
    """08's model has 3 chains (A=CRBN, C=ligand, B=PPIL4) -- drop chain C
    so this is a straight 2-chain protein-protein comparison against the
    reference, which never had the ligand in it either."""
    opener = gzip.open if model_path.endswith(".gz") else open
    with opener(model_path, "rt") as f:
        lines = f.readlines()
    out_path = os.path.join(SCRATCH_DIR, "positive_control_no_ligand.pdb")
    with open(out_path, "w") as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and line[21] == "C":
                continue
            f.write(line)
    return out_path


def main():
    download_cif()
    build_reference_pdb()
    os.makedirs(SCRATCH_DIR, exist_ok=True)

    model_path = find_top_model_path()
    if model_path is None:
        sys.exit(
            f"No rank-1 model found for {POSITIVE_CONTROL_NAME} under {RUN_DIR_BASE}. "
            f"Run the setup steps in this script's docstring first: "
            f"04 with MANUAL_CANDIDATE = ({POSITIVE_CONTROL_NAME!r}, {POSITIVE_CONTROL_SMILES!r}), "
            f"then 08 with CANDIDATE_NAME = {POSITIVE_CONTROL_NAME!r}."
        )
    no_ligand_path = strip_ligand_chain(model_path)

    from pathlib import Path
    from haddock.libs.libcapri import CAPRI
    from haddock.gear.yaml2cfg import read_from_yaml_config
    from haddock.modules.analysis.caprieval import DEFAULT_CONFIG
    capri_params = read_from_yaml_config(DEFAULT_CONFIG)

    print(f"Scoring {POSITIVE_CONTROL_NAME}'s rank-1 model (ligand chain stripped) "
          f"against the real FPFT-2216/9DWV reference structure...")
    capri = CAPRI(
        identificator=1,
        model=Path(no_ligand_path),
        path=Path(SCRATCH_DIR),
        reference=Path(REFERENCE_PDB_PATH),
        params=capri_params,
        ref_id=1,
        ff="aa",
    )
    result = capri.run()
    if result is None:
        sys.exit("Alignment against the FPFT-2216 reference failed.")

    print(f"\ndockq={result.dockq:.3f}  irmsd={result.irmsd:.2f}  fnat={result.fnat:.3f}  "
          f"lrmsd={result.lrmsd:.2f}  ilrmsd={result.ilrmsd:.2f}  rmsd={result.rmsd:.2f}")

    # Persist it -- until now this number only ever reached stdout, which is why
    # 08's dockq_vs_9dwv column stayed empty even after 09 had been run.
    print()
    write_dockq_to_controls_csv(result.dockq)
    write_reference_json(result)

    if result.dockq >= 0.23:
        print("\nPASS-ish: dockq >= 0.23 (CAPRI 'acceptable' threshold) -- the pipeline's own "
              "docking protocol can get at least roughly the right ternary architecture when "
              "given the real drug. This supports treating 08's candidate rankings as meaningful, "
              "but is still one data point, not proof every candidate's pose is correct.")
    else:
        print("\nFAIL: dockq < 0.23 (below CAPRI 'acceptable') -- the pipeline can't reproduce the "
              "known real complex even with the real drug. This is a problem with the methodology "
              "itself (restraints/sampling/force field), not with any particular candidate -- "
              "worth root-causing before trusting 08's candidate rankings.")

    total = time.time() - SCRIPT_START_TIME
    print(f"\nTotal script runtime: {total:.0f}s ({total / 60:.1f} min)")


if __name__ == "__main__":
    main()
