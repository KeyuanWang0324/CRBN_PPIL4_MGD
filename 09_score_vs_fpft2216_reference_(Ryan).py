"""
Positive control: does this project's docking METHODOLOGY (04's Vina pose
-> 08's ligand-inclusive 3-body HADDOCK3 protocol) actually reproduce a
KNOWN real CRBN-PPIL4 ternary complex, when given the real drug that
complex was solved with?

This DOES score every molecule 08 docked -- all candidates plus all
controls -- but read the next two paragraphs before using those numbers,
because what a per-candidate dockq means changed when the protocol was
locked, and it still does NOT mean what it looks like it means.

THE ABSOLUTE VALUES ARE INFLATED AND ARE NOT VALIDATION. 05/06/08's PPIL4
restraint ([249, 250, 273, 275, 276, 277, 278, 279]) is ITSELF derived
from 9DWV's own real contact residues. Scoring a molecule's dockq against
9DWV after handing its restraints the answer key isn't independent
validation -- PPIL4 lands close to the real interface almost by
construction (see 00_validate_docking_interface_05_(Ryan).py: 100% "found
the spot" on PPIL4 across the whole 05 population, mean 6.3/8 residues).
A high dockq here has never been evidence that a molecule is a good
binder; on its own it mostly confirms the restraint was applied. Ranking
still comes from 08's own HADDOCK score, NOT from this column.

WHAT THE LOCK CHANGED, AND WHY IT STILL ISN'T ENOUGH. When the earlier
version of this script was removed, 08 derived the CRBN-side active
residues from EACH candidate's own Vina pose, so every molecule was docked
under a DIFFERENT restraint set and the dockq spread was mostly restraint
noise. Under ppil4_lock_v1 every molecule and control is docked against
ONE committed restraint file (ambig_FIXED.tbl, PROTOCOL_LOCK.md §2), so
the "answer key" is now identical for everyone and the variation between
molecules is no longer explained by the restraints. That argument is what
motivated running the full sweep.

The sweep then falsified the useful half of it. Measured results (24
molecules, ppil4_lock_v1):
  - FPFT-2216, the real drug, tops the list at dockq 0.376 -- the
    methodology check in main() passes, and passes cleanly.
  - But negctrl_no_ppil4_arm (0.294) and negctrl_no_crbn_ppil4_arm (0.290)
    rank 2nd and 3rd, above EVERY candidate. The second of those is the
    fragment with both functional handles stripped and the worst HADDOCK
    score of anything docked (-92.19).
  - Z6466608628, the strongest positive control by HADDOCK score and the
    buy list's rule-C bar, scores 0.061.
  - Spearman vs 08's HADDOCK score: +0.055. Vs ligand heavy-atom count:
    +0.066. So it is neither a restatement of the score nor a size
    artifact -- it is close to independent of both.
  - The distribution is bimodal (irmsd ~2-3 A vs ~5-9 A), i.e. it mostly
    records whether THIS ONE rank-1 model out of ~160 happened to land in
    the native-like basin. That is a single-model sampling outcome, not a
    smooth property of the molecule.

So: a molecule cannot be preferred or rejected on this column. It does not
separate glues from degron-dead fragments -- the negatives beat the
candidates. Treat it as a per-run geometry diagnostic (did the top model
land on the real interface?), useful for spotting which models are worth
looking at in PyMOL, and nothing more. Ranking stays with 08's HADDOCK
score; the buy list's structure axis is built from BSA / cluster
convergence / AIR / sd instead, none of which use this column.

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

No re-docking is done here -- this only re-evaluates the rank-1 model 08
already produced for each molecule, scoring it with HADDOCK3's own CAPRI
class (haddock.libs.libcapri.CAPRI) against the real FPFT-2216/9DWV
reference. 08's model has 3 chains (A=CRBN, C=ligand, B=PPIL4); the
ligand chain is stripped before scoring, since the reference (and
dockq/irmsd/fnat/lrmsd) are about the CRBN-PPIL4 protein-protein
interface specifically -- matching how the reference itself was built
(DDB1/ligand/Zn dropped, see build_reference_pdb()). Each molecule takes
about a second, so the whole sweep is well under a minute.

Outputs (this script used to only print, which is why 08's dockq_vs_9dwv
column stayed empty even after 09 had been run):
  - fills dockq_vs_9dwv for EVERY molecule present in 08's two result
    CSVs -- 08_final_ternary_with_ligand_results_(Ryan).csv (candidates)
    and 08_controls_results_(Ryan).csv (controls). No other column is
    touched: PROTOCOL_LOCK.md §6 fixes the schema, so the caveats above
    live in this docstring rather than in an extra CSV column.
  - writes REFERENCE_FPFT.json, PROTOCOL_LOCK.md §5's golden reference.

08's runs are split across two roots -- ~/haddock_runs (where the
candidates were run, off the iCloud-synced Desktop) and the in-project
docking_tmp/ (where the controls were) -- so RUN_DIR_BASES searches both.

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

# 08's runs live in two places: the candidates were moved to ~/haddock_runs to
# get them off the iCloud-synced Desktop (which corrupts in-flight run dirs),
# while the controls were run before that move and are still under docking_tmp/.
# Searched in order; the first base holding a completed 9_caprieval wins.
RUN_DIR_BASES = [
    os.path.join(os.path.expanduser("~"), "haddock_runs", "haddock3_ternary_with_ligand_run"),
    os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run"),
]
SCRATCH_DIR = os.path.join(SCRIPT_DIR, "docking_tmp", "09_capri_vs_fpft2216_scratch")

REFERENCE_DIR = os.path.join(SCRIPT_DIR, "reference_structures")
CIF_CACHE_PATH = os.path.join(REFERENCE_DIR, "FPFT-2216_9DWV.cif")
REFERENCE_PDB_PATH = os.path.join(REFERENCE_DIR, "FPFT-2216_9DWV_reference_(Ryan).pdb")

# Where this script's results go -- both of 08's result CSVs. 08 writes every
# row with dockq_vs_9dwv deliberately blank (its own caprieval dockq is
# self-referential: reference = the run's own best model), and its locked_stats()
# docstring says "09 fills the real value against 9DWV". This is that fill.
# Only this column is touched; PROTOCOL_LOCK.md §6 fixes the rest of the schema.
FINALISTS_CSV = os.path.join(SCRIPT_DIR, "08_final_ternary_with_ligand_results_(Ryan).csv")
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


def find_caprieval_dir(molecule):
    """08's final caprieval directory for one molecule, looked up across both
    run roots (see RUN_DIR_BASES). 08's STEP_PLAN puts its final caprieval at
    step index 9, same lookup 07/08 use. Returns None if no root has a
    completed run for this molecule."""
    for base in RUN_DIR_BASES:
        caprieval_dir = os.path.join(base, molecule, "run1", "9_caprieval")
        if os.path.exists(os.path.join(caprieval_dir, "capri_ss.tsv")):
            return caprieval_dir
    return None


def find_top_model_path(molecule):
    """The rank-1 model from 08's ligand-inclusive run for one molecule."""
    caprieval_dir = find_caprieval_dir(molecule)
    if caprieval_dir is None:
        return None
    tsv_path = os.path.join(caprieval_dir, "capri_ss.tsv")
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


def locked_stats(molecule):
    """One molecule's locked-schema stats, recomputed here from the same
    capri_ss.tsv 08 reads (mean/best/sd over the 10 best-scoring models, plus
    the best model's cluster and the number of models scored). Duplicated
    rather than imported because 08's module name isn't a valid identifier and
    importing it would fire its whole docking main(). Returns None if no models
    were scored -- callers just skip the JSON in that case."""
    caprieval_dir = find_caprieval_dir(molecule)
    if caprieval_dir is None:
        return None
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


def read_molecules(csv_path):
    """The molecule names 08 wrote into one of its result CSVs, in file order.
    Rows 08 marked incomplete ("-") are skipped -- there's no model to score."""
    if not os.path.exists(csv_path):
        return []
    with open(csv_path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [r["molecule"] for r in rows if r.get("mean_best_10") not in (None, "", "-")]


def write_dockq_to_csv(csv_path, dockq_by_molecule):
    """Fill dockq_vs_9dwv for every scored molecule in one of 08's result CSVs,
    leaving every other row and column byte-identical -- PROTOCOL_LOCK.md §6
    fixes this schema, so no column is added or reordered. 08 preserves whatever
    it finds in this column on re-runs (it only ever writes ""), so these values
    survive a later 08 invocation."""
    name = os.path.basename(csv_path)
    if not os.path.exists(csv_path):
        print(f"NOTE: {name} not found -- skipping.")
        return 0
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if "dockq_vs_9dwv" not in (fieldnames or []):
        print(f"NOTE: {name} has no dockq_vs_9dwv column (pre-lock schema?) -- skipping.")
        return 0
    written = 0
    for row in rows:
        dockq = dockq_by_molecule.get(row["molecule"])
        if dockq is None:
            continue
        row["dockq_vs_9dwv"] = f"{dockq:.3f}"
        written += 1
    if not written:
        print(f"NOTE: nothing to write into {name}.")
        return 0
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote dockq_vs_9dwv for {written}/{len(rows)} molecules -> {name}")
    return written


def write_reference_json(result):
    """PROTOCOL_LOCK.md §5's golden reference: the numbers a second wrapper must
    reproduce (same cluster, DockQ within ~0.05) before its candidate runs are
    trusted. Commit this file."""
    stats = locked_stats(POSITIVE_CONTROL_NAME)
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


def strip_ligand_chain(model_path, molecule):
    """08's model has 3 chains (A=CRBN, C=ligand, B=PPIL4) -- drop chain C
    so this is a straight 2-chain protein-protein comparison against the
    reference, which never had the ligand in it either. Written per molecule
    so a sweep doesn't have every model overwriting one scratch file."""
    opener = gzip.open if model_path.endswith(".gz") else open
    with opener(model_path, "rt") as f:
        lines = f.readlines()
    out_path = os.path.join(SCRATCH_DIR, f"{molecule}_no_ligand.pdb")
    with open(out_path, "w") as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and line[21] == "C":
                continue
            f.write(line)
    return out_path


def score_molecule(molecule, capri_params, identificator):
    """CAPRI metrics for one molecule's rank-1 model against the real 9DWV
    reference. Returns the result object, or None if the run is missing or the
    alignment fails -- one bad molecule must not abort the whole sweep."""
    from pathlib import Path
    from haddock.libs.libcapri import CAPRI

    model_path = find_top_model_path(molecule)
    if model_path is None:
        return None
    no_ligand_path = strip_ligand_chain(model_path, molecule)
    capri = CAPRI(
        identificator=identificator,
        model=Path(no_ligand_path),
        path=Path(SCRATCH_DIR),
        reference=Path(REFERENCE_PDB_PATH),
        params=capri_params,
        ref_id=1,
        ff="aa",
    )
    try:
        return capri.run()
    except Exception as exc:                       # noqa: BLE001 -- keep the sweep going
        print(f"  {molecule}: CAPRI failed ({type(exc).__name__}: {exc})")
        return None


def main():
    download_cif()
    build_reference_pdb()
    os.makedirs(SCRATCH_DIR, exist_ok=True)

    from haddock.gear.yaml2cfg import read_from_yaml_config
    from haddock.modules.analysis.caprieval import DEFAULT_CONFIG
    capri_params = read_from_yaml_config(DEFAULT_CONFIG)

    targets = [(FINALISTS_CSV, read_molecules(FINALISTS_CSV)),
               (CONTROLS_CSV, read_molecules(CONTROLS_CSV))]
    total_molecules = sum(len(mols) for _, mols in targets)
    if not total_molecules:
        sys.exit(
            f"No molecules found in {os.path.basename(FINALISTS_CSV)} or "
            f"{os.path.basename(CONTROLS_CSV)}. Run 08 first."
        )

    print(f"Scoring {total_molecules} molecule(s) from 08 (rank-1 model, ligand chain "
          f"stripped) against the real FPFT-2216/9DWV reference structure...\n")

    results = {}
    identificator = 0
    for csv_path, molecules in targets:
        for molecule in molecules:
            identificator += 1
            result = score_molecule(molecule, capri_params, identificator)
            if result is None:
                print(f"  {molecule:<32} no rank-1 model found -- skipped")
                continue
            results[molecule] = result
            print(f"  {molecule:<32} dockq={result.dockq:.3f}  irmsd={result.irmsd:5.2f}  "
                  f"fnat={result.fnat:.3f}  lrmsd={result.lrmsd:6.2f}")

    if not results:
        sys.exit(
            "Nothing could be scored. 08's run directories were not found under any of: "
            + ", ".join(RUN_DIR_BASES)
        )

    # Persist -- these numbers used to only reach stdout, which is why 08's
    # dockq_vs_9dwv column stayed empty even after 09 had been run.
    print()
    dockq_by_molecule = {m: r.dockq for m, r in results.items()}
    for csv_path, _ in targets:
        write_dockq_to_csv(csv_path, dockq_by_molecule)

    reference = results.get(POSITIVE_CONTROL_NAME)
    if reference is None:
        print(f"\nNOTE: {POSITIVE_CONTROL_NAME} was not scored -- skipping "
              f"{os.path.basename(REFERENCE_JSON_PATH)} and the acceptance check below. "
              f"Run 04 + 08 for it (see this script's docstring).")
        total = time.time() - SCRIPT_START_TIME
        print(f"\nTotal script runtime: {total:.0f}s ({total / 60:.1f} min)")
        return
    write_reference_json(reference)

    # The methodology check is about the REAL drug only: it is the one molecule
    # whose correct answer is independently known. Everything else in the sweep
    # is a relative descriptor measured against this value (see module docstring).
    if reference.dockq >= 0.23:
        print("\nPASS-ish: dockq >= 0.23 (CAPRI 'acceptable' threshold) -- the pipeline's own "
              "docking protocol can get at least roughly the right ternary architecture when "
              "given the real drug. This supports treating 08's candidate rankings as meaningful, "
              "but is still one data point, not proof every candidate's pose is correct.")
    else:
        print("\nFAIL: dockq < 0.23 (below CAPRI 'acceptable') -- the pipeline can't reproduce the "
              "known real complex even with the real drug. This is a problem with the methodology "
              "itself (restraints/sampling/force field), not with any particular candidate -- "
              "worth root-causing before trusting 08's candidate rankings.")

    # Calibration check: if the negative controls aren't at the BOTTOM of this
    # column, the column doesn't separate glues from degron-dead fragments and
    # must not be used to prefer or reject a molecule. As of ppil4_lock_v1 they
    # are not -- two of the three outrank every candidate (see module docstring).
    negatives = {m: r.dockq for m, r in results.items() if m.startswith("negctrl_")}
    candidates = {m: r.dockq for m, r in results.items() if m.startswith("cand_")}
    if negatives and candidates:
        best_negative = max(negatives, key=negatives.get)
        beaten = [m for m, d in candidates.items() if d < negatives[best_negative]]
        print(f"\nCALIBRATION: the best negative control ({best_negative}, "
              f"dockq={negatives[best_negative]:.3f}) outranks {len(beaten)}/{len(candidates)} "
              f"candidates on this column.")
        if len(beaten) > len(candidates) / 2:
            print("=> dockq_vs_9dwv does NOT separate glues from degron-dead fragments. Use it "
                  "only as a per-run geometry diagnostic (did the top model land on the real "
                  "interface?). Do NOT rank, prefer, or reject molecules on it -- ranking stays "
                  "with 08's HADDOCK score.")

    total = time.time() - SCRIPT_START_TIME
    print(f"\nTotal script runtime: {total:.0f}s ({total / 60:.1f} min)")


if __name__ == "__main__":
    main()
