"""
Write a trimmed, viewer-ready PDB for every ternary model on the project page --
20 candidates, 6 controls, and the experimental 9DWV reference -- into
CRBN_Project_site/pdb/, for the page's interactive 3Dmol viewer.

WHY PDB RATHER THAN A CUSTOM FORMAT. The page previously shipped hand-rolled
coordinate arrays and drew the cartoon itself. Reproducing what PyMOL does --
ribbon framing, helix smoothing, strand arrows, sheet twist -- by hand gets you
something that is recognisably not a cartoon. Handing a real PDB to a real
molecular viewer gets the real thing. So this writes standard PDB and the page
renders it with 3Dmol.js.

WHAT IS TRIMMED, AND WHY IT IS SAFE. Cartoons are built from the backbone, so
side chains and hydrogens are dropped and only N, CA, C and O are kept for the
protein. The ligand is kept in full, because it is drawn as ball-and-stick and
every atom matters there. That takes a model from ~6,000 atoms to ~2,400 --
about 200 KB per structure, fetched only when a record is opened.

SECONDARY STRUCTURE IS WRITTEN EXPLICITLY. PyMOL assigns it with `dss` but does
not write HELIX/SHEET records when saving PDB, and a viewer left to guess will
disagree with the static renders about what is a helix. So the assignment is
read back out of PyMOL and emitted as proper HELIX/SHEET records, which 3Dmol
honours -- making the interactive cartoon and the rendered images agree.

EVERYTHING IS SUPERPOSED ON CRBN, so PPIL4's position is comparable between
structures rather than an artifact of each model's arbitrary frame.

The reference is the deposited 9DWV entry, with the real FPFT-2216 read from the
deposited chemical component rather than from any docking run. Its PPIL4 chain
is short because only the RRM domain was resolved experimentally.

Run under PyMOL, from the project directory:
    /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq "15_extract_viewer_structures_(Ryan).py"
"""
import csv
import hashlib
import json
import os
import sys

BUY_LIST_NAME = "final_buy_list_lock_v1_top20_(Ryan).csv"


def _project_dir():
    """PyMOL's `run` rebinds __file__ to its own package directory."""
    for candidate in (os.path.dirname(os.path.abspath(__file__)),
                      os.getcwd(),
                      os.environ.get("ISEF_DIR", "")):
        if candidate and os.path.exists(os.path.join(candidate, BUY_LIST_NAME)):
            return candidate
    sys.exit(f"Cannot locate {BUY_LIST_NAME}. Run from the project directory or set $ISEF_DIR.")


SCRIPT_DIR = _project_dir()
SITE_DIR = os.path.join(SCRIPT_DIR, "CRBN_Project_site")
PDB_DIR = os.path.join(SITE_DIR, "pdb")
INDEX_JSON = os.path.join(SITE_DIR, "structures.json")

BUY_LIST = os.path.join(SCRIPT_DIR, BUY_LIST_NAME)
CONTROLS_CSV = os.path.join(SCRIPT_DIR, "08_controls_results_(Ryan).csv")
REFERENCE_PDB = os.path.join(SCRIPT_DIR, "reference_structures", "FPFT-2216_9DWV_reference_(Ryan).pdb")
REFERENCE_CIF = os.path.join(SCRIPT_DIR, "reference_structures", "FPFT-2216_9DWV.cif")
REFERENCE_LIGAND_CODE = "A1BC8"

RUN_DIR_BASES = [
    os.path.join(os.path.expanduser("~"), "haddock_runs", "haddock3_ternary_with_ligand_run"),
    os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run"),
]


def top_model_path(run_name):
    for base in RUN_DIR_BASES:
        capri = os.path.join(base, run_name, "run1", "9_caprieval")
        tsv = os.path.join(capri, "capri_ss.tsv")
        if not os.path.exists(tsv):
            continue
        with open(tsv, newline="") as f:
            top = next((r for r in csv.DictReader(f, delimiter="\t")
                        if r["caprieval_rank"] == "1"), None)
        if top is None:
            continue
        resolved = os.path.normpath(os.path.join(capri, top["model"]))
        for path in (resolved, resolved + ".gz"):
            if os.path.exists(path):
                return path
    return None


def _col(buf, start, text):
    """Place `text` at 1-based PDB column `start`."""
    i = start - 1
    buf[i:i + len(text)] = list(text)


def helix_record(serial, chain, name1, res1, name2, res2):
    buf = [" "] * 80
    _col(buf, 1, "HELIX ")
    _col(buf, 8, f"{serial:>3d}")
    _col(buf, 12, f"{('H' + str(serial))[:3]:>3s}")
    _col(buf, 16, f"{name1[:3]:>3s}")
    _col(buf, 20, chain)
    _col(buf, 22, f"{res1:>4d}")
    _col(buf, 28, f"{name2[:3]:>3s}")
    _col(buf, 32, chain)
    _col(buf, 34, f"{res2:>4d}")
    _col(buf, 39, f"{1:>2d}")                    # helixClass 1 = right-handed alpha
    _col(buf, 72, f"{res2 - res1 + 1:>5d}")
    return "".join(buf).rstrip()


def sheet_record(strand, chain, name1, res1, name2, res2, total):
    buf = [" "] * 80
    _col(buf, 1, "SHEET ")
    _col(buf, 8, f"{strand:>3d}")
    _col(buf, 12, f"{('S' + str(strand))[:3]:>3s}")
    _col(buf, 15, f"{total:>2d}")
    _col(buf, 18, f"{name1[:3]:>3s}")
    _col(buf, 22, chain)
    _col(buf, 23, f"{res1:>4d}")
    _col(buf, 29, f"{name2[:3]:>3s}")
    _col(buf, 33, chain)
    _col(buf, 34, f"{res2:>4d}")
    _col(buf, 39, f"{0:>2d}")                    # sense 0 = first strand / unknown
    return "".join(buf).rstrip()


def ss_records(cmd, obj):
    """HELIX/SHEET records for PyMOL's own dss assignment.

    Inside iterate's namespace `ss` is the secondary structure but `s` is the
    settings object -- reading `s` raises "unknown setting", hence the local
    dict being named something else."""
    runs, rows = [], []
    cmd.iterate(f"{obj} and name CA and polymer",
                "rows.append((chain, int(resi), resn, ss))",
                space={"rows": rows, "int": int})
    if not rows:
        return []
    current = None
    for chain, resi, resn, sec in rows:
        kind = (sec or "L")[:1].upper()
        if current and current[0] == kind and current[1] == chain and resi == current[4] + 1:
            current[4], current[5] = resi, resn
            continue
        if current and current[0] in ("H", "S"):
            runs.append(current)
        current = [kind, chain, resi, resn, resi, resn]
    if current and current[0] in ("H", "S"):
        runs.append(current)

    helices = [r for r in runs if r[0] == "H" and r[4] > r[2]]
    strands = [r for r in runs if r[0] == "S" and r[4] > r[2]]
    out = []
    for i, r in enumerate(helices, 1):
        out.append(helix_record(i, r[1], r[3], r[2], r[5], r[4]))
    for i, r in enumerate(strands, 1):
        out.append(sheet_record(i, r[1], r[3], r[2], r[5], r[4], len(strands)))
    return out


def reference_ligand_atoms():
    if not os.path.exists(REFERENCE_CIF):
        return []
    atoms, header, in_loop = [], [], False
    with open(REFERENCE_CIF) as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                header.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop:
                if not s or s.startswith(("#", "loop_", "_")):
                    break
                parts = s.split()
                if len(parts) < len(header):
                    continue
                rec = dict(zip(header, parts))
                if rec.get("label_comp_id") != REFERENCE_LIGAND_CODE:
                    continue
                if rec.get("type_symbol", "C").upper() == "H":
                    continue
                atoms.append((rec.get("label_atom_id", "C")[:4], rec["type_symbol"].upper(),
                              (float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"]))))
    return atoms


def write_structure(cmd, obj, out_path):
    """Trim to backbone + ligand and save, with SS records at the top."""
    cmd.remove(f"{obj} and hydro")
    cmd.remove(f"{obj} and polymer and not name N+CA+C+O")
    records = ss_records(cmd, obj)
    body = cmd.get_pdbstr(obj)
    with open(out_path, "w") as f:
        if records:
            f.write("\n".join(records) + "\n")
        f.write(body)
    with open(out_path, "rb") as f:
        version = hashlib.sha256(f.read()).hexdigest()[:10]
    return os.path.getsize(out_path), len(records), version


def main():
    from pymol import cmd

    entries = []
    with open(BUY_LIST, newline="") as f:
        for r in csv.DictReader(f):
            if r["role"] in ("candidate", "crosscheck"):
                entries.append((r["name"], r["run_name"], r["role"]))
    with open(CONTROLS_CSV, newline="") as f:
        for r in csv.DictReader(f):
            entries.append((r["molecule"], r["molecule"], "control"))

    os.makedirs(PDB_DIR, exist_ok=True)
    index, reference_obj, total = {}, None, 0

    for display, run_name, role in entries:
        path = top_model_path(run_name)
        if path is None:
            print(f"  {display}: no rank-1 model -- skipped")
            continue
        obj = "s%d" % len(index)
        cmd.load(path, obj)
        cmd.dss(obj)
        if reference_obj is None:
            reference_obj = obj
        else:
            cmd.super(f"{obj} and chain A and name CA", f"{reference_obj} and chain A and name CA")
        size, nrec, ver = write_structure(cmd, obj, os.path.join(PDB_DIR, f"{display}.pdb"))
        index[display] = {"role": role, "pdb": f"pdb/{display}.pdb?v={ver}"}
        total += size
        print(f"  {display:<32} {size // 1024:>4} KB, {nrec} SS records")

    if os.path.exists(REFERENCE_PDB) and reference_obj:
        obj = "sref"
        cmd.load(REFERENCE_PDB, obj)
        cmd.dss(obj)
        cmd.super(f"{obj} and chain A and name CA", f"{reference_obj} and chain A and name CA")
        # The deposited ligand is in the original frame; move it by the same
        # transform PyMOL just applied to the reference's protein chains.
        atoms = reference_ligand_atoms()
        if atoms:
            tmp = os.path.join(PDB_DIR, "_lig.pdb")
            with open(tmp, "w") as f:
                for i, (name, el, (x, y, z)) in enumerate(atoms, 1):
                    f.write(f"HETATM{i:5d} {name:<4s} LIG C 900    "
                            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el[:2]:>2s}\n")
                f.write("END\n")
            cmd.load(tmp, "sreflig")
            cmd.transform_object("sreflig", cmd.get_object_matrix(obj))
            cmd.create(obj, f"{obj} or sreflig", zoom=0)
            cmd.delete("sreflig")
            os.remove(tmp)
        size, nrec, ver = write_structure(cmd, obj, os.path.join(PDB_DIR, "9DWV_reference.pdb"))
        index["9DWV_reference"] = {"role": "reference", "pdb": f"pdb/9DWV_reference.pdb?v={ver}"}
        total += size
        print(f"  {'9DWV_reference':<32} {size // 1024:>4} KB, {nrec} SS records, "
              f"{len(atoms)} ligand atoms (real FPFT-2216)")

    with open(INDEX_JSON, "w") as f:
        json.dump(index, f, separators=(",", ":"), indent=0)
    print(f"\nWrote {len(index)} PDBs into {os.path.relpath(PDB_DIR, SCRIPT_DIR)}/ "
          f"({total / 1024 / 1024:.1f} MB total, fetched one at a time on demand)")
    print("Superposed on CRBN; HELIX/SHEET written from PyMOL's dss so the interactive "
          "cartoon matches the static renders.")


main()
