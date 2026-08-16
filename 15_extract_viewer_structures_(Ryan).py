"""
Extract a compact, browser-ready form of every ternary model on the project page
-- 20 candidates, 6 controls, and the experimental 9DWV reference -- into
CRBN_Project_site/structures.json, for the page's interactive cartoon viewer.

WHAT IS KEPT, AND WHY IT IS ENOUGH. Full models are ~6,000 atoms each; shipping
27 of them would be tens of megabytes. A cartoon representation only needs, per
residue: where the C-alpha is, which way the peptide plane faces (so the ribbon
can be given a width direction rather than being a bare tube), and what
secondary structure the residue is in (so helices and strands can be drawn wide
and loops thin). That is three numbers, three numbers and one character per
residue, plus the ligand's heavy atoms and bonds -- a few hundred KB for all 27
instead of tens of megabytes.

WHY THIS RUNS UNDER PYMOL. Secondary structure is not in the model files; it has
to be assigned from geometry. Rather than reimplement DSSP, this uses PyMOL's
own `dss`, which is the same assignment the static renders use -- so the
interactive cartoon and the rendered image agree about what is a helix.
Superposition uses PyMOL's `super` for the same reason.

EVERYTHING IS SUPERPOSED ON CRBN. CRBN is the fixed part of the complex, so
putting it in one place across all structures makes PPIL4's position genuinely
comparable between molecules rather than an artifact of each model's arbitrary
frame.

The reference is handled differently in one respect: 9DWV is an experimental
structure, so its ligand is the real FPFT-2216 read from the deposited mmCIF
rather than anything this pipeline produced. Its PPIL4 chain is also short --
only the RRM domain was resolved -- while the docked models carry a full
AlphaFold PPIL4, so the reference looks sparser. That is the experiment, not
the extraction.

Coordinates are rounded to 0.1 A, far below any resolution this supports.

Run under PyMOL, from the project directory:
    /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq "15_extract_viewer_structures_(Ryan).py"
"""
import csv
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
OUT_JSON = os.path.join(SITE_DIR, "structures.json")

BUY_LIST = os.path.join(SCRIPT_DIR, BUY_LIST_NAME)
CONTROLS_CSV = os.path.join(SCRIPT_DIR, "08_controls_results_(Ryan).csv")
REFERENCE_PDB = os.path.join(SCRIPT_DIR, "reference_structures", "FPFT-2216_9DWV_reference_(Ryan).pdb")
REFERENCE_CIF = os.path.join(SCRIPT_DIR, "reference_structures", "FPFT-2216_9DWV.cif")
REFERENCE_LIGAND_CODE = "A1BC8"          # FPFT-2216's PDB chemical component id

RUN_DIR_BASES = [
    os.path.join(os.path.expanduser("~"), "haddock_runs", "haddock3_ternary_with_ligand_run"),
    os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run"),
]

BOND_CUTOFF = 1.95      # longest plausible heavy-atom covalent bond


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


def reference_ligand_atoms():
    """The real FPFT-2216 from the deposited mmCIF -- not a docked pose."""
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
                atoms.append((rec["type_symbol"].upper(),
                              (float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"]))))
    return atoms


def bonds_for(atoms):
    out = []
    for i in range(len(atoms)):
        xi, yi, zi = atoms[i][1]
        for j in range(i + 1, len(atoms)):
            xj, yj, zj = atoms[j][1]
            if (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2 <= BOND_CUTOFF ** 2:
                out.append([i, j])
    return out


def chain_cartoon(cmd, obj, chain):
    """Per-residue (CA, carbonyl O, secondary structure) for one chain.

    The O is what gives the ribbon a width direction: CA alone defines a tube
    but not which way the peptide plane faces, so a CA-only cartoon twists
    arbitrarily. Residues missing either atom are dropped rather than guessed."""
    # NB: inside iterate_state's namespace `ss` is the residue's secondary
    # structure but `s` is PyMOL's settings object -- reading `s` raises
    # "unknown setting". Keep the local dict under a different name.
    ca, ox, sec = {}, {}, {}
    cmd.iterate_state(1, f"{obj} and chain {chain} and name CA and polymer",
                      "ca[int(resi)] = (x, y, z); sec[int(resi)] = ss",
                      space={"ca": ca, "sec": sec, "int": int})
    cmd.iterate_state(1, f"{obj} and chain {chain} and name O and polymer",
                      "ox[int(resi)] = (x, y, z)", space={"ox": ox, "int": int})
    out = []
    for resi in sorted(ca):
        if resi not in ox:
            continue
        out.append((resi, ca[resi], ox[resi], (sec.get(resi) or "L")[:1].upper()))
    return out


def flat(points, r=1):
    return [round(float(v), r) for p in points for v in p]


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

    out, reference_obj = {}, None
    for display, run_name, role in entries:
        path = top_model_path(run_name)
        if path is None:
            print(f"  {display}: no rank-1 model -- skipped")
            continue
        obj = "s" + str(len(out))
        cmd.load(path, obj)
        cmd.dss(obj)                                   # assign secondary structure
        if reference_obj is None:
            reference_obj = obj
        else:
            cmd.super(f"{obj} and chain A and name CA", f"{reference_obj} and chain A and name CA")

        lig = []
        cmd.iterate_state(1, f"{obj} and chain C and not elem H",
                          "lig.append((elem.upper(), (x, y, z)))", space={"lig": lig})
        a, b = chain_cartoon(cmd, obj, "A"), chain_cartoon(cmd, obj, "B")
        out[display] = {
            "role": role,
            "a": {"ca": flat([r[1] for r in a]), "o": flat([r[2] for r in a]),
                  "ss": "".join(r[3] for r in a)},
            "b": {"ca": flat([r[1] for r in b]), "o": flat([r[2] for r in b]),
                  "ss": "".join(r[3] for r in b)},
            "el": [x[0] for x in lig],
            "lig": flat([x[1] for x in lig]),
            "bonds": bonds_for(lig),
        }
        cmd.delete(obj) if False else None             # keep loaded for superposition frame

    # The experimental reference, superposed into the same CRBN frame.
    if os.path.exists(REFERENCE_PDB) and reference_obj:
        obj = "sref"
        cmd.load(REFERENCE_PDB, obj)
        cmd.dss(obj)
        cmd.super(f"{obj} and chain A and name CA", f"{reference_obj} and chain A and name CA")
        # The mmCIF ligand is in the deposited frame; move it by the same
        # transform PyMOL just applied to the reference's protein chains.
        cmd.create("sref_lig", "none")
        lig_atoms = reference_ligand_atoms()
        if lig_atoms:
            pdb = "".join(
                f"HETATM{i + 1:5d} {('X' + str(i))[:4]:<4s} LIG X 900    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el[:2]:>2s}\n"
                for i, (el, (x, y, z)) in enumerate(lig_atoms))
            tmp = os.path.join(SITE_DIR, "_ref_lig.pdb")
            os.makedirs(SITE_DIR, exist_ok=True)
            with open(tmp, "w") as f:
                f.write(pdb + "END\n")
            cmd.load(tmp, "sref_lig")
            matrix = cmd.get_object_matrix(obj)
            cmd.transform_object("sref_lig", matrix)
            os.remove(tmp)
        moved = []
        cmd.iterate_state(1, "sref_lig", "moved.append((elem.upper(), (x, y, z)))",
                          space={"moved": moved})
        a, b = chain_cartoon(cmd, obj, "A"), chain_cartoon(cmd, obj, "B")
        out["9DWV_reference"] = {
            "role": "reference",
            "a": {"ca": flat([r[1] for r in a]), "o": flat([r[2] for r in a]),
                  "ss": "".join(r[3] for r in a)},
            "b": {"ca": flat([r[1] for r in b]), "o": flat([r[2] for r in b]),
                  "ss": "".join(r[3] for r in b)},
            "el": [x[0] for x in moved],
            "lig": flat([x[1] for x in moved]),
            "bonds": bonds_for(moved),
        }
        print(f"  9DWV reference: CRBN {len(a)} res, PPIL4 {len(b)} res, "
              f"ligand {len(moved)} heavy atoms (real FPFT-2216 from the mmCIF)")

    os.makedirs(SITE_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, separators=(",", ":"))

    residues = sum(len(v["a"]["ss"]) + len(v["b"]["ss"]) for v in out.values())
    helix = sum(v["a"]["ss"].count("H") + v["b"]["ss"].count("H") for v in out.values())
    sheet = sum(v["a"]["ss"].count("S") + v["b"]["ss"].count("S") for v in out.values())
    size = os.path.getsize(OUT_JSON) / 1024
    print(f"\nWrote {os.path.relpath(OUT_JSON, SCRIPT_DIR)}: {len(out)} structures, "
          f"{residues:,} residues ({helix:,} helix, {sheet:,} strand), {size:.0f} KB")
    print("All superposed on CRBN; secondary structure assigned by PyMOL's dss, "
          "so the interactive cartoon and the static renders agree.")


main()
