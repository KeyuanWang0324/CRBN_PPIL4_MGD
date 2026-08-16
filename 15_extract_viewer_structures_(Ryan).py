"""
Extract a compact, browser-ready form of every ternary model on the project
page -- 20 candidates, 6 controls, and the experimental 9DWV reference -- into
CRBN_Project_site/structures.json, for the page's interactive 3D viewer.

WHAT IS KEPT, AND WHY IT IS ENOUGH. Full models are ~6,000 atoms each; shipping
26 of them would be tens of megabytes. What the viewer actually draws is a
backbone trace plus the ligand, so only the C-alpha of each residue is kept for
the two protein chains, along with every heavy atom of the ligand and the bonds
between them. That is ~630 points per structure instead of ~6,000, and the whole
set lands in a few hundred KB.

EVERYTHING IS SUPERPOSED ON CRBN, for the same reason the static renders are:
CRBN is the fixed part of the complex, so putting it in one place across all
structures makes PPIL4's position genuinely comparable between molecules rather
than an artifact of each model's arbitrary frame. Superposition is by Kabsch on
the C-alpha atoms of the residue numbers the two chains have in common.

The reference is handled differently in one respect: the 9DWV entry is an
experimental structure, so its ligand is the real FPFT-2216 taken from the
deposited mmCIF rather than anything this pipeline produced. Its PPIL4 chain is
also short (74 residues) because only the RRM domain was resolved -- the docked
models carry a full AlphaFold PPIL4, so the reference will look sparser. That is
a property of the experiment, not of the extraction.

Coordinates are rounded to 0.1 A. That is far below the resolution any of this
supports and roughly halves the payload.

Run with the SYSTEM python (needs numpy), after 08/09:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "15_extract_viewer_structures_(Ryan).py"
"""
import csv
import gzip
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SITE_DIR = os.path.join(SCRIPT_DIR, "CRBN_Project_site")
OUT_JSON = os.path.join(SITE_DIR, "structures.json")

BUY_LIST = os.path.join(SCRIPT_DIR, "final_buy_list_lock_v1_top20_(Ryan).csv")
CONTROLS_CSV = os.path.join(SCRIPT_DIR, "08_controls_results_(Ryan).csv")
REFERENCE_PDB = os.path.join(SCRIPT_DIR, "reference_structures", "FPFT-2216_9DWV_reference_(Ryan).pdb")
REFERENCE_CIF = os.path.join(SCRIPT_DIR, "reference_structures", "FPFT-2216_9DWV.cif")
REFERENCE_LIGAND_CODE = "A1BC8"          # FPFT-2216's PDB chemical component id

RUN_DIR_BASES = [
    os.path.join(os.path.expanduser("~"), "haddock_runs", "haddock3_ternary_with_ligand_run"),
    os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run"),
]

# Longest plausible heavy-atom covalent bond; above this two atoms are not bonded.
BOND_CUTOFF = 1.95


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


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
        if os.path.exists(resolved):
            return resolved
        if os.path.exists(resolved + ".gz"):
            return resolved + ".gz"
    return None


def parse_model(path, ligand_chain="C"):
    """(chain A CA by resnum, chain B CA list, ligand heavy atoms) from a PDB."""
    ca_a, ca_b, lig = {}, [], []
    with open_maybe_gz(path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21]
            element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if chain == ligand_chain and element != "H":
                lig.append((element, xyz))
            elif line[12:16].strip() == "CA" and line.startswith("ATOM"):
                if chain == "A":
                    ca_a[int(line[22:26])] = xyz
                elif chain == "B":
                    ca_b.append(xyz)
    return ca_a, ca_b, lig


def parse_reference_ligand():
    """The real FPFT-2216 from the deposited mmCIF -- not a docked pose."""
    if not os.path.exists(REFERENCE_CIF):
        return []
    atoms = []
    with open(REFERENCE_CIF) as f:
        header, in_loop = [], False
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                header.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop:
                if s.startswith(("#", "loop_", "_")) or not s:
                    break
                parts = s.split()
                if len(parts) < len(header):
                    continue
                rec = dict(zip(header, parts))
                if rec.get("label_comp_id") != REFERENCE_LIGAND_CODE:
                    continue
                element = rec.get("type_symbol", "C").upper()
                if element == "H":
                    continue
                atoms.append((element, (float(rec["Cartn_x"]), float(rec["Cartn_y"]),
                                        float(rec["Cartn_z"]))))
    return atoms


def kabsch(mobile, target):
    """Rotation+translation putting `mobile` onto `target` (both N x 3)."""
    import numpy as np
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, tc - rot @ mc


def bonds_for(atoms):
    import numpy as np
    if not atoms:
        return []
    xyz = np.array([a[1] for a in atoms])
    d2 = ((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1)
    out = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if d2[i, j] <= BOND_CUTOFF ** 2:
                out.append([i, j])
    return out


def main():
    try:
        import numpy as np
    except ImportError:
        sys.exit("numpy required -- use the system python.")

    entries = []
    for r in csv.DictReader(open(BUY_LIST, newline="")):
        if r["role"] in ("candidate", "crosscheck"):
            entries.append((r["name"], r["run_name"], r["role"]))
    for r in csv.DictReader(open(CONTROLS_CSV, newline="")):
        entries.append((r["molecule"], r["molecule"], "control"))

    out, reference_ca = {}, None
    for display, run_name, role in entries:
        path = top_model_path(run_name)
        if path is None:
            print(f"  {display}: no rank-1 model -- skipped")
            continue
        ca_a, ca_b, lig = parse_model(path)
        if not ca_a:
            print(f"  {display}: no chain A CA -- skipped")
            continue

        if reference_ca is None:
            reference_ca = ca_a
            rot, trans = np.eye(3), np.zeros(3)
        else:
            common = sorted(set(ca_a) & set(reference_ca))
            if len(common) < 3:
                print(f"  {display}: only {len(common)} shared CRBN residues -- left unaligned")
                rot, trans = np.eye(3), np.zeros(3)
            else:
                rot, trans = kabsch(np.array([ca_a[i] for i in common]),
                                    np.array([reference_ca[i] for i in common]))

        def place(points):
            if not points:
                return []
            arr = (np.array(points) @ rot.T) + trans
            return [round(float(v), 1) for v in arr.reshape(-1)]

        out[display] = {
            "role": role,
            "a": place([ca_a[k] for k in sorted(ca_a)]),
            "b": place(ca_b),
            "el": [a[0] for a in lig],
            "lig": place([a[1] for a in lig]),
            "bonds": bonds_for(lig),
        }

    # The experimental reference, superposed onto the same CRBN frame.
    if os.path.exists(REFERENCE_PDB):
        ca_a, ca_b, _ = parse_model(REFERENCE_PDB, ligand_chain=None)
        lig = parse_reference_ligand()
        common = sorted(set(ca_a) & set(reference_ca or {}))
        if len(common) >= 3:
            rot, trans = kabsch(np.array([ca_a[i] for i in common]),
                                np.array([reference_ca[i] for i in common]))
        else:
            rot, trans = np.eye(3), np.zeros(3)
            print(f"  9DWV reference: only {len(common)} shared CRBN residues -- left unaligned")

        def place(points):
            if not points:
                return []
            arr = (np.array(points) @ rot.T) + trans
            return [round(float(v), 1) for v in arr.reshape(-1)]

        out["9DWV_reference"] = {
            "role": "reference",
            "a": place([ca_a[k] for k in sorted(ca_a)]),
            "b": place(ca_b),
            "el": [a[0] for a in lig],
            "lig": place([a[1] for a in lig]),
            "bonds": bonds_for(lig),
        }
        print(f"  9DWV reference: CRBN {len(ca_a)} CA, PPIL4 {len(ca_b)} CA, "
              f"ligand {len(lig)} heavy atoms (real FPFT-2216 from the mmCIF), "
              f"superposed on {len(common)} shared CRBN residues")

    os.makedirs(SITE_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    size = os.path.getsize(OUT_JSON) / 1024
    points = sum((len(v["a"]) + len(v["b"]) + len(v["lig"])) // 3 for v in out.values())
    print(f"\nWrote {os.path.relpath(OUT_JSON, SCRIPT_DIR)}: {len(out)} structures, "
          f"{points:,} points, {size:.0f} KB")
    print("All superposed on CRBN, so PPIL4's position is comparable across structures.")


if __name__ == "__main__":
    main()
