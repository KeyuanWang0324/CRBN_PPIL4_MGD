"""
Export PyMOL's own cartoon geometry for every structure on the project page as
binary glTF, so the page shows the actual PyMOL mesh rather than a viewer's
reconstruction of it.

WHAT THIS GIVES YOU THAT A PDB DOES NOT. Handing a PDB to a molecular viewer
means that viewer builds its own cartoon -- its own ribbon framing, its own
helix smoothing, its own idea of where a strand ends. Exporting glTF instead
ships the triangles PyMOL itself generated, with PyMOL's vertex colours baked
in. What renders in the browser is then the same surface PyMOL would ray-trace,
not an approximation of it.

SIZE IS THE WHOLE ENGINEERING PROBLEM, and it is why the settings below are what
they are. Measured on one structure:

    default (cartoon + sticks + spheres)      5.69 MB
    cartoon_sampling 3                        4.62 MB
    ... + sphere_quality 0                    2.08 MB
    ... + sticks only, no spheres             1.58 MB
    ... written as binary glTF instead of     1.00 MB
        base64-embedded .gltf

Two of those matter most. Sphere geometry dominates a scene like this, so the
ligand is drawn as sticks alone. And PyMOL writes .gltf with its buffer inlined
as a base64 data URI, which inflates binary by a third for nothing -- this
script repacks the same data as .glb, which is what the 1.58 -> 1.00 MB step is.
Nothing is resampled in that step; it is purely a container change.

`cartoon_sampling 3` is the one genuine quality trade: it halves the number of
segments PyMOL interpolates along the backbone. At the sizes these are viewed
it is not visible, and it is the difference between a 27 MB site and a 43 MB one.

EVERYTHING IS SUPERPOSED ON CRBN, as with the renders and the PDBs, so PPIL4's
position stays comparable between structures.

Run under PyMOL, from the project directory:
    /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq "16_export_pymol_gltf_(Ryan).py"
"""
import base64
import csv
import hashlib
import json
import os
import struct
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
GLB_DIR = os.path.join(SITE_DIR, "glb")
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

# Page palette, shared with 13's renders and the legend swatches.
COLOURS = {
    # CRBN is the ANCHOR: every structure is superposed on it, so it is identical
    # in all 27 views and carries no information. It is therefore neutral grey --
    # context, then out of the way. PPIL4's position and the ligand are what vary,
    # so they carry the colour.
    #
    # Deliberately dark and low-chroma. These are large filled ribbons: a value
    # that looks right on a legend chip is glaring across a whole protein domain.
    # Blue and orange still sit on the standard colour-vision-safe axis, so they
    # stay separable under deuteranopia and protanopia at this chroma. They are
    # baked into the mesh as vertex colours, so one set has to hold up on both the
    # light and the dark theme -- which is the floor on how dark these can go.
    "crbn":   (0.486, 0.522, 0.514),   # #7C8583 deep grey        - the fixed anchor
    "ppil4":  (0.275, 0.408, 0.549),   # #46688C deep slate blue  - varies between structures
    "ligand": (0.663, 0.424, 0.243),   # #A96C3E deep terracotta  - the molecule being judged
}


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


def gltf_to_glb(gltf_path, glb_path):
    """Repack a base64-embedded .gltf as binary .glb.

    PyMOL inlines the mesh buffer as a data URI, so a third of the file is
    base64 padding. This decodes it back to bytes and writes the standard GLB
    container -- same triangles, same colours, ~36% smaller and faster to parse
    because the browser skips the base64 decode."""
    with open(gltf_path) as f:
        doc = json.load(f)

    blobs = []
    for buf in doc.get("buffers", []):
        uri = buf.get("uri", "")
        if not uri.startswith("data:"):
            return None                      # external .bin: leave it alone
        blobs.append(base64.b64decode(uri.split(",", 1)[1]))
    if len(blobs) != 1:
        return None                          # multi-buffer glTF needs re-indexing

    blob = blobs[0]
    doc["buffers"][0].pop("uri", None)
    doc["buffers"][0]["byteLength"] = len(blob)

    js = json.dumps(doc, separators=(",", ":")).encode()
    js += b" " * ((4 - len(js) % 4) % 4)          # chunks must be 4-byte aligned
    blob += b"\0" * ((4 - len(blob) % 4) % 4)

    with open(glb_path, "wb") as f:
        f.write(b"glTF" + struct.pack("<II", 2, 12 + 8 + len(js) + 8 + len(blob)))
        f.write(struct.pack("<II", len(js), 0x4E4F534A) + js)     # JSON chunk
        f.write(struct.pack("<II", len(blob), 0x004E4942) + blob)  # BIN chunk
    return os.path.getsize(glb_path)


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

    os.makedirs(GLB_DIR, exist_ok=True)
    for name, rgb in COLOURS.items():
        cmd.set_color(f"pg_{name}", list(rgb))

    # See the docstring: sampling 3 halves the backbone interpolation, and the
    # ligand is sticks-only because sphere geometry dominates the mesh.
    cmd.set("cartoon_sampling", 3)
    cmd.set("stick_quality", 8)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_opaque_background", 0)
    # Lighting is kept flat and matte on purpose. Cartoon ribbons are large
    # smooth surfaces, and PyMOL's default specular turns each one into a
    # highlight that reads as glare regardless of how dark the base colour is.
    cmd.set("specular", 0.05)
    cmd.set("shininess", 10)
    cmd.set("ambient", 0.10)
    cmd.set("direct", 0.28)
    cmd.set("reflect", 0.28)
    cmd.set("light_count", 2)

    index, reference_obj, total = {}, None, 0
    tmp_gltf = os.path.join(GLB_DIR, "_tmp.gltf")

    def export(obj, display, role):
        nonlocal total
        cmd.hide("everything")
        cmd.show("cartoon", f"{obj} and (chain A or chain B)")
        cmd.show("sticks", f"{obj} and chain C")
        cmd.color("pg_crbn", f"{obj} and chain A")
        cmd.color("pg_ppil4", f"{obj} and chain B")
        cmd.color("pg_ligand", f"{obj} and chain C")
        cmd.util.cnc(f"{obj} and chain C and not elem C")
        cmd.disable("all")
        cmd.enable(obj)
        if os.path.exists(tmp_gltf):
            os.remove(tmp_gltf)
        cmd.save(tmp_gltf)
        glb_path = os.path.join(GLB_DIR, f"{display}.glb")
        size = gltf_to_glb(tmp_gltf, glb_path)
        if size is None:
            print(f"  {display}: unexpected glTF layout -- keeping .gltf")
            os.replace(tmp_gltf, os.path.join(GLB_DIR, f"{display}.gltf"))
            size = os.path.getsize(os.path.join(GLB_DIR, f"{display}.gltf"))
            url = f"glb/{display}.gltf"
        else:
            os.remove(tmp_gltf)
            url = f"glb/{display}.glb"
        with open(os.path.join(GLB_DIR, os.path.basename(url)), "rb") as f:
            version = hashlib.sha256(f.read()).hexdigest()[:10]
        index[display] = {"role": role, "mesh": f"{url}?v={version}"}
        total += size
        print(f"  {display:<32} {size / 1e6:5.2f} MB")

    loaded = []
    for display, run_name, role in entries:
        path = top_model_path(run_name)
        if path is None:
            print(f"  {display}: no rank-1 model -- skipped")
            continue
        obj = "s%d" % len(loaded)
        cmd.load(path, obj)
        cmd.dss(obj)
        if reference_obj is None:
            reference_obj = obj
        else:
            cmd.super(f"{obj} and chain A and name CA", f"{reference_obj} and chain A and name CA")
        loaded.append((obj, display, role))

    if os.path.exists(REFERENCE_PDB) and reference_obj:
        obj = "sref"
        cmd.load(REFERENCE_PDB, obj)
        cmd.dss(obj)
        cmd.super(f"{obj} and chain A and name CA", f"{reference_obj} and chain A and name CA")
        atoms = reference_ligand_atoms()
        if atoms:
            tmp = os.path.join(GLB_DIR, "_lig.pdb")
            with open(tmp, "w") as f:
                for i, (aname, el, (x, y, z)) in enumerate(atoms, 1):
                    f.write(f"HETATM{i:5d} {aname:<4s} LIG C 900    "
                            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el[:2]:>2s}\n")
                f.write("END\n")
            cmd.load(tmp, "sreflig")
            cmd.transform_object("sreflig", cmd.get_object_matrix(obj))
            cmd.create(obj, f"{obj} or sreflig", zoom=0)
            cmd.delete("sreflig")
            os.remove(tmp)
        loaded.append((obj, "9DWV_reference", "reference"))

    # Put the ligand at the origin so the page's viewer orbits the molecule
    # rather than the midpoint of a sprawling AlphaFold model. Crucially this is
    # ONE translation -- the reference's ligand centre -- applied to every
    # structure, not each structure's own centre. Re-centring each on its own
    # ligand would shift them relative to each other by a couple of angstroms
    # and quietly destroy the shared CRBN frame that makes them comparable.
    if reference_obj:
        cx, cy, cz = cmd.centerofmass(f"{reference_obj} and chain C")
        shift = [1, 0, 0, -cx, 0, 1, 0, -cy, 0, 0, 1, -cz, 0, 0, 0, 1]
        for obj, _, _ in loaded:
            cmd.transform_object(obj, shift)
        print(f"Ligand centred: all structures shifted by "
              f"({-cx:.1f}, {-cy:.1f}, {-cz:.1f}) A -- one shared translation.")

    for obj, display, role in loaded:
        export(obj, display, role)

    with open(INDEX_JSON, "w") as f:
        json.dump(index, f, separators=(",", ":"), indent=0)
    print(f"\nWrote {len(index)} meshes into {os.path.relpath(GLB_DIR, SCRIPT_DIR)}/ "
          f"({total / 1e6:.1f} MB total, fetched one at a time on demand)")
    print("These are PyMOL's own triangles and vertex colours, superposed on CRBN.")


main()
