"""
Render one PyMOL image of the rank-1 ternary model for every molecule on the
buy list, for the project page's per-candidate detail panels.

THE POINT OF THIS SCRIPT IS THE SHARED CAMERA. Rendering each complex with its
own `orient` would produce 20 images in 20 arbitrary frames, which look
comparable but are not -- a candidate whose PPIL4 sits somewhere quite
different would simply be re-oriented until it looked like the others. So every
model is first superposed onto ONE reference by its CRBN chain, then all 20 are
rendered through a single stored view matrix. CRBN is therefore in the same
place in every image, and where PPIL4 lands is a real difference you can see.

Chains follow this project's convention throughout: A = CRBN, B = PPIL4,
C = the ligand (resnum 900).

Colours match the project page's palette so the images sit in the design
rather than next to it: CRBN teal, PPIL4 indigo, ligand ochre. Backgrounds are
rendered transparent, so one image works on both the light and dark theme.

A note on what these images honestly show: PPIL4 here is an AlphaFold model and
its low-confidence loops sprawl. They are NOT trimmed for prettiness -- what is
rendered is the whole chain that was actually docked.

Run under PyMOL (the pip wheel ships broken rpaths on macOS; the app works):
    /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq "13_render_ternary_views_(Ryan).py"
"""
import csv
import glob
import gzip
import os
import shutil
import sys
import tempfile

BUY_LIST_NAME = "final_buy_list_lock_v1_top20_(Ryan).csv"


def _project_dir():
    """PyMOL's `run` rebinds __file__ to its own package directory, so the usual
    dirname(__file__) lands inside PyMOL.app. Fall back to the working directory
    (and then to $ISEF_DIR) -- whichever actually holds the buy list."""
    for candidate in (os.path.dirname(os.path.abspath(__file__)),
                      os.getcwd(),
                      os.environ.get("ISEF_DIR", "")):
        if candidate and os.path.exists(os.path.join(candidate, BUY_LIST_NAME)):
            return candidate
    sys.exit(f"Cannot locate {BUY_LIST_NAME}. Run this from the project directory, "
             "or set $ISEF_DIR to it.")


SCRIPT_DIR = _project_dir()
BUY_LIST = os.path.join(SCRIPT_DIR, BUY_LIST_NAME)
OUT_DIR = os.path.join(SCRIPT_DIR, "CRBN_Project_site", "img")

RUN_DIR_BASES = [
    os.path.join(os.path.expanduser("~"), "haddock_runs", "haddock3_ternary_with_ligand_run"),
    os.path.join(SCRIPT_DIR, "docking_tmp", "haddock3_ternary_with_ligand_run"),
]

WIDTH, HEIGHT, DPI = 620, 460, 120

# Page palette (see CRBN_Project_site/index.html).
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
    """The rank-1 model of this molecule's final caprieval, decompressed to a
    temp file if it is gzipped (PyMOL will not read .pdb.gz directly here)."""
    for base in RUN_DIR_BASES:
        capri = os.path.join(base, run_name, "run1", "9_caprieval")
        tsv = os.path.join(capri, "capri_ss.tsv")
        if not os.path.exists(tsv):
            continue
        with open(tsv, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        top = next((r for r in rows if r["caprieval_rank"] == "1"), None)
        if top is None:
            continue
        resolved = os.path.normpath(os.path.join(capri, top["model"]))
        if os.path.exists(resolved):
            return resolved
        if os.path.exists(resolved + ".gz"):
            tmp = os.path.join(tempfile.mkdtemp(), os.path.basename(resolved))
            with gzip.open(resolved + ".gz", "rb") as src, open(tmp, "wb") as dst:
                shutil.copyfileobj(src, dst)
            return tmp
    return None


def main():
    from pymol import cmd

    with open(BUY_LIST, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["role"] in ("candidate", "crosscheck")]
    if not rows:
        sys.exit(f"No candidates in {os.path.basename(BUY_LIST)} -- run 11 first.")
    os.makedirs(OUT_DIR, exist_ok=True)

    for name, rgb in COLOURS.items():
        cmd.set_color(f"pg_{name}", list(rgb))

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
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("stick_radius", 0.22)
    cmd.set("sphere_scale", 0.26)
    cmd.set("orthoscopic", 1)

    loaded, reference = [], None
    for row in rows:
        run_name, display = row["run_name"], row["name"]
        path = top_model_path(run_name)
        if path is None:
            print(f"  {display}: no rank-1 model found -- skipped")
            continue
        obj = "m_" + display
        cmd.load(path, obj)
        # Superpose on CRBN only. `super` is sequence-independent, which matters
        # because the ligand differs between models and would otherwise drag the fit.
        if reference is None:
            reference = obj
        else:
            cmd.super(f"{obj} and chain A and name CA", f"{reference} and chain A and name CA")
        loaded.append((display, obj))

    if not loaded:
        sys.exit("Nothing loaded -- no rank-1 models found under any run root.")

    cmd.hide("everything")
    cmd.show("cartoon", "chain A or chain B")
    cmd.show("sticks", "chain C")
    cmd.show("spheres", "chain C")
    cmd.color("pg_crbn", "chain A")
    cmd.color("pg_ppil4", "chain B")
    cmd.color("pg_ligand", "chain C")
    cmd.util.cnc("chain C and not elem C")   # heteroatoms keep their element colours

    # One camera for all of them: framed on the reference's CRBN + ligand (the
    # degron pocket, the part every model shares), then pulled back far enough
    # that the whole superposed ensemble stays inside the frame. Because the
    # models are already superposed, this same view is valid for every one.
    # Frame on the part that is actually shared and ordered: CRBN, the ligand,
    # and the PPIL4 residues near the ligand. Zooming to `all` instead makes the
    # frame the union bounding box of 20 superposed AlphaFold models, whose
    # low-confidence tails sprawl in different directions -- that box is mostly
    # empty space and leaves each complex small and off-centre. The selection
    # below spans every loaded object, so the camera still fits them all; the
    # far tails simply run past the edge of the frame, as they would in any
    # published figure of this complex.
    # The ligand goes dead centre. Zooming on a wider selection centres the frame
    # on that selection's bounding box instead, which puts the ligand off to one
    # side -- and the ligand is the subject of every one of these images. Zooming
    # on chain C with a large buffer keeps it centred while pulling back far
    # enough to show both proteins around it. The far AlphaFold tails run past
    # the edge, as they would in any published figure of this complex.
    cmd.orient(f"{reference} and (chain A or chain C)")
    cmd.zoom(f"{reference} and chain C", buffer=26.0, complete=1)
    view = cmd.get_view()

    print(f"Rendering {len(loaded)} views at {WIDTH}x{HEIGHT} into "
          f"{os.path.relpath(OUT_DIR, SCRIPT_DIR)}/ ...")
    for display, obj in loaded:
        cmd.disable("all")
        cmd.enable(obj)
        cmd.set_view(view)          # re-assert after enable/disable
        out = os.path.join(OUT_DIR, f"{display}.png")
        cmd.ray(WIDTH, HEIGHT)
        cmd.png(out, dpi=DPI)
        size = os.path.getsize(out) // 1024 if os.path.exists(out) else 0
        print(f"  {display:<12} {size:>4} KB")

    total = sum(os.path.getsize(p) for p in glob.glob(os.path.join(OUT_DIR, "*.png")))
    print(f"\nWrote {len(loaded)} images, {total / 1024 / 1024:.1f} MB total.")
    print("All rendered through one shared camera after superposition on CRBN -- "
          "so CRBN is in the same place in every image and PPIL4's position is comparable.")


main()
