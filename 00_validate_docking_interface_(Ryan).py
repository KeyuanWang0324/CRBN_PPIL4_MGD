"""
Compare cand_5's own top HADDOCK3 docking model's CRBN/PPIL4 contact
residues against the real interface from PDB 9DWV -- same idea as
00_validate_reference_smiles_(Ryan).py (compare to a known answer instead
of trusting the software's own self-grade), but for the POSE instead of
the molecule.

All four residue sets below are pasted verbatim from session4_handson.md
(Part 2a) -- not recomputed here, not looked up, on purpose: the whole
point of agree() is to compare against a known answer, so this script
must not get to invent that answer.

  REAL_CRBN / REAL_PPIL4: where CRBN and PPIL4 actually touch in the real
      cryo-EM structure (PDB 9DWV).
  YOUR_CRBN / YOUR_PPIL4: where cand_5's own top docking model made them
      touch (07_best_model_cand_5_(Ryan).pdb).

Run with the SYSTEM python:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "00_validate_docking_interface_(Ryan).py"
"""

# Where the two proteins actually touch in the real structure (PDB 9DWV):
REAL_CRBN = {351, 353, 355, 357, 372, 373, 386, 388, 397, 400}
REAL_PPIL4 = {249, 250, 273, 275, 276, 277, 278, 279}

# Where cand_5's own top docking model made them touch:
YOUR_CRBN = {325, 351, 352, 353, 355, 369, 370, 371, 372, 373, 375, 377, 378, 379,
             380, 386, 387, 388, 390, 392, 393, 394, 395, 396, 397, 400}
YOUR_PPIL4 = {46, 49, 50, 52, 60, 61, 62, 63, 64, 74, 77, 79, 80, 82, 85, 97, 98, 99,
              100, 101, 102, 103, 107, 118, 119, 122, 123, 146, 147, 148}


def agree(your_set, real_set, label):
    found = your_set & real_set
    print(f"{label}: found {len(found)} of {len(real_set)} real contact residues")
    print(f"  matched: {sorted(found)}")
    if len(found) >= len(real_set) / 2:
        print("FOUND THE SPOT")
    else:
        print("WRONG SPOT")


agree(YOUR_CRBN, REAL_CRBN, "CRBN")
agree(YOUR_PPIL4, REAL_PPIL4, "PPIL4")
