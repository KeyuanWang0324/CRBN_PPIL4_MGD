"""
Sanity-check a reference glue's SMILES against known-good identifiers
before it's trusted anywhere downstream (e.g. as a positive control in
02/03's classifier training, or a scaffold seed in 01) -- confirms RDKit
parses it into the expected molecule instead of silently trusting a
possibly-corrupted copy-paste.

SMILES_IN is kept exactly as pasted, uncorrected -- the whole point of
check() is to catch corruption like stray spaces or a "0" (digit) swapped
in for an "O" (oxygen), so "fixing" it here would defeat the check.

Ground truth (from the constants below only -- no lookups, no invented
values):
  - EXPECTED_SKELETON: lenalidomide's InChIKey first block, PubChem CID
    216326, full key GOTYRUGSSMKFNF-UHFFFAOYSA-N.
  - EXPECTED_FORMULA: lenalidomide's molecular formula.
  (Thalidomide's full InChIKey block starts UEJJHQNACJXSKW -- a different
  skeleton, noted only for contrast, not used in the check.)

Run with the SYSTEM python (has rdkit installed):
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "00_validate_reference_smiles_(Ryan).py"
"""
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors

SMILES_IN = "O=C1CCC(N2Cc3c (N) cccc3C2=O)C(=O) N1"
EXPECTED_SKELETON = "GOTYRUGSSMKFNF"
EXPECTED_FORMULA = "C13H13N3O3"


def check(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("INVALID SMILES - did not parse")
        return None

    inchikey = Chem.MolToInchiKey(mol)
    skeleton = inchikey.split("-")[0]
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)

    print(f"InChIKey: {inchikey}")
    print(f"Skeleton (first block): {skeleton}")
    print(f"Formula: {formula}")
    print(f"Molecular weight: {mw:.2f}")

    if skeleton == EXPECTED_SKELETON:
        print("MATCH")
    else:
        print("MISMATCH")
        print(f"  got:      {skeleton}")
        print(f"  expected: {EXPECTED_SKELETON}")

    if formula == EXPECTED_FORMULA:
        print("FORMULA MATCH")
    else:
        print("FORMULA MISMATCH")
        print(f"  got:      {formula}")
        print(f"  expected: {EXPECTED_FORMULA}")

    return mol


def show_if_parsed(mol):
    if mol is not None:
        # PIL's own .show() instead of IPython.display.display -- this is
        # a plain script (not a notebook), and .show() opens the image in
        # the OS's default viewer either way without an extra IPython
        # dependency.
        Draw.MolToImage(mol).show()
    else:
        print("No image to display -- check() returned None.")


print("=== SMILES_IN, as pasted (uncorrected) ===")
show_if_parsed(check(SMILES_IN))

# Not one of the three original constants -- added separately, on request,
# purely to confirm the corrected string actually resolves to lenalidomide
# once the copy-paste corruption (stray spaces, "0" for "O") is removed.
CORRECTED_SMILES = "O=C1CCC(N2Cc3c(N)cccc3C2=O)C(=O)N1"
print()
print("=== CORRECTED_SMILES (stray spaces removed, 'O' not '0') ===")
show_if_parsed(check(CORRECTED_SMILES))
