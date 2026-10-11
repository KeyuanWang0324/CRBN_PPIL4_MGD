# Session 4 — hands-on cheat sheet (instructor)

Everything you paste into the class chat is here, in the order you paste it. You do
not need to type anything off the slides. Each part has a **reference solution** —
if a student's AI ties itself in knots, paste that and move on. Every printed output
below is verified (RDKit for Part 1; the real 9DWV structure for Part 2).

**Two parts, same one idea — compare your result to a known answer:**
- **Part 1 (1:13–1:29) — is it the right MOLECULE?**  identity via InChIKey.
- **Part 2 (1:29–1:50) — is it in the right PLACE?**  interface vs the real structure.

The single rule that ties them together: **never let the software supply its own
answer key.** Part 1 breaks if the AI invents the InChIKey; Part 2 is the whole
point (HADDOCK graded itself 0.86; the real structure said 0.003).

---
---

# PART 1 — is it the right MOLECULE?

## 1a — paste the constants first

*"Copy these three lines exactly. Do not let the AI make up its own — that is the
whole point."*

```
SMILES_IN = "O=C1CCC(N2Cc3c(N)cccc3C2=O)C(=O)N1"
EXPECTED_SKELETON = "GOTYRUGSSMKFNF"
EXPECTED_FORMULA  = "C13H13N3O3"
```

- `SMILES_IN` is lenalidomide (PubChem CID 216326). Full InChIKey
  `GOTYRUGSSMKFNF-UHFFFAOYSA-N`.
- We compare only the **first block** — `GOTYRUGSSMKFNF`, the 14 characters before
  the first dash — exactly like `src/crbn/chem.py` line 71.

## 1b — paste the prompt

```
Using RDKit in Google Colab, write ONE cell.

Rules: do NOT look anything up online, and do NOT invent a SMILES or an
InChIKey. Use only the three constants I pasted above.

Start the cell with:  !pip install rdkit

Write check(smiles) that:
 1. mol = Chem.MolFromSmiles(smiles)
    If mol is None: print "INVALID SMILES - did not parse", return None.
    Never call anything else on a None mol.
 2. Print the full InChIKey, its first block (before the first '-'),
    the molecular formula, and the molecular weight.
 3. Print MATCH if that first block == EXPECTED_SKELETON, else MISMATCH
    and show the two blocks side by side.
 4. return mol.

Make the LAST line of the cell display the drawing:
    from IPython.display import display
    display(Draw.MolToImage(check(SMILES_IN)))
```

**A correct run prints** (your "did it work?" key):
```
full InChIKey : GOTYRUGSSMKFNF-UHFFFAOYSA-N
first block   : GOTYRUGSSMKFNF
formula       : C13H13N3O3
mol weight    : 259.26
MATCH
```
…and a picture of one connected molecule.

## 1c — the 2-minute proof (everyone, the moment they see MATCH)

A check that never fails is worthless. Run `check()` on this **wrong** molecule:

```
O=C1CCC(N2Cc3cccc(N)c3C2=O)C(=O)N1
```
Same formula, same weight to the decimal — the NH₂ just moved one carbon. It **must**
print **MISMATCH** (`XKAYAFBLGLCWSY`, not `GOTYRUGSSMKFNF`). Only the identity key
catches it — not the formula, not the weight, not the picture. *That* is the lesson.

## Optional Part-1 extensions
- **Not a molecule at all:** `O=C1CCC(N2Cc3c(N)ccc3C2=O)C(=O)N1` → must print
  `INVALID SMILES - did not parse` (a deleted ring atom), **not** a red crash.
- **Early finisher:** add **thalidomide** `O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1`,
  expected skeleton `UEJJHQNACJXSKW` (formula C13H10N2O4, weight 258.23 — note it is
  *not* an isomer of lenalidomide).
- **Handedness question (most likely in the room):** if a student's AI wrote the
  stereocentre, their tail differs — `GOTYRUGSSMKFNF-JTQLQIEISA-N` or `-SNVBAGLBSA-N`
  instead of `-UHFFFAOYSA-N`. Same first block = **right molecule**. Show them the
  split on the dash; that is exactly why we compare block 1 only. Do not "fix" it.

## Part-1 safety net (paste if an AI version is stuck)

```python
!pip install rdkit -q
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, inchi
from IPython.display import display

SMILES_IN = "O=C1CCC(N2Cc3c(N)cccc3C2=O)C(=O)N1"
EXPECTED_SKELETON = "GOTYRUGSSMKFNF"
EXPECTED_FORMULA  = "C13H13N3O3"

def check(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("INVALID SMILES - did not parse")
        return None
    key      = inchi.MolToInchiKey(mol)
    skeleton = key.split("-")[0]          # first block = identity
    print("full InChIKey :", key)
    print("first block   :", skeleton)
    print("formula       :", rdMolDescriptors.CalcMolFormula(mol))
    print("mol weight    :", round(Descriptors.MolWt(mol), 2))
    print("MATCH" if skeleton == EXPECTED_SKELETON else "MISMATCH")
    return mol

display(Draw.MolToImage(check(SMILES_IN)))
# check("O=C1CCC(N2Cc3cccc(N)c3C2=O)C(=O)N1")   # wrong molecule -> MISMATCH
```

## Part-1 snag triage
| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: rdkit` | forgot install | first line `!pip install rdkit` |
| `NameError: Draw` | AI didn't import it | `from rdkit.Chem import Draw` |
| nothing shows | used `MolToFile`/`img.show()` | last line `display(Draw.MolToImage(...))` |
| correct molecule prints MISMATCH | comparing FULL InChIKey | compare `key.split("-")[0]` only |
| big red `ArgumentError` | called RDKit on a `None` mol | `if mol is None:` guard first |
| AI fetches from PubChem | "helpful" web lookup | delete it — check must not depend on the thing it checks |
| everything MATCHes/MISMATCHes | AI wrapped a bare `except:` | delete it — it hides the real result |

---
---

# PART 2 — is it in the right PLACE?

*Same move — compare to the answer key — but for the POSE, not the molecule. These
are the interface residues (which amino acids touch) from your OWN top docking model,
`cand_5`, versus the real 9DWV structure.*

## 2a — paste the four residue sets

*"Copy these exactly — same rule as before, the AI does not get to invent them."*

```
# Where the two proteins actually touch in the real structure (PDB 9DWV):
REAL_CRBN  = {351, 353, 355, 357, 372, 373, 386, 388, 397, 400}
REAL_PPIL4 = {249, 250, 273, 275, 276, 277, 278, 279}

# Where YOUR top docking model (cand_5) made them touch:
YOUR_CRBN  = {325, 351, 352, 353, 355, 369, 370, 371, 372, 373, 375, 377, 378, 379,
              380, 386, 387, 388, 390, 392, 393, 394, 395, 396, 397, 400}
YOUR_PPIL4 = {46, 49, 50, 52, 60, 61, 62, 63, 64, 74, 77, 79, 80, 82, 85, 97, 98, 99,
              100, 101, 102, 103, 107, 118, 119, 122, 123, 146, 147, 148}
```

## 2b — paste the prompt

```
Using the 4 residue sets I pasted above, write ONE Google Colab cell.

Write agree(your_set, real_set, label) that:
 1. found = your_set & real_set
 2. print the label and how many your model found:
    len(found) out of len(real_set)
 3. print the matched residues, sorted
 4. print FOUND THE SPOT if found is at least half of real_set,
    otherwise WRONG SPOT

Then call it twice:
    agree(YOUR_CRBN,  REAL_CRBN,  "CRBN")
    agree(YOUR_PPIL4, REAL_PPIL4, "PPIL4")
```

**A correct run prints:**
```
CRBN : found 9 of 10 real contact residues -> FOUND THE SPOT
   matched: [351, 353, 355, 372, 373, 386, 388, 397, 400]
PPIL4: found 0 of 8 real contact residues -> WRONG SPOT
   matched: []
```
Celebrate the CRBN 9/10 first — that is a real, correct result. Then the PPIL4 0/8.

## 2c — the self-grade trap (do this LIVE — it's the heart of the day)

```
agree(YOUR_CRBN, YOUR_CRBN, "CRBN vs ITSELF")   ->  found 26 of 26 -> FOUND THE SPOT
```
*"You compared your answer to itself and scored a perfect 100%. Useless, right? That
is exactly what HADDOCK did when it graded itself DockQ 0.86. A program grading its
own homework always passes. Only the real structure gave an honest grade."*

## 2d — the root cause (why exactly zero) — this is Session 5

```
REFERENCE_HAS = set(range(240, 314))   # 9DWV only resolves PPIL4 240-313 (RRM domain)
print("your PPIL4 contacts that even EXIST in 9DWV:",
      sorted(YOUR_PPIL4 & REFERENCE_HAS) or "NONE")
print("your PPIL4 contacts span residues",
      min(YOUR_PPIL4), "to", max(YOUR_PPIL4))
```
Prints `NONE` and `46 to 148`. The real glue grabs PPIL4's **RRM domain (240–313)**;
your model used residues **46–148**. You docked the **whole** 492-residue protein and
it grabbed the wrong end — the real structure does not even contain those residues,
so nothing could match. **Fix (Session 5): dock only the RRM domain, then re-score.**

## Part-2 safety net (paste if an AI version is stuck)

```python
REAL_CRBN  = {351, 353, 355, 357, 372, 373, 386, 388, 397, 400}
REAL_PPIL4 = {249, 250, 273, 275, 276, 277, 278, 279}
YOUR_CRBN  = {325, 351, 352, 353, 355, 369, 370, 371, 372, 373, 375, 377, 378, 379,
              380, 386, 387, 388, 390, 392, 393, 394, 395, 396, 397, 400}
YOUR_PPIL4 = {46, 49, 50, 52, 60, 61, 62, 63, 64, 74, 77, 79, 80, 82, 85, 97, 98, 99,
              100, 101, 102, 103, 107, 118, 119, 122, 123, 146, 147, 148}

def agree(your_set, real_set, label):
    found = your_set & real_set
    verdict = "FOUND THE SPOT" if len(found) >= len(real_set) / 2 else "WRONG SPOT"
    print(f"{label}: found {len(found)} of {len(real_set)} real contact residues -> {verdict}")
    print(f"   matched: {sorted(found)}")

agree(YOUR_CRBN,  REAL_CRBN,  "CRBN ")
agree(YOUR_PPIL4, REAL_PPIL4, "PPIL4")
```

## Part-2 snag triage
| Symptom | Cause | Fix |
|---|---|---|
| AI writes a for-loop over residues | didn't use sets | fine if it works, but nudge to `your_set & real_set` — shared members IS the idea |
| AI recomputes sets from a structure file | invented data | stop it — use the numbers we pasted, like Part 1 |
| nonzero PPIL4 overlap | a mistyped set | re-paste the sets from this file |
| `TypeError: unsupported operand & ` | pasted lists, not sets | wrap in `set(...)` or keep the `{ }` braces |

---

**Where these numbers came from** (if a student asks): `figures/make_ternary_figures.py`
measured the interface residues (heavy atoms within 5 Å across the two chains) directly
from the students' own `09_best_model_cand_5_(Ryan).pdb` and the 9DWV reference. The
PASS/MARGINAL/FAIL bands on slide 8 come from `docking/4_overlap.py` line 31; the
first-block InChIKey comparison from `src/crbn/chem.py` line 71.
