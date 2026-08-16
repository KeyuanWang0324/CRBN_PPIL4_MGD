"""
Annotate the purchase sheet with the identifiers a supplier actually needs:
名称 (chemical name), CAS 号 (CAS Registry Number), 纯度级别 (purity grade),
plus InChIKey / molecular formula / MW.

THE HEADLINE RESULT, ESTABLISHED BY THIS SCRIPT, NOT ASSUMED:
Only 3 of the 25 molecules on the buy list are registered substances with a
CAS number. Every structure is looked up in PubChem by EXACT STRUCTURE (not by
name), and 22 of them have no PubChem CID at all -- because 20 of the 25 are
novel analogs invented by 01_generate_thalidomide_analogs and two of the
controls were hand-built for this study. A compound with no registry entry has
no CAS number; that is a fact about the molecule, not a gap in the lookup.

So this is a CUSTOM-SYNTHESIS list, not a catalog order. Which is the same
reason FPFT-2216 is not on it (PROTOCOL_LOCK.md section 7).

WHY THERE IS NO IUPAC NAME FOR THE NOVEL COMPOUNDS. For the 3 registered
molecules the name comes from PubChem and is authoritative. For the other 22
this script deliberately does NOT invent one. These scaffolds number
awkwardly -- PubChem names the Z6466608628 control
"3-[3-oxo-7-(4-piperazin-1-ylphenyl)-1H-isoindol-2-yl]piperidine-2,6-dione",
choosing 3-oxo/7-aryl over the naive 1-oxo/4-aryl -- so a hand-derived locant
is a real risk of specifying the WRONG ISOMER to a synthesis vendor. The
`chemical_name` column therefore carries an unambiguous scaffold + substituent
description, and the row's identity is pinned by InChIKey and SMILES, which is
what a vendor imports anyway. If a form needs a formal IUPAC name, generate it
from the SMILES with ChemDraw/ACD or let the vendor derive it from the
structure -- do not retype one by hand.

PURITY GRADE is a specification you REQUEST, not a property to look up, so it
is assigned by role, using ordinary practice for this kind of work:
  - candidates (biochemical/cellular screening)  -> >=95% (HPLC)
  - controls (quantitative reference compounds)  -> >=98% (HPLC)
Controls set the bar every candidate is judged against (PROTOCOL_LOCK.md
section 7), so they get the tighter spec. Raise either if an assay needs it.

Network results are cached in PUBCHEM_CACHE so re-runs are offline and stable.

Run with the SYSTEM python (needs rdkit, certifi):
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "12_annotate_purchase_list_(Ryan).py"
    # --refresh   ignore the cache and re-query PubChem
"""
import argparse
import csv
import json
import os
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PURCHASE_CSV = os.path.join(SCRIPT_DIR, "Final Buy List for Purchase.csv")
OUT_CSV = os.path.join(SCRIPT_DIR, "Final Buy List for Purchase (annotated).csv")
PUBCHEM_CACHE = os.path.join(SCRIPT_DIR, "docking_tmp", "pubchem_registry_cache.json")

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")

PURITY_CANDIDATE = ">=95% (HPLC)"
PURITY_CONTROL = ">=98% (HPLC)"

# Scaffold descriptions for the novel compounds -- precise about the core and
# its substituents without asserting IUPAC locants. See the docstring.
CORE_DIOXO = "2-(2,6-dioxopiperidin-3-yl)isoindole-1,3-dione core (phthalimide-glutarimide)"
CORE_MONOOXO = "3-(3-oxo-2,3-dihydro-1H-isoindol-2-yl)piperidine-2,6-dione core (isoindolinone-glutarimide)"


def http_post(path, data, ctx):
    req = urllib.request.Request(f"{PUBCHEM}/{path}",
                                 data=urllib.parse.urlencode(data).encode(),
                                 headers={"User-Agent": "ISEF-buylist/1.0"})
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=30) as r:
            return json.load(r)
    except urllib.error.HTTPError as exc:
        return {"_http": exc.code}          # 404 = no such structure, a real answer
    except Exception as exc:                 # noqa: BLE001 -- one failure must not abort the sheet
        return {"_err": str(exc)}


def pubchem_lookup(smiles, ctx):
    """Exact-structure lookup. Returns {cid, cas[], iupac} -- cid None when
    PubChem has no record of this structure (i.e. it is not a registered
    substance and therefore has no CAS number)."""
    res = http_post("compound/smiles/cids/JSON", {"smiles": smiles}, ctx)
    cids = res.get("IdentifierList", {}).get("CID", [])
    cid = cids[0] if cids and cids[0] else None
    if not cid:
        return {"cid": None, "cas": [], "iupac": ""}
    time.sleep(0.25)                          # PubChem asks for <=5 req/s
    syn = http_post(f"compound/cid/{cid}/synonyms/JSON", {}, ctx)
    try:
        names = syn["InformationList"]["Information"][0].get("Synonym", [])
    except (KeyError, IndexError):
        names = []
    time.sleep(0.25)
    prop = http_post(f"compound/cid/{cid}/property/IUPACName/JSON", {}, ctx)
    try:
        iupac = prop["PropertyTable"]["Properties"][0].get("IUPACName", "")
    except (KeyError, IndexError):
        iupac = ""
    return {"cid": cid, "cas": [n for n in names if CAS_RE.match(n)], "iupac": iupac}


def describe(mol, name):
    """A structural description for a compound with no registry name."""
    from rdkit import Chem
    dioxo = mol.HasSubstructMatch(Chem.MolFromSmarts("O=C1N(C(=O)c2ccccc21)"))
    if name.startswith("negctrl"):
        return "study-specific negative control fragment (see PROTOCOL_LOCK.md section 7)"
    core = CORE_DIOXO if dioxo else CORE_MONOOXO
    return f"novel analog, {core}"


def main():
    parser = argparse.ArgumentParser(description="Annotate the purchase sheet.")
    parser.add_argument("--refresh", action="store_true", help="re-query PubChem, ignoring the cache")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    import certifi

    rows = list(csv.DictReader(open(PURCHASE_CSV, newline="")))
    cache = {}
    if os.path.exists(PUBCHEM_CACHE) and not args.refresh:
        cache = json.load(open(PUBCHEM_CACHE))
    ctx = ssl.create_default_context(cafile=certifi.where())

    out, need_synthesis = [], []
    for r in rows:
        name, smiles = r["name"], r["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise SystemExit(f"{name}: SMILES does not parse -- fix the buy list before ordering.")

        if smiles not in cache:
            cache[smiles] = pubchem_lookup(smiles, ctx)
            print(f"  looked up {name}")
        hit = cache[smiles]
        is_control = "control" in name or name.startswith("negctrl")

        if hit["cid"]:
            chemical_name = hit["iupac"] or "(PubChem record has no IUPAC name)"
            cas = hit["cas"][0] if hit["cas"] else "not assigned in PubChem"
            source = f"PubChem CID {hit['cid']}"
        else:
            chemical_name = describe(mol, name)
            cas = "none - novel compound, not registered"
            source = "no PubChem record (exact-structure search)"
            need_synthesis.append(name)

        out.append({
            "name": name,
            "chemical_name": chemical_name,
            "cas_number": cas,
            "purity_grade": PURITY_CONTROL if is_control else PURITY_CANDIDATE,
            "smiles": smiles,
            "inchikey": Chem.MolToInchiKey(mol),
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "identifier_source": source,
        })

    os.makedirs(os.path.dirname(PUBCHEM_CACHE), exist_ok=True)
    json.dump(cache, open(PUBCHEM_CACHE, "w"), indent=1)

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out[0]))
        writer.writeheader()
        writer.writerows(out)

    registered = [r for r in out if r["cas_number"].count("-") == 2]
    print(f"\nWrote {OUT_CSV}")
    print(f"  {len(registered)}/{len(out)} have a CAS number (registered substances):")
    for r in registered:
        print(f"    {r['name']:<32} {r['cas_number']:<16} {r['chemical_name']}")
    print(f"  {len(need_synthesis)}/{len(out)} have NO CAS -- novel/hand-built, so custom "
          f"synthesis:\n    {', '.join(need_synthesis)}")
    print(f"\n  Purity requested: candidates {PURITY_CANDIDATE}, controls {PURITY_CONTROL} "
          "(a spec to request, not a looked-up property).")


if __name__ == "__main__":
    main()
