"""
Build the compound-order package for the buy list: the sheet a synthesis vendor
quotes from, plus the structures in the form they actually import.

WHAT A SYNTHESIS VENDOR NEEDS, AND WHY THE SHEET LOOKS LIKE THIS. Only 3 of the
25 molecules on this list are registered substances -- the rest are novel designs
with no CAS number and no name -- so a vendor cannot look any of them up. The
specification therefore has to be the structure itself. That means two things
must be unambiguous on every row: 分子式 (molecular formula) and 结构式
(structure), the latter as SMILES so it can be pasted straight into a drawing
package. Everything else is noise on an order form and has been dropped.

The list also contains five pairs of constitutional isomers -- same formula,
different molecule, including one pair with the Z6466608628 control -- so the
formula alone is NOT a specification here. That is why every row also carries an
InChIKey in the internal file, and why the SDF below exists at all.

OUTPUTS
  Final Buy List for Purchase.csv        the order sheet. Bilingual headers,
                                         UTF-8 with BOM so Excel opens the
                                         Chinese correctly rather than as mojibake.
  Final Buy List for Purchase.sdf        the same 25 with real 2D structures and
                                         the order fields as SD properties. This
                                         is the file to send -- a vendor imports
                                         it directly and never retypes a SMILES.
  Final Buy List Structures.svg          a printable structure sheet, for the
                                         humans reading the quote.
  Final Buy List for Purchase (annotated).csv
                                         internal superset, keeps InChIKey and
                                         provenance; feeds the project page (14).

NO INVENTED NAMES. For the 22 novel compounds there is no registry name, and
these scaffolds number awkwardly -- PubChem calls the Z6466608628 control
"3-[3-oxo-7-(4-piperazin-1-ylphenyl)-1H-isoindol-2-yl]piperidine-2,6-dione",
choosing 3-oxo/7-aryl over the naive 1-oxo/4-aryl. A hand-derived locant would
risk specifying the wrong isomer, so the name column is left empty for them and
the structure carries the specification. Only the 3 registry hits get a name,
taken from PubChem.

PURITY is a specification you REQUEST, not a property to look up: candidates get
>=95% (HPLC) for screening, controls >=98% (HPLC) because they are the
quantitative reference every candidate is judged against.

Run with the SYSTEM python (needs rdkit, certifi), after 11:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "12_annotate_purchase_list_(Ryan).py"
    # --refresh   re-query PubChem, ignoring the cache
"""
import argparse
import csv
import json
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BUY_LIST = os.path.join(SCRIPT_DIR, "final_buy_list_lock_v1_top20_(Ryan).csv")

OUT_ORDER_CSV = os.path.join(SCRIPT_DIR, "Final Buy List for Purchase.csv")
OUT_SDF = os.path.join(SCRIPT_DIR, "Final Buy List for Purchase.sdf")
OUT_SHEET_SVG = os.path.join(SCRIPT_DIR, "Final Buy List Structures.svg")
OUT_ANNOTATED = os.path.join(SCRIPT_DIR, "Final Buy List for Purchase (annotated).csv")
PUBCHEM_CACHE = os.path.join(SCRIPT_DIR, "docking_tmp", "pubchem_registry_cache.json")

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")

PURITY_CANDIDATE = ">=95% (HPLC)"
PURITY_CONTROL = ">=98% (HPLC)"
NOT_REGISTERED = "-"       # no CAS exists; the structure is the specification

ORDER_COLUMNS = [
    "编号 Code",
    "分子式 Molecular Formula",
    "结构式 Structure (SMILES)",
    "分子量 MW",
    "CAS号 CAS No.",
    "纯度级别 Purity",
]


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


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
    """Exact-structure lookup. cid None means PubChem has no record of this
    structure -- i.e. it is not a registered substance and has no CAS number."""
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


def main():
    parser = argparse.ArgumentParser(description="Build the compound-order package.")
    parser.add_argument("--refresh", action="store_true", help="re-query PubChem, ignoring the cache")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    import certifi

    rows = read_csv(BUY_LIST)
    if not rows:
        sys.exit(f"{os.path.basename(BUY_LIST)} not found or empty -- run 11 first.")

    cache = {}
    if os.path.exists(PUBCHEM_CACHE) and not args.refresh:
        cache = json.load(open(PUBCHEM_CACHE))
    ctx = ssl.create_default_context(cafile=certifi.where())

    order, annotated, mols, unregistered = [], [], [], []
    for row in rows:
        name, smiles = row["name"], row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            sys.exit(f"{name}: SMILES does not parse -- fix the buy list before ordering.")

        if smiles not in cache:
            cache[smiles] = pubchem_lookup(smiles, ctx)
            print(f"  looked up {name}")
        hit = cache[smiles]
        is_control = "control" in row["role"] or name.startswith("negctrl")

        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = round(Descriptors.MolWt(mol), 2)
        purity = PURITY_CONTROL if is_control else PURITY_CANDIDATE
        cas = hit["cas"][0] if hit["cid"] and hit["cas"] else NOT_REGISTERED
        if cas == NOT_REGISTERED:
            unregistered.append(name)

        order.append({
            "编号 Code": name,
            "分子式 Molecular Formula": formula,
            "结构式 Structure (SMILES)": smiles,
            "分子量 MW": mw,
            "CAS号 CAS No.": cas,
            "纯度级别 Purity": purity,
        })
        annotated.append({
            "name": name,
            "role": row["role"],
            "molecular_formula": formula,
            "smiles": smiles,
            "molecular_weight": mw,
            "cas_number": cas,
            "registered": "yes" if cas != NOT_REGISTERED else "no",
            "purity_grade": purity,
            "inchikey": Chem.MolToInchiKey(mol),
            # Only real registry names. See NO INVENTED NAMES in the docstring.
            "chemical_name": hit["iupac"] if hit["cid"] else "",
            "pubchem_cid": hit["cid"] or "",
        })

        AllChem.Compute2DCoords(mol)
        mol.SetProp("_Name", name)
        for key, value in (("Code", name), ("Molecular_Formula", formula), ("MW", str(mw)),
                           ("CAS", cas), ("Purity", purity),
                           ("InChIKey", Chem.MolToInchiKey(mol)),
                           ("Role", "control" if is_control else "candidate")):
            mol.SetProp(key, value)
        mols.append(mol)

    os.makedirs(os.path.dirname(PUBCHEM_CACHE), exist_ok=True)
    json.dump(cache, open(PUBCHEM_CACHE, "w"), indent=1)

    # Order sheet. utf-8-sig: without the BOM Excel reads the Chinese headers as
    # mojibake on a default Windows install, which is where these get opened.
    with open(OUT_ORDER_CSV, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=ORDER_COLUMNS)
        writer.writeheader()
        writer.writerows(order)

    with open(OUT_ANNOTATED, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(annotated[0]))
        writer.writeheader()
        writer.writerows(annotated)

    with Chem.SDWriter(OUT_SDF) as w:
        for mol in mols:
            w.write(mol)

    # Printable structure sheet: the drawn 结构式, for whoever reads the quote.
    legends = [f"{m.GetProp('Code')}  {m.GetProp('Molecular_Formula')}" for m in mols]
    svg = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(340, 290),
                               legends=legends, useSVG=True)
    svg = svg.data if hasattr(svg, "data") else str(svg)
    with open(OUT_SHEET_SVG, "w") as f:
        f.write(svg)

    print(f"\nWrote the order package for {len(order)} compounds:")
    for path in (OUT_ORDER_CSV, OUT_SDF, OUT_SHEET_SVG, OUT_ANNOTATED):
        print(f"  {os.path.basename(path):<46} {os.path.getsize(path) // 1024:>4} KB")
    print(f"\n  {len(order) - len(unregistered)}/{len(order)} have a CAS number; the other "
          f"{len(unregistered)} are novel and specified by structure alone.")
    print("  Five pairs of constitutional isomers are on this list, so the formula is not a "
          "specification by itself -- send the SDF, not just the CSV.")


if __name__ == "__main__":
    main()
