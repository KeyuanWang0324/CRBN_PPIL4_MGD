# minimal xlsx reader - openpyxl is not installed
import zipfile, re
import xml.etree.ElementTree as ET
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def col_to_idx(ref):
    c=re.match(r'([A-Z]+)', ref).group(1); n=0
    for ch in c: n=n*26+(ord(ch)-64)
    return n-1
def read(path):
    z=zipfile.ZipFile(path)
    shared=[]
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall('m:si',NS):
            shared.append("".join(t.text or "" for t in si.iter('{%s}t'%NS['m'])))
    rels={r.get('Id'):r.get('Target') for r in
          ET.fromstring(z.read("xl/_rels/workbook.xml.rels")).findall('*')}
    out={}
    for sh in ET.fromstring(z.read("xl/workbook.xml")).find('m:sheets',NS):
        name=sh.get('name'); tgt=rels[sh.get('{%s}id'%NS['r'])].lstrip('/')
        if not tgt.startswith('xl/'): tgt='xl/'+tgt
        rows=[]
        for row in ET.fromstring(z.read(tgt)).iter('{%s}row'%NS['m']):
            cells={}
            for c in row.findall('m:c',NS):
                v=c.find('m:v',NS); t=c.get('t')
                if t=='inlineStr':
                    isn=c.find('m:is',NS)
                    val="".join(x.text or "" for x in isn.iter('{%s}t'%NS['m'])) if isn is not None else ""
                elif v is None: continue
                elif t=='s': val=shared[int(v.text)]
                else: val=v.text
                cells[col_to_idx(c.get('r'))]=val
            rows.append(cells)
        width=max([max(r)+1 for r in rows if r] or [0])
        out[name]=[[ (r.get(i) or "") for i in range(width)] for r in rows]
    return out
