import os
import json, re, sys
SP = os.path.dirname(os.path.abspath(__file__))
html=open(f"{SP}/tpl.html").read()
figs=json.load(open(f"{SP}/figs.json"))
pngs=json.load(open(f"{SP}/png.json"))
pts=json.load(open(f"{SP}/pts.json"))

data=[{"n":p["n"],"v":round(p["v"],3),"role":p["role"],"ha":p["ha"],"dq":p["dq"]} for p in pts]
html=html.replace("__DATA__", json.dumps(data, separators=(",",":")))

miss=[]
for m in sorted(set(re.findall(r"__FIG_([A-Za-z0-9_]+)__", html))):
    if m in figs: html=html.replace(f"__FIG_{m}__", figs[m])
    else: miss.append("FIG:"+m)
for m in sorted(set(re.findall(r"__PNG_([A-Za-z0-9_]+)__", html))):
    if m in pngs: html=html.replace(f"__PNG_{m}__", pngs[m])
    else: miss.append("PNG:"+m)

# fix the duplicate style attribute written by hand
html=html.replace('<div class="stack" style="gap:9px" style="padding-top:26px">',
                  '<div class="stack" style="gap:9px; padding-top:22px">')

left=re.findall(r"__[A-Z][A-Za-z0-9_]*__", html)
open(f"{SP}/deck.html","w").write(html)
print("missing:",miss or "none")
print("unresolved placeholders:", set(left) or "none")
print("slides:", html.count('<section class="slide"'))
print("size: %.1f KB" % (len(html.encode())/1024))
