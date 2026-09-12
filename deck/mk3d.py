# -*- coding: utf-8 -*-
"""Backbone traces rendered straight from the project's own PDB files, as theme-safe SVG."""
import os
import json
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def records(path, chain):
    """dedupe: chain B of the SALL4 file has every record written twice"""
    seen, out = set(), []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")): continue
        if l[21] != chain: continue
        key = (l[17:20], l[21], l[22:27], l[12:16])
        if key in seen: continue
        seen.add(key)
        out.append((l[:6].strip(), int(l[22:26]), l[17:20].strip(), l[12:16].strip(),
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), float(l[60:66])))
    return out

def frame(P):
    C = P - P.mean(0)
    _, _, V = np.linalg.svd(C, full_matrices=False)
    R = V
    if np.linalg.det(R) < 0: R[2] *= -1
    return R

def catmull(p, n=4):
    p = np.asarray(p, float)
    if len(p) < 2: return p
    if n <= 1: return p
    ext = np.vstack([p[0], p, p[-1], p[-1]]); out = []
    for i in range(len(p) - 1):
        p0, p1, p2, p3 = ext[i], ext[i+1], ext[i+2], ext[i+3]
        for k in range(n):
            t = k/n; t2, t3 = t*t, t*t*t
            out.append(.5*((2*p1) + (-p0+p2)*t + (2*p0-5*p1+4*p2-p3)*t2 + (-p0+3*p1-3*p2+p3)*t3))
    out.append(p[-1])
    return np.array(out)

def emit(strands, spheres, W, H, pad, wide, nb=8):
    # strands: (pts, colourvar, opacity_scale, width_scale)
    """painter's algorithm on depth buckets; consecutive same-bucket segments merge into one polyline"""
    pts_all = np.vstack([s[0] for s in strands] + ([np.array([s[0] for s in spheres])] if spheres else []))
    lo, hi = pts_all[:, :2].min(0), pts_all[:, :2].max(0)
    sc = (min(W - 2*pad, H - 2*pad)) / (hi - lo).max()
    cx, cy = (lo + hi) / 2
    X = lambda p: (W/2 + (p[0]-cx)*sc, H/2 - (p[1]-cy)*sc)
    zl, zh = pts_all[:, 2].min(), pts_all[:, 2].max()
    bucket = lambda z: min(nb-1, max(0, int((z - zl) / (zh - zl + 1e-9) * nb)))
    prims = []                                    # (bucket, svg)
    for pts, col, osc, wsc in strands:
        run, rb = [X(pts[0])], bucket((pts[0][2] + pts[1][2]) / 2)
        for i in range(len(pts) - 1):
            b = bucket((pts[i][2] + pts[i+1][2]) / 2)
            if b != rb:
                if len(run) > 1: prims.append((rb, col, osc, wsc, run))
                run, rb = [run[-1]], b
            run.append(X(pts[i+1]))
        if len(run) > 1: prims.append((rb, col, osc, wsc, run))
    out = []
    for b in range(nb):
        d = (b + .5) / nb
        w = wide[0] + (wide[1] - wide[0]) * d
        for pb, col, osc, wsc, run in prims:
            if pb != b: continue
            pl = " ".join("%.1f,%.1f" % q for q in run)
            out.append('<polyline points="%s" fill="none" stroke="var(%s)" stroke-width="%.1f" '
                       'stroke-opacity="%.2f" stroke-linecap="round" stroke-linejoin="round"/>'
                       % (pl, col, w*wsc, min(1, (.34 + .66*d) * osc)))
        for p, r, fill in spheres:
            if bucket(p[2]) != b: continue
            out.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="var(%s)"/>' % (*X(p), r*sc, fill))
    return "".join(out), X, sc

# ---------------------------------------------------------- SALL4 zinc finger (C2H2)
zfp = f"{ROOT}/CRBN-Thalidomide-SALL4_(Ryan).pdb"
rec = records(zfp, "B")
ca = [(r[1], r[4]) for r in rec if r[0] == "ATOM" and r[3] == "CA"]
zn = [r[4] for r in rec if r[0] == "HETATM" and r[2] == "ZN"]
coord = []
for r in rec:
    if r[0] != "ATOM": continue
    if (r[2] == "CYS" and r[3] == "SG") or (r[2] == "HIS" and r[3] in ("NE2", "ND1")):
        if any(np.linalg.norm(r[4] - z) < 3.0 for z in zn): coord.append((r[1], r[2], r[4]))
print("ZF: %d residues %d-%d, %d Zn, clamped by %s"
      % (len(ca), ca[0][0], ca[-1][0], len(zn), [f"{c[1]}{c[0]}" for c in coord]))

P = np.array([c[1] for c in ca]); R = frame(P); mu = P.mean(0)
pj = lambda X_: (np.atleast_2d(X_) - mu) @ R.T
nums = np.array([c[0] for c in ca]); Q = pj(P)
# ββα: CA(i)-CA(i+4) is ~12 A through 420 and ~6.3 A from 421 on, so the helix starts at 421.
# The two Cys sit in the hairpin, the two His in the helix - that is what makes it C2H2.
HELIX = 421
strands = [(catmull(Q[nums <= HELIX], 5), "--accent",  1.0, 1.0),   # beta hairpin: what CRBN reads
           (catmull(Q[nums >= HELIX], 5), "--s-ctrl", 1.0, 1.0)]    # alpha helix
for rn, res, xyz in coord:
    b = [c[1] for c in ca if c[0] == rn]
    col = "--accent" if res == "CYS" else "--s-ctrl"                 # colour each side chain like its element
    if b: strands.append((pj(np.vstack([b[0], xyz])), col, .95, .70))
print("  hairpin %d-%d, helix %d-%d" % (nums.min(), HELIX-1, HELIX, nums.max()))
zf_body, Xf, scf = emit(strands, [(pj(z)[0], 1.5, "--s-cross") for z in zn], 214, 152, 15, (2.3, 5.2))
ZNXY = Xf(pj(zn[0])[0])
zf = zf_body

# ---------------------------------------------------------- PPIL4 (AlphaFold model)
rec4 = records(f"{ROOT}/PPIL4_alphafold_(Ryan).pdb", "A")
ca4 = [(r[1], r[4], r[5]) for r in rec4 if r[0] == "ATOM" and r[3] == "CA"]
P4 = np.array([c[1] for c in ca4]); B4 = np.array([c[2] for c in ca4])
core = P4[B4 >= 70]
R4 = frame(core); mu4 = core.mean(0)
Q = (P4 - mu4) @ R4.T
# two folded modules, 67 A apart by 2-means on the confident CAs:
# the N-terminal cyclophilin/PPIase domain and an RRM. Everything else is low-confidence.
N4 = np.array([c[0] for c in ca4])
def label(i):
    if B4[i] < 70: return "dis"
    return "cyp" if N4[i] <= 187 else ("rrm" if N4[i] >= 204 else "dis")
STYLE = {"cyp": ("--s-ctrl", 1.0, 1.0), "rrm": ("--s-dom2", 1.0, 1.0), "dis": ("--muted", .34, .58)}
runs, cur, curl = [], [0], label(0)
for i in range(1, len(Q)):
    l = label(i)
    if l != curl: cur.append(i); runs.append((cur, curl)); cur, curl = [i], l
    else: cur.append(i)
runs.append((cur, curl))
s4 = [(catmull(Q[idx], 4 if lab != "dis" else 1), *STYLE[lab]) for idx, lab in runs if len(idx) > 1]
for lab in ("cyp", "rrm", "dis"):
    k = [N4[i] for i in range(len(N4)) if label(i) == lab]
    print("  %-4s %d residues (%d-%d)" % (lab, len(k), min(k), max(k)))
p4_body, _, _ = emit(s4, [], 234, 152, 11, (2.0, 4.6))
p4 = p4_body
print("PPIL4: %d residues, %d confident (%.0f%%), %d disordered"
      % (len(ca4), (B4 >= 70).sum(), 100*(B4 >= 70).mean(), (B4 < 70).sum()))

figs = json.load(open("figs.json")); figs["zf3d"] = zf; figs["ppil43d"] = p4
figs["zf3d_zn"] = "%.1f,%.1f" % ZNXY
json.dump(figs, open("figs.json", "w"))
print("Zn anchor in the 214x152 frame:", "%.1f,%.1f" % ZNXY)
print("sizes: zf %.1f KB, ppil4 %.1f KB" % (len(zf)/1024, len(p4)/1024))
