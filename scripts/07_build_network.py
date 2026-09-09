#!/usr/bin/env python3
"""
Clean class-network layout: 4,824 ChemOnt classes, laid out for READABILITY.

Lesson from the point-cloud attempts: 73.1M points on a 2M-pixel screen puts 550-1500
points on every occupied pixel, so it saturates to flat white no matter how it is
arranged. Here the compounds are carried by node SIZE (log-scaled) rather than drawn
individually, which is what makes the map legible.

Outputs:
  results/network_layout.json   node positions, radii, per-class point ranges
  results/network_preview.png   a render of exactly what the frontend will draw
"""
import json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont

ROOT = HERE.parent
RES = ROOT / "results"
G = json.loads((RES / "graph_data.json").read_text())
nodes = G["nodes"]
idx = {n["i"]: k for k, n in enumerate(nodes)}
NN = len(nodes)
cnt = np.array([n["dp"] for n in nodes], dtype=np.float64)
TOTAL = int(cnt.sum())

# node radius: log-scaled so a 26M class is ~6x a 100-compound class, not 500x
R = 2.6 + 7.4 * (np.log10(1 + cnt) / math.log10(1 + cnt.max())) ** 1.35
R[cnt == 0] = 1.9
dep = np.array([n["d"] for n in nodes])
R[dep == 1] *= 1.5

tree = [(idx[l["s"]], idx[l["t"]]) for l in G["tree"] if l["s"] in idx and l["t"] in idx]
mw = max(l["w"] for l in G["links"]) or 1
conf = [(idx[l["s"]], idx[l["t"]], l["w"] / mw) for l in G["links"]
        if l["s"] in idx and l["t"] in idx]
ti = np.array([e[0] for e in tree]); tj = np.array([e[1] for e in tree])
ci = np.array([e[0] for e in conf]); cj = np.array([e[1] for e in conf])
cw = np.array([e[2] for e in conf])

# ---------------------------------------------------------------------------
# RADIAL TREE + HIERARCHICAL EDGE BUNDLING
#
# A force layout on this data produces a hairball: 4,822 of the edges form a strict
# tree whose superclass hubs have hundreds of children, so repulsion pushes every node
# to the rim and the spokes fill the middle with white. Trees have an exact readable
# layout instead -- radius = depth, angle = subtree extent -- with zero edge crossings.
# The confusability links are then bundled along the tree (Holten), so cross-links read
# as flowing arcs rather than chords through the centre.
# ---------------------------------------------------------------------------
sc = [n["s"] if n["s"] is not None else -1 for n in nodes]
uniq = sorted(set(sc))

kids = defaultdict(list)
root_ids = []
for k, n in enumerate(nodes):
    p_ = n["p"]
    if p_ is None or p_ < 0 or p_ not in idx:
        root_ids.append(k)
    else:
        kids[idx[p_]].append(k)
print(f"tree roots: {[nodes[r]['n'] for r in root_ids]}")

# leaf angular weight: log-scaled compounds, floored so empty classes stay visible
leafw = np.log10(1 + cnt) + 0.65

wsum = np.zeros(NN)
order_post = []
def post(k):
    for c in kids[k]:
        post(c)
    order_post.append(k)
sys.setrecursionlimit(10000)
for r in root_ids:
    post(r)
for k in order_post:
    wsum[k] = leafw[k] if not kids[k] else sum(wsum[c] for c in kids[k])

TWO_PI = 2 * math.pi
GAPFRAC = 0.014                       # a wedge of blank sky between the two kingdoms
ang0 = np.zeros(NN); ang1 = np.zeros(NN)
total_w = sum(wsum[r] for r in root_ids)
cursor = -math.pi / 2
for r in sorted(root_ids, key=lambda z: -wsum[z]):
    span_r = TWO_PI * (1 - GAPFRAC * len(root_ids)) * wsum[r] / total_w
    ang0[r], ang1[r] = cursor, cursor + span_r
    cursor += span_r + TWO_PI * GAPFRAC

def assign(k):
    if not kids[k]:
        return
    lo, hi = ang0[k], ang1[k]
    tot = sum(wsum[c] for c in kids[k]) or 1.0
    # widest child first keeps big families contiguous and easy to name
    cur = lo
    for c in sorted(kids[k], key=lambda z: -wsum[z]):
        w_ = (hi - lo) * wsum[c] / tot
        ang0[c], ang1[c] = cur, cur + w_
        cur += w_
        assign(c)
for r in root_ids:
    assign(r)

MAXD = int(dep.max())
# evenly-spaced rings: the earlier geometric taper crushed depths 5-11 together and
# produced thousands of radial collisions
RING = np.linspace(0.0, 1.0, 12) ** 0.78
rad_r = RING[np.clip(dep, 0, len(RING) - 1)] * 1000.0
mid = (ang0 + ang1) / 2
P = np.stack([np.cos(mid) * rad_r, np.sin(mid) * rad_r], 1)

# node radius: log-scaled, and never wider than its own angular slot at that ring
R = 2.2 + 8.2 * (np.log10(1 + cnt) / math.log10(1 + cnt.max())) ** 1.30
R[cnt == 0] = 1.6
slot = (ang1 - ang0) * np.maximum(rad_r, 60.0) * 0.46
R = np.minimum(R, np.maximum(slot, 0.9))
ring_gap = np.diff(RING).min() * 1000.0
R = np.minimum(R, ring_gap * 0.40)          # never touch the next ring
R[dep <= 1] = np.maximum(R[dep <= 1], 7.0)

d = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
np.fill_diagonal(d, np.inf)
overlap = int((d < (R[:, None] + R[None, :])).sum() // 2)
print(f"overlapping node pairs: {overlap}")

# ---- bundle the confusability edges along the tree -------------------------
parent = np.full(NN, -1)
for k, n in enumerate(nodes):
    if n["p"] is not None and n["p"] >= 0 and n["p"] in idx:
        parent[k] = idx[n["p"]]
def chain(k):
    out = []
    while k != -1:
        out.append(k); k = parent[k]
    return out
def bundle_path(a, b):
    ca, cb = chain(a), chain(b)
    sa = set(ca)
    lca = next((x for x in cb if x in sa), None)
    if lca is None:
        return [a, b]
    up = ca[:ca.index(lca) + 1]
    dn = cb[:cb.index(lca)][::-1]
    return up + dn

BETA = 0.86
bundles = []
for (a, b, w) in sorted(conf, key=lambda e: -e[2])[:5200]:
    path = bundle_path(a, b)
    if len(path) < 2:
        continue
    pts = P[path]
    # straighten toward the chord by beta, the standard bundling relaxation
    t = np.linspace(0, 1, len(pts))[:, None]
    straight = pts[0] * (1 - t) + pts[-1] * t
    ctrl = BETA * pts + (1 - BETA) * straight
    bundles.append({"c": [[round(float(x), 1), round(float(y), 1)] for x, y in ctrl],
                    "w": float(w),
                    "a": nodes[a]["i"], "b": nodes[b]["i"]})
print(f"bundled {len(bundles):,} confusability arcs")

# radial dendrogram links: an arc along the parent's ring, then a radial spoke.
# Straight chords produced the dark wedge artefacts in the first preview.
treelinks = []
for a, b in tree:
    r0_, r1_ = rad_r[a], rad_r[b]
    t0_, t1_ = mid[a], mid[b]
    steps = max(2, min(14, int(abs(t1_ - t0_) * max(r0_, 1) / 26) + 2))
    arc = [[float(math.cos(t0_ + (t1_ - t0_) * u) * r0_),
            float(math.sin(t0_ + (t1_ - t0_) * u) * r0_)]
           for u in np.linspace(0, 1, steps)]
    arc.append([float(math.cos(t1_) * r1_), float(math.sin(t1_) * r1_)])
    sib = len(kids[a]) or 1
    treelinks.append({"p": [[round(x, 1), round(y, 1)] for x, y in arc],
                      "d": int(dep[b]), "s": sib,
                      "pa": nodes[a]["i"], "ch": nodes[b]["i"]})
print(f"radial dendrogram links: {len(treelinks):,}")

span = np.abs(P).max(0)
# ---- preview render ---------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

byId = {n["i"]: n for n in nodes}
scs = [s for s in uniq if s != -1]
org = [s for s in scs if byId[s]["k"] == 0]
ino = [s for s in scs if byId[s]["k"] == 1]
hue = {}
for i, s in enumerate(org): hue[s] = (348 + i / max(1, len(org) - 1) * 82) % 360
for i, s in enumerate(ino): hue[s] = 166 + i / max(1, len(ino) - 1) * 88
import colorsys
def col(k):
    s = sc[k]
    if s == -1: return (0.95, 0.65, 0.35) if nodes[k]["k"] == 0 else (0.31, 0.76, 0.69)
    return colorsys.hls_to_rgb(hue[s] / 360, 0.63, 0.62)

W = 26
fig = plt.figure(figsize=(W, W), facecolor="#080B12")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor("#080B12"); ax.axis("off")
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
for bd in bundles[:5200]:
    c = np.array(bd["c"])
    if len(c) < 3:
        continue
    verts = [c[0]]
    codes = [MPath.MOVETO]
    for i in range(1, len(c) - 1):
        verts += [c[i], (c[i] + c[i + 1]) / 2]
        codes += [MPath.CURVE3, MPath.CURVE3]
    verts += [c[-1], c[-1]]
    codes += [MPath.CURVE3, MPath.CURVE3]
    ax.add_patch(PathPatch(MPath(verts, codes), fill=False,
                           edgecolor=(0.47, 0.60, 0.85),
                           lw=0.30, alpha=0.085, zorder=1))
# A parent with 400 children stacks 400 near-parallel spokes into a solid white
# wedge. Scale each link's opacity by its sibling count so a dense fan sums to the
# same ink as a sparse one.
for tl in treelinks:
    c = np.array(tl["p"])
    fan = 1.0 / (1.0 + tl["s"] / 9.0)
    ax.plot(c[:, 0], c[:, 1], "-", color=(0.66, 0.71, 0.82),
            lw=0.5 if tl["d"] <= 3 else 0.28,
            alpha=min(0.34, (0.34 if tl["d"] <= 3 else 0.20) * fan + 0.012), zorder=2)
ax.scatter(P[:, 0], P[:, 1], s=(R * 1.55) ** 2, c=[col(k) for k in range(NN)],
           linewidths=0, alpha=0.95, zorder=3)
taken = []
placed = 0
for k in np.argsort(-cnt):
    if placed >= 34:
        break
    a_ = mid[k] % TWO_PI
    if any(min(abs(a_ - t), TWO_PI - abs(a_ - t)) < 0.085 for t in taken):
        continue
    taken.append(a_)
    rr = rad_r[k] + R[k] + 26
    deg = math.degrees(a_)
    flip = 90 < (deg % 360) < 270
    ax.text(math.cos(a_) * rr, math.sin(a_) * rr, nodes[k]["n"][:30],
            ha="right" if flip else "left", va="center",
            rotation=deg + 180 if flip else deg, rotation_mode="anchor",
            fontsize=7.0, color="#E8EDF7", zorder=5)
    placed += 1
m = np.array([span.max(), span.max()]) * 1.20
ax.set_xlim(-m[0], m[0]); ax.set_ylim(-m[1], m[1])
out = RES / "network_preview.png"
fig.savefig(out, dpi=76, facecolor="#080B12")
print(f"preview -> {out}")

# ---- normalise + emit -------------------------------------------------------
half = span.max() * 1.06
Pn = P / half
Rn = R / half
# bundles and dendrogram links are built in raw layout units; put them in the SAME
# normalised space as the nodes, or the frontend draws them 1000x off.
for bd in bundles:
    bd["c"] = [[round(x / half, 5), round(y / half, 5)] for x, y in bd["c"]]
for tl in treelinks:
    tl["p"] = [[round(x / half, 5), round(y / half, 5)] for x, y in tl["p"]]
order = np.argsort(-cnt)
off, ranges = 0, {}
for k in order:
    c = int(cnt[k])
    ranges[nodes[k]["i"]] = [off, c]
    off += c
assert off == TOTAL

layout = {
    "meta": {**G["meta"], "points": TOTAL, "aspect": float(span[0] / span[1]),
             "bytes_per_point": 4, "bin": "points.bin"},
    "bundles": bundles,
    "treelinks": treelinks,
    "nodes": [{**n, "px": round(float(Pn[idx[n["i"]], 0]), 5),
               "py": round(float(Pn[idx[n["i"]], 1]), 5),
               "nr": round(float(Rn[idx[n["i"]]]), 5),
               "off": ranges[n["i"]][0], "cnt": ranges[n["i"]][1]} for n in nodes],
    "tree": G["tree"], "links": G["links"],
}
(RES / "network_layout.json").write_text(json.dumps(layout, separators=(",", ":")))
print(f"wrote network_layout.json  aspect {span[0]/span[1]:.2f}:1")
