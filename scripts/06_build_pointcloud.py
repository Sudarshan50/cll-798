#!/usr/bin/env python3
"""
Generate one screen position for every one of the 73,105,281 compounds.

Honest description of what a point is: each point is one real compound, placed inside
the region of the ChemOnt class it was actually assigned to. Counts per class are exact,
taken from the full 73.1M-row scan -- the totals reconcile to 73,105,281 with no
sampling. Within a class, points are laid out on a deterministic Vogel (sunflower)
spiral, because we have no 2D structural embedding: position within a class carries no
meaning, position between classes carries the taxonomy and its confusability overlay.

Outputs:
  results/layout.json    4,824 class positions + metadata (drives the interactive layer)
  public/points.bin      int16 x,y per compound, 4 bytes each  (~292 MB)
"""
import json, math, struct, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont

ROOT = HERE.parent
RES = ROOT / "results"
OUTDIR = ROOT / "frontend" / "public"
OUTDIR.mkdir(parents=True, exist_ok=True)
tax = chemont.from_dictionary_tsv(ROOT / "data" / "chemont_dictionary.tsv")

G = json.loads((RES / "graph_data.json").read_text())
nodes = [n for n in G["nodes"]]
idx = {n["i"]: k for k, n in enumerate(nodes)}
NN = len(nodes)

# ---- layout: superclass islands packed into a wide strip ---------------------
# A single global force layout piles every class into one dense disc, which is
# unreadable. ChemOnt already carries the grouping we want: 31 superclasses. Each
# becomes an island laid out on its own, then islands are shelf-packed into a wide
# landscape band so the map is panned rather than squinted at.

cnt_all = np.array([n["dp"] for n in nodes], dtype=np.float64)

groups = defaultdict(list)
for k, n in enumerate(nodes):
    g = n["s"] if n["s"] is not None else ("KINGDOM_%d" % n["k"])
    groups[g].append(k)

tree_pairs = [(idx[l["s"]], idx[l["t"]]) for l in G["tree"] if l["s"] in idx and l["t"] in idx]
conf_pairs = [(idx[l["s"]], idx[l["t"]], l["w"]) for l in G["links"] if l["s"] in idx and l["t"] in idx]
mw = max((w for _, _, w in conf_pairs), default=1) or 1

rng = np.random.default_rng(7)


def fr_layout(m, edges, iters=260, seed=0):
    """Damped Fruchterman-Reingold on a small subgraph; returns unit-disc coords."""
    if m == 1:
        return np.zeros((1, 2))
    r = np.random.default_rng(seed)
    ang = r.uniform(0, 2 * math.pi, m)
    rad = np.sqrt(r.uniform(0, 1, m))
    P = np.stack([np.cos(ang) * rad, np.sin(ang) * rad], 1) * (1.4 * math.sqrt(m))
    if edges:
        ei = np.array([e[0] for e in edges]); ej = np.array([e[1] for e in edges])
        ew = np.array([e[2] for e in edges], dtype=np.float64)
    K = 1.9
    temp = 1.6 * math.sqrt(m)
    for t in range(iters):
        d = P[:, None, :] - P[None, :, :]
        dist = np.sqrt((d ** 2).sum(-1))
        np.fill_diagonal(dist, np.inf)
        dist = np.maximum(dist, 0.02)
        disp = (d / dist[:, :, None] * (K * K / dist)[:, :, None]).sum(1)
        if edges:
            dv = P[ej] - P[ei]
            L = np.maximum(np.sqrt((dv ** 2).sum(1)), 1e-4)
            f = dv / L[:, None] * (ew * L * L / K)[:, None]
            att = np.zeros_like(P)
            np.add.at(att, ei, f); np.add.at(att, ej, -f)
            disp -= att
        dl = np.maximum(np.sqrt((disp ** 2).sum(1)), 1e-9)
        P += disp / dl[:, None] * np.minimum(dl, temp)[:, None]
        temp *= 0.982
        if not np.isfinite(P).all():
            raise SystemExit("island layout diverged")
    P -= P.mean(0)
    m2 = np.sqrt((P ** 2).sum(1)).max()
    return P / (m2 if m2 > 1e-9 else 1.0)


islands = []
print(f"laying out {len(groups)} superclass islands ...")
for gi, (g, members) in enumerate(sorted(groups.items(), key=lambda kv: -sum(cnt_all[k] for kv2 in [kv[1]] for k in kv2))):
    loc = {k: i for i, k in enumerate(members)}
    ed = [(loc[a_], loc[b_], 1.0) for a_, b_ in tree_pairs if a_ in loc and b_ in loc]
    ed += [(loc[a_], loc[b_], 0.10 * w / mw) for a_, b_, w in conf_pairs if a_ in loc and b_ in loc]
    L = fr_layout(len(members), ed, seed=gi + 1)
    pop = float(cnt_all[members].sum())
    nm = nodes[members[0]]["sc"] if not str(g).startswith("KINGDOM") else nodes[members[0]]["n"]
    islands.append({"key": g, "members": members, "local": L, "pop": pop, "name": nm,
                    "k": nodes[members[0]]["k"]})
    print(f"   {nm[:46]:46} {len(members):5} classes  {int(pop):>12,} compounds")

# ---- squarified treemap ------------------------------------------------------
# Circle packing left 96% of the canvas empty and crushed 73.1M points into 3.8% of it,
# so every occupied pixel saturated (median 1,494 points/px). Rectangles tile without
# gaps, which is the only way to give the cloud enough area to form a readable gradient.
#
# Cell AREA is proportional to sqrt(compounds), not compounds: the top two superclasses
# hold 68% of all compounds, so true-proportional area would reduce the other 31 to
# slivers. Every cell is labelled with its true count, so the number carries the truth
# and the area carries a compressed, declared encoding.

STRIP_W, STRIP_H = 5.6, 2.0                     # world units -> 2.8:1
AREA = STRIP_W * STRIP_H
wt = np.array([max(i["pop"], 1.0) ** 0.5 for i in islands], dtype=np.float64)
wt = np.maximum(wt, wt.sum() * 0.006)           # floor so nothing vanishes
wt = wt / wt.sum() * AREA
for isl, w_ in zip(islands, wt):
    isl["area"] = float(w_)

def squarify(items, x, y, w, h, out):
    """Standard squarified treemap; appends (item, x, y, w, h)."""
    if not items:
        return
    if len(items) == 1:
        out.append((items[0], x, y, w, h)); return
    total = sum(i["area"] for i in items)
    if w >= h:
        best, cut = None, 1
        for k in range(1, len(items)):
            a_ = sum(i["area"] for i in items[:k])
            cw = w * a_ / total
            worst = max(max(cw / (a_ / cw * (i["area"] / a_) * 0 + max(1e-9, i["area"] / cw)),
                            max(1e-9, i["area"] / cw) / cw) for i in items[:k])
            if best is None or worst < best:
                best, cut = worst, k
        a_ = sum(i["area"] for i in items[:cut])
        cw = w * a_ / total
        yy = y
        for i in items[:cut]:
            ih = h * i["area"] / a_
            out.append((i, x, yy, cw, ih)); yy += ih
        squarify(items[cut:], x + cw, y, w - cw, h, out)
    else:
        best, cut = None, 1
        for k in range(1, len(items)):
            a_ = sum(i["area"] for i in items[:k])
            ch = h * a_ / total
            worst = max(max(ch / max(1e-9, i["area"] / ch), max(1e-9, i["area"] / ch) / ch)
                        for i in items[:k])
            if best is None or worst < best:
                best, cut = worst, k
        a_ = sum(i["area"] for i in items[:cut])
        ch = h * a_ / total
        xx = x
        for i in items[:cut]:
            iw = w * i["area"] / a_
            out.append((i, xx, y, iw, ch)); xx += iw
        squarify(items[cut:], x, y + ch, w, h - ch, out)

cells = []
squarify(sorted(islands, key=lambda i: -i["area"]), 0.0, 0.0, STRIP_W, STRIP_H, cells)
assert len(cells) == len(islands), f"treemap lost cells: {len(cells)} != {len(islands)}"

PAD = 0.035
P = np.zeros((NN, 2))
rad = np.zeros(NN)
for isl, cx0, cy0, cw, ch in cells:
    iw, ih = max(cw - 2 * PAD, 0.05), max(ch - 2 * PAD, 0.05)
    cx, cy = cx0 + cw / 2, cy0 + ch / 2
    isl.update(x=cx0 + PAD, y=cy0 + PAD, w=iw, h=ih, cx=cx, cy=cy)
    mem = isl["members"]
    L = isl["local"]                                  # unit disc
    P[mem] = np.stack([L[:, 0] * iw * 0.46, L[:, 1] * ih * 0.46], 1) + np.array([cx, cy])
    tot = max(cnt_all[mem].sum(), 1.0)
    share = cnt_all[mem] / tot
    # fill ~62% of the cell with class discs so the cloud has room to form a gradient
    rad[mem] = np.sqrt(np.maximum(share, 1e-9) * 0.62 * iw * ih / math.pi)
    rad[mem] = np.minimum(rad[mem], min(iw, ih) * 0.44)
    rad[mem] = np.maximum(rad[mem], 0.004)

fill = sum(math.pi * rad[i["members"]].__pow__(2).sum() for i in islands) / AREA
print(f"\ntreemap: {len(cells)} cells over {STRIP_W}x{STRIP_H}  "
      f"disc coverage {fill*100:.1f}% of canvas")
assert fill > 0.25, f"coverage still too low: {fill*100:.1f}%"

assert np.isfinite(P).all(), "layout contains NaN/inf"
# Centre and scale ONCE, and put the islands through the identical transform -- they
# were previously scaled but not shifted, so every halo and label sat off its contents.
shift_x, shift_y = STRIP_W / 2, STRIP_H / 2
P[:, 0] -= shift_x
P[:, 1] -= shift_y
half_h = STRIP_H / 2
P /= half_h
rad /= half_h
for isl in islands:
    for k_ in ("cx", "x"):
        isl[k_] = (isl[k_] - shift_x) / half_h
    for k_ in ("cy", "y"):
        isl[k_] = (isl[k_] - shift_y) / half_h
    isl["w"] /= half_h
    isl["h"] /= half_h
    isl["R"] = min(isl["w"], isl["h"]) / 2   # w/h are already normalised; do not divide again
    isl["R"] = isl["R"] / half_h
SHIFT_X = (P[:, 0].max() + P[:, 0].min()) / 2
ASPECT = np.abs(P[:, 0]).max()
print(f"strip extent: x +/-{ASPECT:.2f}  y +/-{np.abs(P[:,1]).max():.2f}  "
      f"-> aspect {ASPECT:.2f}:1")
assert ASPECT > 2.4, f"strip is not wide enough: {ASPECT:.2f}:1"

# (per-class disc radii are set by the treemap above; the old circle-packing
#  radius block used to run here and silently clobbered them back to the floor)

# every class disc must sit inside its treemap cell
for isl in islands:
    mem = isl["members"]
    ox = np.abs(P[mem, 0] - isl["cx"]) + rad[mem] - isl["w"] / 2
    oy = np.abs(P[mem, 1] - isl["cy"]) + rad[mem] - isl["h"] / 2
    if ox.max() > 1e-9 or oy.max() > 1e-9:
        sh = min(1.0, min(isl["w"] / 2 / (np.abs(P[mem, 0] - isl["cx"]) + rad[mem]).max(),
                          isl["h"] / 2 / (np.abs(P[mem, 1] - isl["cy"]) + rad[mem]).max()))
        P[mem] = (P[mem] - np.array([isl["cx"], isl["cy"]])) * sh + np.array([isl["cx"], isl["cy"]])
        rad[mem] *= sh

# ---- exact per-class compound counts ---------------------------------------
cnt = np.array([n["dp"] for n in nodes], dtype=np.int64)
TOTAL = int(cnt.sum())
print(f"compounds to place: {TOTAL:,}")

# ---- write the point cloud --------------------------------------------------
# int16 must span the FULL extent, not just y. The strip runs wider than it is tall,
# so a fixed 32000 divisor silently clipped x to the int16 limit and flattened the map.
EXT = float(max(np.abs(P).max(), (np.abs(P[:, 0]) + rad).max(), (np.abs(P[:, 1]) + rad).max()))
S = 32600.0 / EXT
print(f"int16 scale: {S:.1f} counts per world unit (extent {EXT:.3f})")
out = OUTDIR / "points.bin"
GA = math.pi * (3 - math.sqrt(5))                  # golden angle
CH = 2_000_000
written = 0
with out.open("wb") as fh:
    for k in range(NN):
        n = int(cnt[k])
        if n == 0:
            continue
        cx, cy, rr = P[k, 0], P[k, 1], rad[k]
        done = 0
        while done < n:
            m = min(CH, n - done)
            j = np.arange(done, done + m, dtype=np.float64)
            th = j * GA
            rq = rr * np.sqrt((j + 0.5) / n)
            xf = (cx + np.cos(th) * rq) * S
            yf = (cy + np.sin(th) * rq) * S
            if np.abs(xf).max() > 32767 or np.abs(yf).max() > 32767:
                raise SystemExit("int16 overflow: scale is wrong")
            x = xf.astype(np.int16); y = yf.astype(np.int16)
            np.stack([x, y], 1).tofile(fh)
            done += m
            written += m
        if k % 600 == 0:
            print(f"   class {k:4d}/{NN}  {written:>12,} points")
print(f"\nwrote {out}  {out.stat().st_size/1e6:.1f} MB  ({written:,} points)")
assert written == TOTAL, f"{written} != {TOTAL}"

# verify what actually landed on disk, rather than trusting the writer
chk = np.fromfile(out, dtype=np.int16, count=4_000_000).reshape(-1, 2)
frac_origin = float((np.abs(chk).sum(1) == 0).mean())
print(f"sample check: {frac_origin*100:.3f}% of the first 2M points sit at the origin, "
      f"x range {chk[:,0].min()}..{chk[:,0].max()}")
assert frac_origin < 0.02, "points collapsed to the origin -- layout is broken"
full = np.fromfile(out, dtype=np.int16).reshape(-1, 2)
assert np.abs(full).max() < 32767, "points hit the int16 limit -- the map is clipped"
print(f"full-file check: x {full[:,0].min()}..{full[:,0].max()}  "
      f"y {full[:,1].min()}..{full[:,1].max()}  (no clipping)")
del full

layout = {
    "meta": {**G["meta"], "points": TOTAL, "scale": S, "aspect": float(ASPECT),
             "bytes_per_point": 4, "bin": "points.bin"},
    "islands": [{"name": i["name"], "k": i["k"], "cx": round(i["cx"], 5),
                 "cy": round(i["cy"], 5), "x": round(i["x"], 5), "y": round(i["y"], 5),
                 "w": round(i["w"], 5), "h": round(i["h"], 5), "r": round(i["R"], 5),
                 "pop": int(i["pop"]), "classes": len(i["members"])}
                for i in sorted(islands, key=lambda z: -z["pop"])],
    "nodes": [{**n, "px": round(float(P[idx[n["i"]], 0]), 5),
               "py": round(float(P[idx[n["i"]], 1]), 5),
               "pr": round(float(rad[idx[n["i"]]]), 5)} for n in nodes],
    "tree": G["tree"], "links": G["links"],
}
(RES / "layout.json").write_text(json.dumps(layout, separators=(",", ":")))
(OUTDIR / "layout.json").write_text(json.dumps(layout, separators=(",", ":")))
print(f"wrote layout.json ({(RES/'layout.json').stat().st_size/1e6:.2f} MB)")
