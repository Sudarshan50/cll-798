#!/usr/bin/env python3
"""
Phase 2, step 3 -- is the descriptive/competing separation arbitrary? (plan §3.3)

Sweeps theta x tau x n_min (36 cells) of A(theta) ∩ B(tau, n_min) and reports, per cell,
the three prior criteria the plan fixes before any validation: agreement of the two
independent criteria (Jaccard, and overlap as a share of each side), stability under a
Poisson bootstrap of the (dp, alt) counts, and face validity on the three cases Phase 1
identified independently. Nothing here reads the expert-flagged classes; the H1 columns of
the surface are added afterwards by 12_ and 13_ and merged by 13c_.

The bootstrap resamples the 1.05M-row sufficient statistic, not the 73.1M compounds: each
n(d,a) ~ Poisson(n(d,a)), marginals and PMI recomputed, 200 replicates.
"""
import csv, json, sys, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase2_common as C  # noqa: E402

REPS = 200
t0 = time.time()
tax = C.taxonomy()
A = C.anc_matrix(tax, C.K)
d, a, n = C.load_s1()
pmi, dtot, amarg = C.pmi_directed(d, a, n)
U = C.Undirected(d, a)
n_uv, n_vu, psym = U.collapse(n, pmi)
w_raw = n_uv + n_vu
lca = C.lca_depth(A, U.u, U.v)
keys = U.u * C.K + U.v
meta = json.loads((C.OUT / "overlay_meta.json").read_text())


def idx(x, y):
    q = min(x, y) * C.K + max(x, y)
    i = int(np.searchsorted(keys, q))
    return i if i < len(keys) and keys[i] == q else None


by = {v: k for k, v in tax.name.items()}
FACE = {
    "aryl_thioethers": [idx(by["Aryl thioethers"], by["Alkylarylthioethers"])],
    "oxosteroids": [idx(by["Oxosteroids"], by["3-oxosteroids"])],
    # the axis of the dominant PATH shape and R6 pair: Depsipeptides vs Peptidomimetics /
    # 'Amino acids, peptides, and analogues'; none of these pairs occurs in S1 at all
    "depsipeptides_axis": [idx(by["Depsipeptides"], by["Peptidomimetics"]),
                           idx(by["Depsipeptides"], by["Amino acids, peptides, and analogues"]),
                           idx(by["Carboxylic acids and derivatives"], by["Peptidomimetics"])],
}
face_in_s1 = {k: [i for i in v if i is not None] for k, v in FACE.items()}

Amask = {th: lca >= th for th in C.THETAS}
Bmask = {(t, m): C.criterion_b(psym, w_raw, t, m) for t in C.TAUS for m in C.NMINS}
cells = [(th, t, m) for th in C.THETAS for t in C.TAUS for m in C.NMINS]
members = {c: np.flatnonzero(Amask[c[0]] & Bmask[c[1:]]) for c in cells}
assert len(members[(C.THETA0, meta["separation"]["tau"], meta["separation"]["nmin"])]) \
    == meta["separation"]["A_and_B"], "sweep disagrees with 10_build_overlay"

# ---- Poisson bootstrap -------------------------------------------------------
rng = np.random.default_rng(C.SEED)
kept = {c: np.zeros(len(members[c])) for c in cells}
for r in range(REPS):
    ns = rng.poisson(n).astype(np.float64)
    pm, _, _ = C.pmi_directed(d, a, ns)
    nu, nv, ps = U.collapse(ns, pm)
    wr = nu + nv
    for t in C.TAUS:
        for m in C.NMINS:
            b = (ps >= t) & (wr >= m)
            for th in C.THETAS:
                kept[(th, t, m)] += b[members[(th, t, m)]]
print(f"bootstrap {REPS} reps: {time.time()-t0:.0f}s")

rows = []
for c in cells:
    th, t, m = c
    Am, Bm, mem = Amask[th], Bmask[(t, m)], members[c]
    inter, union = len(mem), int((Am | Bm).sum())
    frac = kept[c] / REPS
    wp = np.maximum(0, psym[mem])
    face = {}
    for k, ids in face_in_s1.items():
        hit = [i for i in ids if Am[i] and Bm[i]]
        face[k] = int(bool(hit))
        face[k + "_pctile"] = (round(float((wp <= np.maximum(0, psym[hit[0]])).mean() * 100), 1)
                               if hit else None)
        face[k + "_pctile_raw"] = (round(float((w_raw[mem] <= w_raw[hit[0]]).mean() * 100), 1)
                                   if hit else None)
    nodes = int(len(np.unique(np.concatenate([U.u[mem], U.v[mem]])))) if inter else 0
    rows.append({"theta": th, "tau": t, "nmin": m, "A": int(Am.sum()), "B": int(Bm.sum()),
                 "edges": inter, "mass": int(w_raw[mem].sum()), "nodes": nodes,
                 "jaccard": round(inter / union, 4) if union else 0.0,
                 "share_of_A": round(inter / Am.sum(), 4) if Am.sum() else 0.0,
                 "share_of_B": round(inter / Bm.sum(), 4) if Bm.sum() else 0.0,
                 "enrichment": round(inter / (Am.sum() * Bm.sum() / U.m), 2) if Am.sum() and Bm.sum() else 0.0,
                 "stable95": round(float((frac >= 0.95).mean()), 4) if inter else None,
                 "mean_retention": round(float(frac.mean()), 4) if inter else None,
                 **face, "primary": int(c == (C.THETA0, meta["separation"]["tau"], meta["separation"]["nmin"]))})

with (C.OUT / "separation_sweep.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("" if v is None else v) for k, v in r.items()})

# face validity against the S3 tier as well: the cases are proved violations there
s3 = set()
with (C.OUT / "overlay_edges.tsv").open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["src"] == "3":
            s3.add((int(r["u"]), int(r["v"])))
face_s3 = {k: any((min(p), max(p)) in s3 for p in pairs) for k, pairs in {
    "aryl_thioethers": [(by["Aryl thioethers"], by["Alkylarylthioethers"])],
    "oxosteroids": [(by["Oxosteroids"], by["3-oxosteroids"])],
    "depsipeptides_axis": [(by["Depsipeptides"], by["Peptidomimetics"]),
                           (by["Depsipeptides"], by["Amino acids, peptides, and analogues"]),
                           (by["Carboxylic acids and derivatives"], by["Peptidomimetics"])]}.items()}
prim = next(r for r in rows if r["primary"])
summary = {"replicates": REPS, "primary_cell": prim,
           "face_cases_present_in_s1": {k: bool(v) for k, v in face_in_s1.items()},
           "face_cases_in_s3_tier": face_s3,
           "cells_stable95_ge_0.9": sum(1 for r in rows if (r["stable95"] or 0) >= 0.9),
           "cells_jaccard_ge_0.1": sum(1 for r in rows if r["jaccard"] >= 0.1),
           "cells_share_of_smaller_ge_0.1": sum(1 for r in rows if max(r["share_of_A"], r["share_of_B"]) >= 0.1),
           "runtime_s": round(time.time() - t0, 1)}
C.dump(summary, "separation_summary.json")
print(json.dumps(summary, indent=1))
