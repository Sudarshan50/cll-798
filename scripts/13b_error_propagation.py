#!/usr/bin/env python3
"""
Phase 2, extra (project abstract, objective d) -- Monte Carlo error propagation.

If ClassyFire had picked a competitor instead of the direct parent it actually chose, how
far up the published 5-slot label would the mistake reach? Compounds are sampled in
proportion to their direct parent's population, the direct parent is replaced by a
competitor a ~ p(a|d) restricted to an overlay, and the shallowest of kingdom / superclass
/ class / subclass whose padded slot changes is recorded ('none' = only the direct-parent
slot moves). Two nulls draw the wrong class instead uniformly over the other direct-parent
classes, or in proportion to their size.

Overlays: the tree-free B-only overlay is the informative one. On the primary A(4) ∩ B the
answer is fixed by construction (LCA depth >= 4 means the first four slots agree), so it is
run only as a check. The unfiltered S1 overlay (every alternative parent treated as a
competitor) is the upper bound the descriptive alternative parents would imply.

Uncertainty: 200 replicates, each re-drawing the (dp, alt) counts and class sizes from a
Poisson bootstrap and then 20,000 compounds; 95 % percentile intervals.
"""
import json, sys, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase2_common as C  # noqa: E402

REPS, DRAWS = 200, 20000
t0 = time.time()
tax = C.taxonomy()
A, depth = C.anc_matrix(tax, C.K), C.depth_array(tax, C.K)
meta = json.loads((C.OUT / "overlay_meta.json").read_text())
TAU, NMIN = meta["separation"]["tau"], meta["separation"]["nmin"]
d, a, n = C.load_s1()
pmi, dtot, _ = C.pmi_directed(d, a, n)
U = C.Undirected(d, a)
nu, nv, ps = U.collapse(n, pmi)
lca = C.lca_depth(A, U.u, U.v)
size = C.compounds_as_dp()
dps = np.flatnonzero(dtot > 0)
level = C.first_changed_level(A, depth, d, a)          # per directed row
NAMES = list(C.LEVELS) + ["none"]
overlays = {
    "B_only": C.directed_mask(U, C.criterion_b(ps, nu + nv, TAU, NMIN)),
    "primary_A_and_B": C.directed_mask(U, (lca >= C.THETA0) & C.criterion_b(ps, nu + nv, TAU, NMIN)),
    "unfiltered_S1": np.ones(len(d), bool),
}


def sampler(mask, nn):
    """Row sampler for a ~ p(a|d) on the masked rows: one searchsorted over
    rank(d) + within-d cumulative probability."""
    rows = np.flatnonzero(mask & (nn > 0))
    rows = rows[np.argsort(d[rows], kind="stable")]
    dd, w = d[rows], nn[rows]
    tot = np.bincount(dd, w, C.K)
    cp = np.zeros(len(rows))
    start = np.r_[0, np.flatnonzero(np.diff(dd)) + 1]
    cs = np.cumsum(w)
    base = np.repeat(cs[start] - w[start], np.diff(np.r_[start, len(rows)]))
    cp = (cs - base) / tot[dd]
    rank = np.full(C.K, -1); rank[dd[start]] = np.arange(len(start))
    key = rank[dd] + cp
    return rows, key, rank


def profile(levels):
    return np.bincount(levels, minlength=5) / len(levels)


rng = np.random.default_rng(C.SEED)
res = {k: [] for k in list(overlays) + ["uniform", "size"]}
cover = {k: [] for k in overlays}
expected = {}
for name, m in overlays.items():                      # exact expectation, no Monte Carlo
    rows = np.flatnonzero(m)
    has = np.bincount(d[rows], minlength=C.K) > 0
    pd_ = size * has / (size * has).sum()
    pl = np.zeros((C.K, 5))
    np.add.at(pl, (d[rows], level[rows]), n[rows])
    pl = pl / np.maximum(pl.sum(1, keepdims=True), 1)
    expected[name] = (pd_ @ pl).round(5).tolist()

for r in range(REPS):
    nn = rng.poisson(n).astype(float)
    sz = rng.poisson(size).astype(float)
    for name, m in overlays.items():
        rows, key, rank = sampler(m, nn)
        has = np.zeros(C.K, bool); has[d[rows]] = True
        w = sz * has
        cover[name].append(float(w.sum() / sz[dps].sum()))
        dd = rng.choice(C.K, DRAWS, p=w / w.sum())
        pick = rows[np.searchsorted(key, rank[dd] + rng.random(DRAWS) * (1 - 1e-12))]
        assert (d[pick] == dd).all()
        res[name].append(profile(level[pick]))
        if name == "B_only":                           # nulls on the same compound draw
            j = rng.integers(0, len(dps) - 1, DRAWS)       # skip the true class's slot
            alt_u = dps[j + (j >= np.searchsorted(dps, dd))]
            assert (alt_u != dd).all()
            res["uniform"].append(profile(C.first_changed_level(A, depth, dd, alt_u)))
            ps_ = sz[dps] / sz[dps].sum()
            alt_s = dps[rng.choice(len(dps), DRAWS, p=ps_)]
            for _ in range(20):
                bad = alt_s == dd
                if not bad.any():
                    break
                alt_s[bad] = dps[rng.choice(len(dps), bad.sum(), p=ps_)]
            res["size"].append(profile(C.first_changed_level(A, depth, dd, alt_s)))
print(f"{REPS} replicates x {DRAWS} draws: {time.time()-t0:.0f}s")


def summ(x):
    x = np.array(x)
    return {"mean": x.mean(0).round(5).tolist(), "ci95_lo": np.percentile(x, 2.5, 0).round(5).tolist(),
            "ci95_hi": np.percentile(x, 97.5, 0).round(5).tolist()}


out = {"levels": NAMES, "replicates": REPS, "draws_per_replicate": DRAWS,
       "config": {"theta": C.THETA0, "tau": TAU, "nmin": NMIN},
       "profiles": {k: summ(v) for k, v in res.items()},
       "expected_exact": expected,
       "coverage_compound_share_with_competitor": {k: {"mean": float(np.mean(v)),
                                                       "ci95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]}
                                                   for k, v in cover.items()},
       "notes": {"primary_A_and_B": "fixed by construction: LCA depth >= 4 forces the four "
                                    "label slots to agree, so every swap lands in 'none'",
                 "nulls": "drawn for the same compounds as B_only (direct parents that have at "
                          "least one B-competitor), uniform over the other 3,603 direct-parent "
                          "classes, or proportional to compounds-as-direct-parent"}}
b, u_, s_ = (np.array(out["profiles"][k]["mean"]) for k in ("B_only", "uniform", "size"))
out["headline"] = {"B_only_superclass_or_above": float(b[:2].sum()),
                   "uniform_superclass_or_above": float(u_[:2].sum()),
                   "size_superclass_or_above": float(s_[:2].sum()),
                   "B_only_within_subclass_or_none": float(b[3:].sum()),
                   "uniform_within_subclass_or_none": float(u_[3:].sum())}
C.dump(out, "error_propagation.json")
print(json.dumps({k: out["profiles"][k]["mean"] for k in res}, indent=1))
print(json.dumps(out["headline"], indent=1), json.dumps(out["coverage_compound_share_with_competitor"]))
