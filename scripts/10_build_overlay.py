#!/usr/bin/env python3
"""
Phase 2, step 2 -- the weighted confusability overlay on ChemOnt (plan §3-§4).

Supersedes the unfiltered prototype in 05_build_graph_data.py without touching it, so the
Phase-1 frontend payload stays byte-reproducible. Three evidence channels are kept apart
(src 1 = dp/alt co-assignment, 2 = alt/alt co-occurrence, 3 = proved rule violations) and
every S1 pair carries all three weights (w_raw, w_cond, w_pmi), its LCA depth and sibling
flag, so the separation in 11_ can be re-cut from this file alone.

The PMI threshold and support floor are not free: they are set here by the calibration
rule pre-registered in docs/phase2_prereg.md §3 (R5 pairs + PATH-shape pairs only; the
expert-flagged classes are never read by this script).
"""
import csv, json, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase2_common as C  # noqa: E402

tax = C.taxonomy()
A, depth, par = C.anc_matrix(tax, C.K), C.depth_array(tax, C.K), C.parent_array(tax, C.K)
sc = C.superclass_label(tax)

# ---- S1 --------------------------------------------------------------------
d, a, n = C.load_s1()
assert len(np.unique(d * C.K + a)) == len(d), "duplicate (dp, alt) rows"
pmi, dtot, amarg = C.pmi_directed(d, a, n)
U = C.Undirected(d, a)
n_uv, n_vu, psym = U.collapse(n, pmi)
w_raw = n_uv + n_vu
lca = C.lca_depth(A, U.u, U.v)
sib = par[U.u] == par[U.v]
with np.errstate(divide="ignore", invalid="ignore"):
    w_cond = np.fmax(np.where(n_uv > 0, n_uv / dtot[U.u], np.nan),
                     np.where(n_vu > 0, n_vu / dtot[U.v], np.nan))
w_pmi = np.maximum(0, psym)
keys = U.u * C.K + U.v


def pair_index(pairs):
    q = np.array([x * C.K + y for x, y in pairs], dtype=np.int64)
    i = np.clip(np.searchsorted(keys, q), 0, len(keys) - 1)
    return i[keys[i] == q]


# LCA-depth table over directed rows (reproduces plan §1A.5)
lca_dir = C.lca_depth(A, d, a)
bins = [("0 (root)", lca_dir == 0), ("1 (Kingdom)", lca_dir == 1), ("2 (SuperClass)", lca_dir == 2),
        ("3 (Class)", lca_dir == 3), ("4 (SubClass)", lca_dir == 4), ("5-11", lca_dir >= 5),
        ("siblings", par[d] == par[a])]
lca_table = [{"depth": k, "pairs": int(m.sum()), "mass": int(n[m].sum())} for k, m in bins]

# ---- calibration (prereg §3): B alone, 9 cells -------------------------------
r5_rows = C.r5_pairs(tax, d, a)
n_shapes, shape_pairs = C.path_shape_pairs(tax)
cal_r5 = set(U.inv[r5_rows].tolist())
cal_path = set(pair_index(shape_pairs).tolist())
cal = np.array(sorted(cal_r5 | cal_path))
cells = []
for tau in C.TAUS:
    for nmin in C.NMINS:
        B = C.criterion_b(psym, w_raw, tau, nmin)
        rec, ret = B[cal].mean(), B.mean()
        cells.append({"tau": tau, "nmin": nmin, "B_edges": int(B.sum()), "retain": float(ret),
                      "recall": float(rec), "cal_retained": int(B[cal].sum()),
                      "lift": float(rec / ret) if ret else 0.0})
ok = [c for c in cells if c["recall"] >= 0.5]
best = (max(ok, key=lambda c: (c["lift"], c["recall"], -c["tau"])) if ok
        else max(cells, key=lambda c: (c["recall"], -c["tau"])))
TAU, NMIN = best["tau"], best["nmin"]
calib = {
    "rule": "max lift = recall/retain over tau in {1,2,3} x n_min in {10,100,1000}, "
            "criterion B alone, subject to recall >= 0.5 (docs/phase2_prereg.md §3)",
    "r5_directed_pairs": int(len(r5_rows)),
    "r5_direct_parent_child": int(sum(1 for i in r5_rows if par[d[i]] == a[i] or par[a[i]] == d[i])),
    "r5_undirected_pairs": len(cal_r5),
    "path_violation_shapes": n_shapes,
    "path_shape_slot_pairs": len(shape_pairs),
    "path_shape_pairs_present_in_s1": len(cal_path),
    "calibration_pairs_in_s1": int(len(cal)),
    "cells": cells, "chosen": {"tau": TAU, "nmin": NMIN},
}
C.dump(calib, "calibration.json")

inA = lca >= C.THETA0
inB = C.criterion_b(psym, w_raw, TAU, NMIN)
prim = inA & inB
dmask = C.directed_mask(U, prim)
H_unf = C.entropy_by_dp(d, n, np.ones(len(d), bool))
H_f = C.entropy_by_dp(d, n, dmask)
MS_f = C.max_share_score(d, n, dmask)
is_dp = dtot > 0
size = C.compounds_as_dp()
label_freq = size + amarg                   # a direct parent is never also listed as an alt

# ---- S2 ----------------------------------------------------------------------
s2 = C.load_s2(d, a, n)
s2_meta = {"available": s2 is not None}
if s2 is not None:
    x2, y2, c2, s2src = s2
    with np.errstate(divide="ignore", invalid="ignore"):
        p2 = np.log2(c2 * C.N_COMPOUNDS / (amarg[x2] * amarg[y2]))
    l2 = C.lca_depth(A, x2, y2)
    A2, B2 = l2 >= C.THETA0, C.criterion_b(p2, c2, TAU, NMIN)
    P2 = A2 & B2
    s2_meta.update({"source": s2src, "pairs": int(len(c2)), "cooccurrence_increments": int(c2.sum()),
                    "pairs_with_zero_alt_marginal": int(((amarg[x2] == 0) | (amarg[y2] == 0)).sum()),
                    "A4": int(A2.sum()), "B": int(B2.sum()), "primary": int(P2.sum()),
                    "primary_not_in_s1_primary": int(len(set((x2[P2] * C.K + y2[P2]).tolist())
                                                         - set(keys[prim].tolist()))),
                    "written_rows_lca_ge_2": int((l2 >= 2).sum())})
    if (RES_TSV := C.RES / "alt_cooccurrence_pairs.tsv").exists():
        s2_meta["sha256"] = C.sha256(RES_TSV)

# ---- S3 ----------------------------------------------------------------------
s3 = {}
for i in r5_rows.tolist():
    k = (min(d[i], a[i]), max(d[i], a[i]))
    r, kinds = s3.get(k, (0, set()))
    s3[k] = (r + n[i], kinds | {"r5"})
for kind, pairs in (("r6", C.r6_pairs(tax)), ("path", shape_pairs)):
    for k, r in pairs.items():
        rr, kinds = s3.get(k, (0, set()))
        s3[k] = (rr + r, kinds | {kind})

# ---- write -------------------------------------------------------------------
C.OUT.mkdir(parents=True, exist_ok=True)
f6 = lambda x: "" if not np.isfinite(x) else f"{x:.6g}"
f4 = lambda x: "" if not np.isfinite(x) else f"{x:.4f}"
with (C.OUT / "overlay_edges.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["u", "v", "src", "kind", "w_raw", "w_cond", "w_pmi", "lca_depth", "sibling",
                "n_uv", "n_vu", "A4", "B", "primary"])
    for i in range(U.m):
        w.writerow([U.u[i], U.v[i], 1, "dp_alt", int(w_raw[i]), f6(w_cond[i]), f4(w_pmi[i]),
                    lca[i], int(sib[i]), int(n_uv[i]), int(n_vu[i]), int(inA[i]), int(inB[i]),
                    int(prim[i])])
    if s2 is not None:
        for i in np.flatnonzero(l2 >= 2):
            w.writerow([x2[i], y2[i], 2, "alt_alt", int(c2[i]), "", f4(max(0.0, p2[i])), l2[i],
                        int(par[x2[i]] == par[y2[i]]), "", "", int(A2[i]), int(B2[i]), int(P2[i])])
    for (x, y), (r, kinds) in sorted(s3.items()):
        L = int(C.lca_depth(A, np.array([x]), np.array([y]))[0])
        w.writerow([x, y, 3, "+".join(sorted(kinds)), int(r), "", "", L, int(par[x] == par[y]),
                    "", "", int(L >= C.THETA0), "", ""])

deg = np.bincount(U.u[prim], minlength=C.K) + np.bincount(U.v[prim], minlength=C.K)
stren = (np.bincount(U.u[prim], w_pmi[prim], C.K) + np.bincount(U.v[prim], w_pmi[prim], C.K))
with (C.OUT / "overlay_nodes.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["id", "name", "depth", "kingdom", "superclass", "size", "label_freq", "is_dp",
                "dp_alt_total", "alt_marginal", "H_filtered", "H_unfiltered", "maxshare_filtered",
                "degree", "strength"])
    for i in sorted(x for x in tax.parent if x != C.ROOT_N):
        w.writerow([i, tax.name[i], depth[i], tax.anc_at[i][1], sc[i], int(size[i]),
                    int(label_freq[i]), int(is_dp[i]), int(dtot[i]), int(amarg[i]),
                    f"{H_f[i]:.4f}" if is_dp[i] else "", f"{H_unf[i]:.4f}" if is_dp[i] else "",
                    f"{MS_f[i]:.4f}" if is_dp[i] else "", int(deg[i]), f"{stren[i]:.4f}"])

nodes_used = int(((np.bincount(d, minlength=C.K) > 0) | (np.bincount(a, minlength=C.K) > 0)).sum())
nz = int((deg > 0).sum())
meta = {
    "s1": {"directed_rows": int(len(d)), "labels": int(n.sum()), "undirected_pairs": int(U.m),
           "dp_classes": int(is_dp.sum()), "alt_classes": int((amarg > 0).sum()),
           "union_nodes": nodes_used,
           "density_unfiltered_on_union": U.m / (nodes_used * (nodes_used - 1) / 2),
           "both_directions": int(((n_uv > 0) & (n_vu > 0)).sum())},
    "lca_depth_directed": lca_table,
    "separation": {"theta": C.THETA0, "tau": TAU, "nmin": NMIN, "A": int(inA.sum()),
                   "B": int(inB.sum()), "A_and_B": int(prim.sum()),
                   "A_and_B_mass": int(w_raw[prim].sum()),
                   "jaccard": float(prim.sum() / (inA | inB).sum()),
                   "sibling_only": int(sib.sum())},
    "primary_graph": {"nodes_total": len(tax.parent) - 1, "nodes_nonisolated": nz,
                      "edges": int(prim.sum()),
                      "density_on_nonisolated": float(prim.sum() / (nz * (nz - 1) / 2)) if nz > 1 else 0.0},
    "entropy": {"H_unfiltered_mean": float(H_unf[is_dp].mean()),
                "H_unfiltered_median": float(np.median(H_unf[is_dp])),
                "H_unfiltered_max": float(H_unf[is_dp].max()),
                "H_filtered_mean": float(H_f[is_dp].mean()),
                "H_filtered_zero_share": float((H_f[is_dp] == 0).mean()),
                "r_log10size_H_unfiltered": float(np.corrcoef(np.log10(size[is_dp] + 1), H_unf[is_dp])[0, 1]),
                "r_log10size_H_filtered": float(np.corrcoef(np.log10(size[is_dp] + 1), H_f[is_dp])[0, 1])},
    "s2": s2_meta,
    "s3": {"pairs": len(s3), "r5": sum("r5" in k for _, k in s3.values()),
           "r6": sum("r6" in k for _, k in s3.values()),
           "path": sum("path" in k for _, k in s3.values())},
    "calibration_chosen": {"tau": TAU, "nmin": NMIN},
}
C.dump(meta, "overlay_meta.json")
print(json.dumps({k: meta[k] for k in ("s1", "separation", "primary_graph", "entropy", "s2", "s3")}, indent=1))
print("calibration:", json.dumps(cells, indent=0))
