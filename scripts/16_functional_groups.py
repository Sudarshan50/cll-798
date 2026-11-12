#!/usr/bin/env python3
"""
Extension D — the class x functional-group bipartite network.

ChemOnt classifies by skeleton (the largest structural feature decides the direct parent),
while reactivity and process-safety behaviour -- flammability, thermal instability,
toxicity -- are carried by functional groups. D[d, s] counts compounds with direct parent d
that carry substituent s (09_scan_bipartite), so it asks directly: does the structure-based
taxonomy sort compounds by their functional-group chemistry?

    binarise    revealed comparative advantage (Balassa 1965; Hidalgo & Hausmann 2009),
                RCA_ds = (D_ds / sum_s D_ds) / (sum_d D_ds / sum D) >= 1, removes the class-size
                effect; raw presence D > 0 is reported alongside.
    NODF        nestedness (Almeida-Neto et al. 2008) vs a fixed-fixed curveball null
                (Strona et al. 2014): are rare groups found only in classes that have common ones?
    F-C         Fitness-Complexity (Tacchella et al. 2012): class 'fitness' = diversity weighted
                by the exclusivity of its groups; group 'complexity' is bounded by its least fit holder.
    Mantel      cosine similarity of substituent profiles vs tree distance between classes,
                permutation test on a random subset.

Caveat built into the data: the release's substituents are 'functional groups present, with
anything already implied by the classification removed', so the profiles are partly
decorrelated from the class by construction, which biases toward a weaker Mantel signal.
"""

import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bipartite_common as bc  # noqa: E402

NULL_REPS = 30
MANTEL_N = 1000
MANTEL_PERMS = 999
MIN_COMPOUNDS = 100
HAZARD = ["Organic nitro compounds", "Nitroaromatic compounds", "C-nitro compounds", "Azo compounds",
          "Azoxy compounds", "Diazo compounds", "Alpha-diazo esters", "Organic diazonium salts",
          "Organic N-nitroso compounds", "Organic perchlorate salts", "Fluorodinitromethyl compounds",
          "Organometallic compounds", "Organolithium compounds", "Organic peroxides", "Peroxides",
          "Oxiranes", "Epoxides", "Acyl halides", "Isocyanates", "Organic azides", "Azides",
          "Hydrazines and derivatives", "Organic chlorite salts"]


def tree_distance(tax, ids):
    """Path length and LCA depth between all pairs of `ids` from ancestor-at-level arrays."""
    H = 1 + max(tax.depth.values())
    A = np.full((len(ids), H), -1)
    for r, n in enumerate(ids):
        for L, a in tax.anc_at[int(n)].items():
            A[r, L] = a
    same = (A[:, None, :] == A[None, :, :]) & (A[:, None, :] >= 0)
    lca = same.sum(2) - 1                                  # root is level 0
    dep = np.array([tax.depth[int(n)] for n in ids])
    return dep[:, None] + dep[None, :] - 2 * lca, lca


def cosine(X):
    X = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-300)
    return X @ X.T


def mantel(S, Dm, rng, perms):
    iu = np.triu_indices(len(S), 1)
    x, y = S[iu], Dm[iu].astype(np.float64)
    r = np.corrcoef(x, y)[0, 1]
    rs = stats.spearmanr(x, y).statistic
    Rk = np.zeros_like(S); Rk[iu] = stats.rankdata(x); Rk = Rk + Rk.T
    yr = stats.rankdata(y)
    nul, nuls = np.empty(perms), np.empty(perms)
    for t in range(perms):
        p = rng.permutation(len(S))
        nul[t] = np.corrcoef(S[np.ix_(p, p)][iu], y)[0, 1]
        nuls[t] = np.corrcoef(Rk[np.ix_(p, p)][iu], yr)[0, 1]
    return {"pearson_r": float(r), "spearman_r": float(rs), "pairs": int(len(x)),
            "p_one_sided_negative": float((1 + (nul <= r).sum()) / (perms + 1)),
            "p_spearman_one_sided": float((1 + (nuls <= rs).sum()) / (perms + 1)),
            "null_mean": float(nul.mean()), "null_sd": float(nul.std()), "permutations": perms}


def nested(M, rng, label):
    t = time.time()
    obs, rows, cols = bc.nodf(M)
    null = np.array([bc.nodf(bc.curveball(M, rng))[0] for _ in range(NULL_REPS)])
    return {"matrix": label, "shape": list(M.shape), "fill": float(M.mean()), "nodf": obs,
            "nodf_rows": rows, "nodf_cols": cols, "null": "curveball (fixed-fixed), 5*rows trades",
            "null_reps": NULL_REPS, "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "z": float((obs - null.mean()) / max(null.std(), 1e-12)), "elapsed_s": round(time.time() - t, 1)}


def main():
    t0 = time.time()
    rng = np.random.default_rng(42)
    tax, deg = bc.taxonomy(), bc.degrees()
    Dsp = sparse.load_npz(bc.OUT / "dp_substituent.npz").tocsr()
    rs, cs = np.asarray(Dsp.sum(1)).ravel(), np.asarray(Dsp.sum(0)).ravel()
    rid, sid = np.nonzero(rs)[0], np.nonzero(cs)[0]
    gap = int(deg["kS"].sum() - Dsp.sum())
    assert 0 <= gap < 1000, "D must hold every substituent incidence of rows with a direct parent"
    D = Dsp[rid][:, sid].toarray().astype(np.float64)
    tot = D.sum()
    R = (D / D.sum(1, keepdims=True)) / (D.sum(0) / tot)
    Mr = (R >= 1).astype(np.int8)
    Mp = (D > 0).astype(np.int8)

    slots, cnt = bc.path_counts()
    ncomp = np.bincount(slots[:, 4][slots[:, 4] >= 0], weights=cnt[slots[:, 4] >= 0], minlength=bc.N)[rid]

    nest = {"rca": nested(Mr, rng, "RCA >= 1"), "presence": nested(Mp, rng, "D > 0")}
    print(f"NODF rca {nest['rca']['nodf']:.2f} (null {nest['rca']['null_mean']:.2f}+-{nest['rca']['null_sd']:.2f})"
          f"  presence {nest['presence']['nodf']:.2f} (null {nest['presence']['null_mean']:.2f}"
          f"+-{nest['presence']['null_sd']:.2f})  ({time.time()-t0:.0f}s)")

    F, Q = bc.fitness_complexity(Mr)
    F2, Q2 = bc.fitness_complexity(Mr, it=400)
    div, ubi = Mr.sum(1), Mr.sum(0)
    fc = {"iterations": 200,
          "rank_stability_200_vs_400": float(stats.spearmanr(F, F2).statistic),
          "spearman_fitness_diversity": float(stats.spearmanr(F, div).statistic),
          "spearman_fitness_compounds": float(stats.spearmanr(F, ncomp).statistic),
          "spearman_complexity_ubiquity": float(stats.spearmanr(Q, ubi).statistic),
          "fitness_share_top1pct": float(np.sort(F)[::-1][:max(1, len(F) // 100)].sum() / F.sum()),
          "top_classes": [{"id": int(rid[x]), "name": tax.nm(int(rid[x])), "fit": float(F[x]),
                           "diversity": int(div[x]), "compounds": int(ncomp[x])}
                          for x in np.argsort(-F)[:20]],
          "top_groups": [{"id": int(sid[x]), "name": tax.nm(int(sid[x])), "complexity": float(Q[x]),
                          "ubiquity": int(ubi[x]), "compounds": int(cs[sid[x]])}
                         for x in np.argsort(-Q)[:20]],
          "bottom_groups": [{"id": int(sid[x]), "name": tax.nm(int(sid[x])), "complexity": float(Q[x]),
                             "ubiquity": int(ubi[x])} for x in np.argsort(Q)[:10]],
          "ubiquity1_groups": int((ubi == 1).sum()),
          "top20_groups_with_ubiquity1": int((ubi[np.argsort(-Q)[:20]] == 1).sum())}
    # sensitivity: F-C is driven by groups held by a single class; some are the class itself
    selfm = rid[:, None] == sid[None, :]
    fc["self_incidences"] = {"compounds": float(D[selfm].sum()), "classes": int((D[selfm] > 0).sum())}
    for tag, drop in (("no_self", selfm), ("no_self_ubiq_ge5", selfm | (ubi < 5)[None, :])):
        Ms = np.where(drop, 0, Mr)
        Ms = Ms[:, Ms.sum(0) > 0]
        rows = Ms.sum(1) > 0
        Fs, _ = bc.fitness_complexity(Ms[rows])
        fc[tag] = {"classes": int(rows.sum()), "groups": int(Ms.shape[1]),
                   "spearman_vs_primary": float(stats.spearmanr(Fs, F[rows]).statistic),
                   "top_classes": [tax.nm(int(rid[rows][x])) for x in np.argsort(-Fs)[:8]]}

    # --- class-class projection by profile similarity, and the taxonomy
    big = np.nonzero(ncomp >= MIN_COMPOUNDS)[0]
    pick = np.sort(rng.choice(big, min(MANTEL_N, len(big)), replace=False))
    ids = rid[pick]
    dist, lca = tree_distance(tax, ids)
    S_rca = cosine(Mr[pick].astype(np.float64))
    S_share = cosine(D[pick] / D[pick].sum(1, keepdims=True))
    sup = bc.level_map(tax, ids, 2)
    iu = np.triu_indices(len(ids), 1)
    same = (sup[:, None] == sup[None, :])[iu]
    man = {"subset": int(len(ids)), "eligible_classes": int(len(big)),
           "min_compounds": MIN_COMPOUNDS,
           "rca_cosine_vs_path_length": mantel(S_rca, dist, rng, MANTEL_PERMS),
           "share_cosine_vs_path_length": mantel(S_share, dist, rng, MANTEL_PERMS),
           "rca_cosine_vs_lca_depth_spearman": float(stats.spearmanr(S_rca[iu], lca[iu]).statistic),
           "mean_cosine_same_superclass": float(S_rca[iu][same].mean()),
           "mean_cosine_diff_superclass": float(S_rca[iu][~same].mean()),
           "mean_cosine_by_lca_depth": {int(L): float(S_rca[iu][lca[iu] == L].mean())
                                        for L in np.unique(lca[iu]) if (lca[iu] == L).sum() >= 20}}
    m = man["rca_cosine_vs_path_length"]
    print(f"Mantel r {m['pearson_r']:.3f} p {m['p_one_sided_negative']:.4f}  spearman {m['spearman_r']:.3f}"
          f"  same-sc cos {man['mean_cosine_same_superclass']:.3f} vs {man['mean_cosine_diff_superclass']:.3f}")

    # full cosine projection over all direct-parent classes (for the projections table)
    Sall = cosine(Mr.astype(np.float64)); np.fill_diagonal(Sall, 0)
    n = len(rid); pairs = n * (n - 1) / 2
    e_pos = int((Sall > 0).sum() // 2); e_half = int((Sall >= 0.5).sum() // 2)
    A = sparse.csr_matrix(np.triu(Sall >= 0.5, 1))
    ii, jj = A.nonzero()
    avg, tr = bc.clustering(bc.adjacency(rid[ii], rid[jj]))
    del Sall

    # --- functional-group concentration across the taxonomy (process-safety reading)
    supall = bc.level_map(tax, rid, 2)
    usup, inv = np.unique(supall, return_inverse=True)
    Dsup = np.zeros((len(usup), len(sid)))
    np.add.at(Dsup, inv, D)
    share = Dsup / Dsup.sum(0)
    Hn = -(share * np.log(np.where(share > 0, share, 1))).sum(0) / np.log(len(usup))
    comp_share = Dsup.sum(1) / tot
    Hbase = -(comp_share * np.log(comp_share)).sum() / np.log(len(usup))
    name2col = {tax.nm(int(s)): x for x, s in enumerate(sid)}
    qrank = stats.rankdata(-Q)
    haz = []
    for nm in HAZARD:
        if nm not in name2col:
            continue
        x = name2col[nm]
        top = np.argsort(-D[:, x])[:3]
        haz.append({"id": int(sid[x]), "name": nm, "compounds": int(cs[sid[x]]),
                    "ubiquity_rca": int(ubi[x]), "classes_present": int(Mp[:, x].sum()),
                    "complexity": float(Q[x]), "complexity_rank": int(qrank[x]),
                    "top_superclass": tax.nm(int(usup[share[:, x].argmax()])),
                    "top_superclass_share": float(share[:, x].max()),
                    "norm_entropy_over_superclasses": float(Hn[x]),
                    "top_direct_parents": [{"name": tax.nm(int(rid[d])), "share": float(D[d, x] / D[:, x].sum())}
                                           for d in top]})
    conc = {"superclasses": int(len(usup)),
            "median_norm_entropy": float(np.median(Hn)),
            "compound_mass_norm_entropy": float(Hbase),
            "share_groups_top_superclass_gt_0.5": float((share.max(0) > 0.5).mean()),
            "note": "normalised Shannon entropy of a group's compounds over the SuperClasses of their "
                    "direct parents (1 = spread like the compound mass itself would be at most; "
                    "compound_mass_norm_entropy is the entropy of all compounds for reference)"}

    payload = {"matrix": {"direct_parents": int(len(rid)), "substituents": int(len(sid)),
                          "incidences": float(tot), "nonzero": int((D > 0).sum()),
                          "rca_ones": int(Mr.sum()), "gap_vs_kS_sum": gap,
                          "gap_note": "substituent incidences on the 7 rows with no direct parent"},
               "nestedness": nest, "fitness_complexity": fc, "mantel": man,
               "cosine_projection": {"nodes": n, "edges_cos_gt0": e_pos, "density_cos_gt0": e_pos / pairs,
                                     "edges_cos_ge_0.5": e_half, "density_cos_ge_0.5": e_half / pairs,
                                     "avg_clustering_cos_ge_0.5": avg, "transitivity_cos_ge_0.5": tr},
               "concentration": conc, "hazard_groups": haz,
               "rca_choice": "RCA >= 1 (Hidalgo & Hausmann 2009 convention); presence D > 0 as sensitivity",
               "elapsed_s": round(time.time() - t0, 1)}
    fg = {"nodf": nest["rca"]["nodf"], "nodf_null_mean": nest["rca"]["null_mean"],
          "nodf_null_sd": nest["rca"]["null_sd"], "mantel_r": m["pearson_r"],
          "mantel_p": m["p_one_sided_negative"],
          "top_classes": [{k: v for k, v in x.items() if k in ("id", "name", "fit")} for x in fc["top_classes"]],
          "top_groups": [{k: v for k, v in x.items() if k in ("id", "name", "complexity", "ubiquity")}
                         for x in fc["top_groups"]],
          "nodf_presence": nest["presence"]["nodf"], "nodf_presence_null_mean": nest["presence"]["null_mean"],
          "hazard_groups": [{k: x[k] for k in ("id", "name", "compounds", "ubiquity_rca", "complexity_rank",
                                                "top_superclass", "top_superclass_share")} for x in haz]}
    ui = {"fg": fg,
          "nodes": {str(int(r)): {"fit": float(F[x]), "k": int(deg["kM"][r])} for x, r in enumerate(rid)},
          "projections": [{"name": "class x class, cosine of RCA substituent profiles (>= 0.5)",
                           "mode": "class", "nodes": n, "edges": e_half, "density": e_half / pairs,
                           "clustering": avg, "components": None,
                           "note": f"{e_pos:,} pairs have cosine > 0 (density {e_pos / pairs:.3f})"}],
          "headline": [{"label": "Functional-group nestedness (NODF, RCA)", "value": f"{nest['rca']['nodf']:.1f}",
                        "note": f"curveball null {nest['rca']['null_mean']:.1f} +- {nest['rca']['null_sd']:.1f}"},
                       {"label": "Profile similarity vs tree distance (Mantel r)",
                        "value": f"{m['pearson_r']:.3f}",
                        "note": f"one-sided p = {m['p_one_sided_negative']:.3g}, {MANTEL_PERMS} permutations, "
                                f"{len(ids)} classes"}]}
    bc.save_section("D_functional_groups", payload, ui)
    print(f"F-C top class {fc['top_classes'][0]['name']}  top group {fc['top_groups'][0]['name']}  "
          f"cos>=0.5 density {e_half / pairs:.4f}  ({time.time() - t0:.0f}s)")
    for h in haz:
        print(f"   {h['name']:32s} n={h['compounds']:9,} ubiq {h['ubiquity_rca']:4d} Qrank {h['complexity_rank']:3d}"
              f"  top SC {h['top_superclass']} {h['top_superclass_share']:.2f}")


if __name__ == "__main__":
    main()
