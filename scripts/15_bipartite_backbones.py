#!/usr/bin/env python3
"""
Extension C — the class-mode projection, and which of its edges are real.

The raw projection C = B^T B links two ChemOnt categories whenever one compound carries
both. At 73.1M compounds that is almost every pair a generic descriptor touches, and the
weights are dominated by ubiquity ('Hydrocarbon derivatives' sits on 99.5 % of rows).
Network measures on it describe the labelling scheme's hubs, not chemistry, so three
statistically validated backbones are extracted and compared:

    SVN        Tumminello et al. (2011) PLoS ONE 6:e17994. Under random co-assignment of two
               labels with degrees k_a, k_b, C_ab ~ Hypergeom(N_rows, k_a, k_b). Bonferroni
               and BH-FDR at alpha = 0.01 over all K(K-1)/2 pairs.
    BiCM       Saracco et al. (2015) Sci Rep 5:10595; projection test Saracco et al. (2017)
               New J Phys 19:053022. Maximum-entropy null fixing BOTH degree sequences, so
               a compound with 40 labels is expected to co-occur with more of them than a
               compound with 5 -- the heterogeneity the hypergeometric null ignores. C_ab is
               Poisson-binomial; Poisson approximation (conservative in the upper tail, see
               test_bipartite), BH-FDR at 0.01.
    BiCM+lift  the PRIMARY backbone: BiCM-validated AND lift = C_ab / E_BiCM[C_ab] >= 2. At
               N = 7.3e7 a 0.4 % excess is ~30 sd, so significance alone keeps trivial edges.
    disparity  Serrano, Boguna, Vespignani (2009) PNAS 106:6483, on Newman's (2001)
               collaboration weights W = B^T diag(1/(|M|-1)) B; local, no global null.

Each backbone gets clustering against a degree-preserving rewiring null, Louvain
communities (seed 42) and their agreement (NMI, ARI) with the taxonomy's own SuperClass
and Class partitions -- the question being whether label co-assignment recovers the tree.
Partitions are also scored by Guimera-Sales-Pardo-Amaral (2007) bipartite modularity M_B on
the full C, whose null is bipartite and so does not reward the cliques every label induces:
the Louvain partitions, the official SuperClass / Class partitions, and shuffled nulls.
"""

import csv
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bipartite_common as bc  # noqa: E402

ALPHA = 0.01
REWIRE_REPS = 10
REWIRE_BUDGET_S = 90           # per backbone; fewer replicates if a replicate is slow
UI_EDGES = 6000
LIFT = 2.0                     # effect-size floor: C_ab / E_BiCM[C_ab] (extension_research 2.6)
PRIMARY = "bicm_lift2"
MB_NULL = 100
HUBS = 20


def raw_projection(i, j, c, k, nodes):
    K = len(nodes)
    s = np.bincount(i, weights=c, minlength=bc.N) + np.bincount(j, weights=c, minlength=bc.N)
    A = bc.adjacency(i, j)
    avg, tr = bc.clustering(A)
    idx = np.full(bc.N, -1); idx[nodes] = np.arange(K)
    Wd = np.zeros((K, K))
    Wd[idx[i], idx[j]] = c; Wd[idx[j], idx[i]] = c
    wc = bc.weighted_clustering_dense(Wd)
    del Wd
    deg = np.asarray(A.sum(1)).ravel()[nodes]
    return {"nodes": K, "edges": int(len(i)), "density": len(i) / (K * (K - 1) / 2),
            "mean_degree": float(deg.mean()), "avg_clustering": avg, "transitivity": tr,
            "weighted_clustering_onnela": wc, "weight_total": float(c.sum()),
            "weight_max": int(c.max()), "weight_median": float(np.median(c)),
            "strength_hist": bc.log_binned(s[nodes].astype(np.int64)),
            "weight_hist": bc.log_binned(c.astype(np.int64)),
            "degree_hist": bc.log_binned(deg.astype(np.int64))}


def svn(i, j, c, k, m):
    t = time.time()
    lp = bc.hypergeom_logsf(c, bc.ROWS, k[i], k[j])
    ratio = c / (k[i].astype(np.float64) * k[j] / bc.ROWS)
    bonf = lp < np.log(ALPHA / m)
    fdr = bc.bh(lp, ALPHA, m)
    return lp, ratio, bonf, fdr, round(time.time() - t, 1)


def bicm(deg, i, j, c, m):
    t = time.time()
    dh = deg["degM"]
    L = np.nonzero(dh)[0]; L = L[L > 0]
    nodes = np.nonzero(deg["kM"])[0]
    fit = bc.bicm_fit(L, dh[L], deg["kM"][nodes])
    assert fit["err_rows"] < 1e-6 and fit["err_cols"] < 1e-6, (fit["err_rows"], fit["err_cols"])
    P = np.zeros((len(L), bc.N)); P[:, nodes] = fit["P"]
    mom = bc.vmotif_moments(P, dh[L].astype(np.float64), i, j)
    lp = bc.poisson_logsf(c, mom[0])
    lp_rna = bc.rna_logsf(c, mom[0], mom[1], mom[2])
    keep = bc.bh(lp, ALPHA, m)
    keep_rna = bc.bh(lp_rna, ALPHA, m)
    info = {"degree_classes": int(len(L)), "max_rel_err_rows": fit["err_rows"],
            "max_rel_err_cols": fit["err_cols"], "converged": fit["converged"],
            "nfev": fit["nfev"], "x_range": [float(np.exp(fit["a"].min())), float(np.exp(fit["a"].max()))],
            "p_max": float(fit["P"].max()),
            "approximation": "Poisson(lambda = sum_r p_ra p_rb) upper tail; sensitivity: "
                             "refined normal (Volkova 1996) with exact mean/var/skew",
            "validated_rna": int(keep_rna.sum()),
            "agreement_poisson_vs_rna_jaccard": float((keep & keep_rna).sum() / max((keep | keep_rna).sum(), 1)),
            "elapsed_s": round(time.time() - t, 1)}
    return lp, mom[0], keep, info


def disparity(deg):
    i, j, w = bc.upper_pairs("class_projection_newman")
    s = np.bincount(i, weights=w, minlength=bc.N) + np.bincount(j, weights=w, minlength=bc.N)
    kk = np.bincount(i, minlength=bc.N) + np.bincount(j, minlength=bc.N)
    a_i = (1 - w / s[i]) ** (kk[i] - 1)
    a_j = (1 - w / s[j]) ** (kk[j] - 1)
    a = np.minimum(a_i, a_j)
    return i, j, w, a, a < ALPHA


def rewire_null(u, v, rng):
    E = len(u)
    vals, reps, t0, acc = [], 0, time.time(), []
    while reps < REWIRE_REPS:
        ts = time.time()
        ru, rv, done = bc.double_edge_swap(u, v, E, rng)
        vals.append(bc.clustering(bc.adjacency(ru, rv)))
        acc.append(done / E)
        reps += 1
        if time.time() - t0 + (time.time() - ts) > REWIRE_BUDGET_S:
            break
    vals = np.array(vals)
    return {"replicates": reps, "swaps_per_edge": float(np.mean(acc)),
            "avg_clustering_mean": float(vals[:, 0].mean()), "avg_clustering_sd": float(vals[:, 0].std()),
            "transitivity_mean": float(vals[:, 1].mean()), "transitivity_sd": float(vals[:, 1].std()),
            "elapsed_s": round(time.time() - t0, 1)}


def louvain(u, v, nodes, sup_all, cls_all, mb, w=None):
    """Louvain (seed 42) on the backbone, unweighted or with edge weights w; its NG
    modularity on the backbone, agreement with the tree, and bipartite modularity M_B."""
    G = nx.Graph(); G.add_nodes_from(nodes.tolist())
    if w is None:
        G.add_edges_from(zip(u.tolist(), v.tolist()))
    else:
        G.add_weighted_edges_from(zip(u.tolist(), v.tolist(), np.asarray(w, float).tolist()))
    wk = None if w is None else "weight"
    comms = sorted(nx.community.louvain_communities(G, weight=wk, seed=42), key=len, reverse=True)
    lab = np.full(bc.N, -1)
    for ci, cm in enumerate(comms):
        lab[list(cm)] = ci
    sup, cls = sup_all[nodes], cls_all[nodes]
    sup_parts = [set(nodes[sup == s].tolist()) for s in np.unique(sup)]
    return {"modularity": nx.community.modularity(G, comms, weight=wk), "communities": len(comms),
            "communities_ge5": int(sum(len(c) >= 5 for c in comms)),
            "modularity_superclass_partition": nx.community.modularity(G, sup_parts, weight=wk),
            "nmi_superclass": bc.nmi(lab[nodes], sup), "ari_superclass": bc.ari(lab[nodes], sup),
            "nmi_class": bc.nmi(lab[nodes], cls), "ari_class": bc.ari(lab[nodes], cls),
            "bimod": mb(lab)}, lab, G


def analyse(name, u, v, rng, sup_all, cls_all, mb):
    t = time.time()
    A = bc.adjacency(u, v)
    avg, tr = bc.clustering(A)
    nodes = np.nonzero(np.asarray(A.sum(1)).ravel())[0]
    K = len(nodes)
    lv, lab, G = louvain(u, v, nodes, sup_all, cls_all, mb)
    comp = list(nx.connected_components(G))
    null = rewire_null(u, v, rng)
    out = {"method": name, "nodes": K, "edges": int(len(u)), "density": len(u) / (K * (K - 1) / 2),
           "mean_degree": 2 * len(u) / K, "avg_clustering": avg, "transitivity": tr,
           "components": len(comp), "giant_share": max(map(len, comp)) / K, **lv,
           "rewiring_null": null,
           "clustering_z": (avg - null["avg_clustering_mean"]) / max(null["avg_clustering_sd"], 1e-12),
           "elapsed_s": round(time.time() - t, 1)}
    return out, lab, nodes


def dp_alt_edges():
    u, v = [], []
    with (bc.ROOT / "results" / "dp_alt_pairs.tsv").open() as fh:
        r = csv.reader(fh, delimiter="\t"); next(r)
        for a, b, _ in r:
            if a and b and a != b:
                u.append(int(a)); v.append(int(b))
    return np.array(u), np.array(v)


def write_edges(name, u, v, w, score):
    o = np.argsort(-score, kind="stable")
    with (bc.OUT / f"backbone_edges_{name}.tsv").open("w") as fh:
        fh.write("u\tv\tweight\tscore\n")
        for x in o:
            fh.write(f"{u[x]}\t{v[x]}\t{w[x]:.6g}\t{score[x]:.6g}\n")


def write_comms(name, lab, nodes, tax, sup_all):
    with (bc.OUT / f"communities_{name}.tsv").open("w") as fh:
        fh.write("id\tname\tcommunity\tsuperclass_id\tsuperclass\n")
        for n in nodes:
            fh.write(f"{n}\t{tax.nm(int(n))}\t{lab[n]}\t{sup_all[n]}\t{tax.nm(int(sup_all[n]))}\n")


def main():
    t0 = time.time()
    rng = np.random.default_rng(42)
    deg, tax = bc.degrees(), bc.taxonomy()
    k = deg["kM"].astype(np.int64)
    full = sparse.load_npz(bc.OUT / "class_projection.npz")
    assert np.array_equal(full.diagonal(), k), "diagonal of C must equal k_c"
    i, j, c = bc.upper_pairs("class_projection")
    c = c.astype(np.int64)
    nodes = np.nonzero(k)[0]
    K = len(nodes); m = K * (K - 1) // 2
    ids = np.arange(bc.N)
    known = np.array([n in tax.anc_at for n in ids])
    sup_all = np.where(known, [tax.anc_at[n].get(2, n) if n in tax.anc_at else n for n in ids], ids)
    cls_all = np.where(known, [tax.anc_at[n].get(3, n) if n in tax.anc_at else n for n in ids], ids)

    raw = raw_projection(i, j, c, k, nodes)
    print(f"raw  nodes {K} edges {raw['edges']:,} density {raw['density']:.4f} "
          f"C {raw['avg_clustering']:.3f} Cw {raw['weighted_clustering_onnela']:.4f}  ({time.time()-t0:.0f}s)")

    lp_h, ratio, bonf, fdr, th = svn(i, j, c, k, m)
    print(f"SVN  bonferroni {bonf.sum():,}  fdr {fdr.sum():,}  ({th}s)")
    lp_b, lam, keep_b, binfo = bicm(deg, i, j, c, m)
    print(f"BiCM validated {keep_b.sum():,} (rna {binfo['validated_rna']:,}) err "
          f"{binfo['max_rel_err_rows']:.1e}/{binfo['max_rel_err_cols']:.1e}  ({binfo['elapsed_s']}s)")
    di, dj, dw, da, keep_d = disparity(deg)
    print(f"disp kept {keep_d.sum():,}")

    sc_h, sc_b = -lp_h / np.log(10), -lp_b / np.log(10)
    lift = c / lam
    keep_l = keep_b & (lift >= LIFT)
    sets = {"svn_bonferroni": (i[bonf], j[bonf], c[bonf], sc_h[bonf]),
            "svn_fdr": (i[fdr], j[fdr], c[fdr], sc_h[fdr]),
            "bicm": (i[keep_b], j[keep_b], c[keep_b], sc_b[keep_b]),
            PRIMARY: (i[keep_l], j[keep_l], c[keep_l], sc_b[keep_l]),
            "disparity": (di[keep_d], dj[keep_d], dw[keep_d], -np.log10(np.maximum(da[keep_d], 1e-300)))}
    mb = lambda lab: bc.bimod(i, j, c, k, lab, deg["degM"])
    res, labs, keys, bnodes = {}, {}, {}, {}
    for name, (u, v, w, s) in sets.items():
        write_edges(name, u, v, w, s)
        out, lab, bn = analyse(name, u, v, rng, sup_all, cls_all, mb)
        res[name], labs[name], bnodes[name] = out, lab, bn
        keys[name] = bc.edge_keys(u, v)
        write_comms(name, lab, bn, tax, sup_all)
        print(f"{name:15s} E {out['edges']:8,} dens {out['density']:.4f} C {out['avg_clustering']:.3f} "
              f"(null {out['rewiring_null']['avg_clustering_mean']:.3f}, {out['rewiring_null']['replicates']} reps) "
              f"Q {out['modularity']:.3f} comms {out['communities']} NMI_sc {out['nmi_superclass']:.3f} "
              f"ARI_sc {out['ari_superclass']:.3f} NMI_c {out['nmi_class']:.3f} M_B {out['bimod']:.4f}"
              f"  ({out['elapsed_s']}s)")

    # weighted Louvain on the primary backbone, w = ln(lift) > ln 2 > 0
    pu, pv = sets[PRIMARY][:2]
    wres, wlab, _ = louvain(pu, pv, bnodes[PRIMARY], sup_all, cls_all, mb, w=np.log(lift[keep_l]))
    wres.update({"weight": "ln(C_ab / E_BiCM[C_ab])", "nmi_vs_unweighted":
                 bc.nmi(wlab[bnodes[PRIMARY]], labs[PRIMARY][bnodes[PRIMARY]])})
    write_comms(PRIMARY + "_weighted", wlab, bnodes[PRIMARY], tax, sup_all)
    print(f"weighted Louvain Q {wres['modularity']:.3f} comms {wres['communities']} "
          f"NMI_sc {wres['nmi_superclass']:.3f} M_B {wres['bimod']:.4f}")

    # bipartite modularity of the official tree, and random-relabelling nulls
    live = k > 0
    tb = {lev: mb(np.where(live, part, -1)) for lev, part in
          (("kingdom", np.where(known, [tax.anc_at[n].get(1, n) if n in tax.anc_at else n for n in ids], ids)),
           ("superclass", sup_all), ("class", cls_all))}
    plab = labs[PRIMARY]; pn = bnodes[PRIMARY]
    null_p = []
    for _ in range(MB_NULL):
        sh = plab.copy(); sh[pn] = rng.permutation(plab[pn]); null_p.append(mb(sh))
    null_s = []
    ln_ = np.nonzero(live)[0]
    for _ in range(MB_NULL):
        sh = np.full(bc.N, -1); sh[ln_] = rng.permutation(sup_all[ln_]); null_s.append(mb(sh))
    null_p, null_s = np.array(null_p), np.array(null_s)
    taxonomy_bimod = {
        **tb, "primary_louvain": res[PRIMARY]["bimod"], "primary_louvain_weighted": wres["bimod"],
        "primary_null_mean": float(null_p.mean()), "primary_null_sd": float(null_p.std()),
        "primary_z": float((res[PRIMARY]["bimod"] - null_p.mean()) / max(null_p.std(), 1e-12)),
        "superclass_null_mean": float(null_s.mean()), "superclass_null_sd": float(null_s.std()),
        "superclass_z": float((tb["superclass"] - null_s.mean()) / max(null_s.std(), 1e-12)),
        "null_replicates": MB_NULL,
        "primary_beats_superclass": bool(res[PRIMARY]["bimod"] > tb["superclass"]),
        "formula": "M_B = sum_s [ sum_{a!=b in s} C_ab / sum_d n_d d(d-1) - sum_{a!=b in s} k_a k_b "
                   "/ (sum_d n_d d)^2 ], Guimera et al. 2007 Eq. 5; off-diagonal C; ordered pairs",
        "note": "evaluated on the full raw C, not on the backbone; nodes outside a backbone are "
                "singletons (contribute 0). Null: module labels shuffled among the same nodes "
                "(sizes preserved)."}
    print(f"M_B kingdom {tb['kingdom']:.4f} superclass {tb['superclass']:.4f} class {tb['class']:.4f}  "
          f"primary {res[PRIMARY]['bimod']:.4f} (null {null_p.mean():.4f}+-{null_p.std():.4f})")

    # what the lift floor removes
    hubs = np.argsort(-k)[:HUBS]
    ishub = np.isin(i, hubs) & np.isin(j, hubs)
    heavy = np.argsort(-c)[:15]
    lift_floor = {"bicm_edges": int(keep_b.sum()), "kept": int(keep_l.sum()),
                  "removed": int((keep_b & ~keep_l).sum()),
                  "removed_share": float((keep_b & ~keep_l).sum() / keep_b.sum()),
                  "hub_definition": f"top {HUBS} categories by k_c",
                  "hub_hub_pairs_raw": int(ishub.sum()), "hub_hub_in_bicm": int((ishub & keep_b).sum()),
                  "hub_hub_in_primary": int((ishub & keep_l).sum()),
                  "hub_hub_primary_pairs": [{"u": tax.nm(int(i[x])), "v": tax.nm(int(j[x])),
                                             "lift": float(lift[x])} for x in np.nonzero(ishub & keep_l)[0]],
                  "heaviest15_raw_in_bicm": int(keep_b[heavy].sum()),
                  "heaviest15_raw_in_primary": int(keep_l[heavy].sum()),
                  "removed_weight_share": float(c[keep_b & ~keep_l].sum() / c[keep_b].sum()),
                  "median_k_endpoint_removed": float(np.median(np.r_[k[i[keep_b & ~keep_l]], k[j[keep_b & ~keep_l]]])),
                  "median_k_endpoint_kept": float(np.median(np.r_[k[i[keep_l]], k[j[keep_l]]]))}
    print(f"lift floor removes {lift_floor['removed']:,} of {lift_floor['bicm_edges']:,}; hub-hub "
          f"{lift_floor['hub_hub_in_bicm']} -> {lift_floor['hub_hub_in_primary']}; heaviest-15 "
          f"{lift_floor['heaviest15_raw_in_bicm']} -> {lift_floor['heaviest15_raw_in_primary']}")

    du, dv = dp_alt_edges()
    keys["dp_alt_overlay_undirected"] = sorted(set(bc.edge_keys(du, dv)))
    raw_keys = set(bc.edge_keys(i, j))
    names = list(keys)
    jac = {a: {b: bc.jaccard(keys[a], keys[b]) for b in names} for a in names}
    da_set = set(keys["dp_alt_overlay_undirected"])
    overlap = {n: {"share_of_backbone_in_overlay": len(set(keys[n]) & da_set) / max(len(keys[n]), 1),
                   "share_of_overlay_in_backbone": len(set(keys[n]) & da_set) / len(da_set)}
               for n in sets}

    # the strongest over-expressions the BiCM keeps, and the heaviest raw edges it drops
    o = np.argsort(-sc_b)
    top = [{"u": tax.nm(int(i[x])), "v": tax.nm(int(j[x])), "C": int(c[x]), "E_bicm": float(lam[x]),
            "ratio_bicm": float(c[x] / lam[x]), "ratio_hyper": float(ratio[x])}
           for x in o[:15]]
    dropped = [{"u": tax.nm(int(i[x])), "v": tax.nm(int(j[x])), "C": int(c[x]), "E_bicm": float(lam[x]),
                "bicm_validated": bool(keep_b[x]), "primary": bool(keep_l[x]),
                "svn_bonferroni": bool(bonf[x])} for x in heavy]
    lca_same_sup = {n: float(np.mean(sup_all[sets[n][0]] == sup_all[sets[n][1]])) for n in sets}
    lca_same_sup["raw"] = float(np.mean(sup_all[i] == sup_all[j]))

    payload = {"raw_projection": raw, "tests": m,
               "svn": {"bonferroni_threshold_p": ALPHA / m, "bonferroni": int(bonf.sum()),
                       "fdr": int(fdr.sum()), "elapsed_s": th,
                       "ratio_median_validated": float(np.median(ratio[bonf])) if bonf.any() else None,
                       "ratio_median_all": float(np.median(ratio)),
                       "share_of_raw_edges_with_ratio_gt1": float((ratio > 1).mean()),
                       "tail": "exact log-pmf + summed tail with geometric remainder bound "
                               "(conservative); scipy hypergeom.logsf ~1 ms/pair at N=7.3e7"},
               "bicm": {**binfo, "validated": int(keep_b.sum()),
                        "ratio_median_validated": float(np.median(c[keep_b] / lam[keep_b])),
                        "validated_by_effect_size": {f"obs/exp<{t}": int((c[keep_b] / lam[keep_b] < t).sum())
                                                     for t in (1.1, 1.5, 2, 5)},
                        "effect_size_note": "at N = 7.3e7 a 0.4 % excess on a hub pair is ~30 sd; "
                                            "validated != large. Counts above split the backbone by "
                                            "observed / BiCM-expected co-occurrence."},
               "disparity": {"alpha": ALPHA, "kept": int(keep_d.sum()), "on": "Newman-weighted W"},
               "primary": PRIMARY, "lift_floor": lift_floor,
               "primary_weighted_louvain": wres, "taxonomy_bimod": taxonomy_bimod,
               "backbones": res, "jaccard": jac, "dp_alt_overlap": overlap,
               "same_superclass_share": lca_same_sup,
               "bicm_top_edges": top, "heaviest_raw_edges": dropped,
               "notes": {"louvain": "unweighted backbone, networkx louvain_communities seed=42, resolution 1",
                         "partitions": "classes shallower than the level map to themselves",
                         "rewiring": f"bc.double_edge_swap (nx algorithm, reimplemented), 1 accepted "
                                     f"swap per edge per replicate, up to {REWIRE_REPS} replicates or "
                                     f"{REWIRE_BUDGET_S}s per backbone",
                         "dp_alt": "Phase-2 overlay is directed (direct parent -> alternative parent) "
                                   "and excludes alt-alt co-assignment; compared here undirected, "
                                   "self-pairs dropped. It is a subset of the raw projection: "
                                   f"{len(da_set - raw_keys)} overlay pairs absent from C."},
               "elapsed_s": round(time.time() - t0, 1)}

    b = res[PRIMARY]
    bu, bv, bw, bs = sets[PRIMARY]
    o = np.argsort(-bs, kind="stable")[:UI_EDGES]
    lab = labs[PRIMARY]
    comms_ui = []
    for ci in range(b["communities"]):
        mem = np.nonzero(lab == ci)[0]
        vals, cnt = np.unique(sup_all[mem], return_counts=True)
        comms_ui.append({"id": ci, "size": int(len(mem)), "label": tax.nm(int(vals[cnt.argmax()])),
                         "purity": float(cnt.max() / len(mem))})
    ui = {"projections": [{"name": "class x class, raw C = B^T B", "mode": "class", "nodes": K,
                           "edges": raw["edges"], "density": raw["density"],
                           "clustering": raw["avg_clustering"], "components": None,
                           "note": f"weighted (Onnela) clustering {raw['weighted_clustering_onnela']:.4f}; "
                                   f"max weight {raw['weight_max']:,}"}] +
                         [{"name": f"class backbone: {n}", "mode": "class", "nodes": r["nodes"],
                           "edges": r["edges"], "density": r["density"], "clustering": r["avg_clustering"],
                           "components": r["components"], "note": f"Louvain Q {r['modularity']:.3f}"}
                          for n, r in res.items()],
          "backbones": [{"method": n, "edges": r["edges"], "density": r["density"],
                         "clustering": r["avg_clustering"],
                         "clustering_null": r["rewiring_null"]["avg_clustering_mean"],
                         "modularity": r["modularity"], "communities": r["communities"],
                         "nmi_superclass": r["nmi_superclass"], "ari_superclass": r["ari_superclass"],
                         "nmi_class": r["nmi_class"], "ari_class": r["ari_class"], "nodes": r["nodes"],
                         "bimod": r["bimod"], "primary": n == PRIMARY}
                        for n, r in res.items()] +
                       [{"method": PRIMARY + "_weighted", "edges": b["edges"], "density": b["density"],
                         "clustering": b["avg_clustering"],
                         "clustering_null": b["rewiring_null"]["avg_clustering_mean"],
                         "modularity": wres["modularity"], "communities": wres["communities"],
                         "nmi_superclass": wres["nmi_superclass"], "ari_superclass": wres["ari_superclass"],
                         "nmi_class": wres["nmi_class"], "ari_class": wres["ari_class"], "nodes": b["nodes"],
                         "bimod": wres["bimod"], "primary": False,
                         "note": "same edges as the primary; Louvain weighted by ln(lift)"}],
          "edges": [[int(bu[x]), int(bv[x]), int(bw[x]), round(float(bs[x]), 2)] for x in o],
          "nodes": {str(int(n)): {"k": int(k[n]), "comm": int(lab[n])} for n in nodes},
          "communities": comms_ui,
          "headline": [{"label": "Class-projection density (raw)", "value": f"{raw['density']:.3f}",
                        "note": f"{raw['edges']:,} edges among {K:,} classes"},
                       {"label": "BiCM + lift>=2 class edges (primary)", "value": f"{b['edges']:,}",
                        "note": f"{b['edges']/raw['edges']:.1%} of raw; lift floor removes "
                                f"{lift_floor['removed']:,} BiCM edges; NMI vs SuperClass "
                                f"{b['nmi_superclass']:.2f}"},
                       {"label": "Bipartite modularity: co-assignment vs tree",
                        "value": f"{b['bimod']:.3f} vs {tb['superclass']:.3f}",
                        "note": f"M_B of primary Louvain partition vs SuperClass partition "
                                f"(Class {tb['class']:.3f}; shuffled null {null_p.mean():.3f})"}],
          "extra": {"backbone_jaccard": jac, "taxonomy_bimod": taxonomy_bimod,
                    "primary_backbone": PRIMARY, "lift_floor": {x: y for x, y in lift_floor.items()
                                                                if x != "hub_hub_primary_pairs"}}}
    bc.save_section("C_class_projection_backbones", payload, ui)
    print(f"done ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
