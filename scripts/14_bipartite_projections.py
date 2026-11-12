#!/usr/bin/env python3
"""
Extension A+B — the release as a two-mode network, and what projecting it does.

Newman (Networks: An Introduction, ch. 6, "Bipartite networks") warns that the one-mode
projection of an affiliation network discards information and manufactures structure:
every group becomes a clique, so clustering and density are put in by construction.
ClassyFire is an affiliation network in exactly that sense (compounds x categories), and
the claim is tested here quantitatively rather than asserted.

    A   both degree sequences of B (compound x label set M) and S (compound x substituent),
        their sum identities, and a CSN power-law fit of the class degrees.
    B   the compound-mode projection. Per hierarchy level it is a disjoint union of
        cliques, so every statistic is closed-form in the class sizes (path_counts.tsv).
        For the multi-label M it is not, and density / transitivity are estimated on the
        200k sample by random pairs and triples, before and after removing hub labels.
        Compounds with identical M are structurally equivalent: indistinguishable.

Inputs are the scan's sufficient statistics plus the 200k sample; no 2 GB pass.
"""

import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bipartite_common as bc  # noqa: E402

PAIRS = 1_000_000
TRIPLES = 1_000_000
FOCAL = 1000
HUBS = [0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500]


def two_mode(deg, scan, tax):
    out, modes = {}, []
    rows = scan["rows"]
    for tag, dh, kk, nm_c, nm_k in (("M", deg["degM"], deg["kM"], "compounds (labels M)", "classes"),
                                    ("S", deg["degS"], deg["kS"], "compounds (substituents)",
                                     "substituents")):
        L = np.arange(len(dh))
        E = int(kk.sum())
        assert int((L * dh).sum()) == E, f"sum identity fails for {tag}"
        assert int(dh.sum()) == rows, f"compound histogram does not cover every row ({tag})"
        k = kk[kk > 0]
        out[tag] = {"incidences": E, "compounds": rows, "compounds_degree0": int(dh[0]),
                    "compound_mean_degree": E / rows,
                    "compound_degree_max": int(np.nonzero(dh)[0].max()),
                    "compound_degree_median": int(np.searchsorted(np.cumsum(dh), rows / 2)),
                    "other_nodes": int(k.size), "other_mean_degree": E / k.size,
                    "other_degree_median": float(np.median(k)), "other_degree_max": int(k.max()),
                    "sum_identity": "sum_c k_c == sum_L L*deg[L] == %d (asserted)" % E,
                    "top": [{"id": int(i), "name": tax.nm(int(i)), "k": int(kk[i]),
                             "share_of_compounds": kk[i] / rows}
                            for i in np.argsort(-kk)[:10]]}
        modes += [{"name": nm_c, "n": rows, "mean_degree": E / rows,
                   "degree_hist": [[int(l), int(c)] for l, c in zip(L, dh) if c]},
                  {"name": nm_k, "n": int(k.size), "mean_degree": E / k.size,
                   "degree_hist": bc.log_binned(k)}]

    import powerlaw
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fits = {}
        for tag, kk in (("classes", deg["kM"]), ("substituents", deg["kS"])):
            x = kk[kk > 0]
            f = powerlaw.Fit(x, discrete=True, verbose=False)
            cmp = {}
            for alt in ("lognormal", "exponential", "truncated_power_law"):
                R, p = f.distribution_compare("power_law", alt, normalized_ratio=alt != "truncated_power_law")
                cmp[alt] = {"R": R, "p": p, "favours": ("power_law" if R > 0 else alt)
                            if p < 0.1 else "inconclusive"}
            fits[tag] = {"n": int(x.size), "xmin": float(f.xmin), "alpha": float(f.power_law.alpha),
                         "sigma": float(f.power_law.sigma), "n_tail": int((x >= f.xmin).sum()),
                         "tail_share": float((x >= f.xmin).mean()), "KS_D": float(f.power_law.D),
                         "lognormal_mu": float(f.lognormal.mu),
                         "lognormal_sigma": float(f.lognormal.sigma), "compare": cmp,
                         "note": "R > 0 favours the power law; p is Vuong's significance. "
                                 "truncated_power_law is nested, so its p is the nested LR test."}
    out["powerlaw"] = fits
    return out, modes


def level_projections(tax):
    slots, cnt = bc.path_counts()
    out, chart, projs = {}, [], []
    for L, lev in enumerate(bc.LEVELS):
        v = slots[:, L]
        null = int(cnt[v < 0].sum())
        ids, inv = np.unique(v[v >= 0], return_inverse=True)
        sizes = np.bincount(inv, weights=cnt[v >= 0]).astype(np.int64)
        st = bc.clique_stats(sizes, null)
        padded = int(cnt[(v >= 0) & (v == slots[:, L - 1])].sum()) if L in (3,) else 0
        st.update({"labels": int(ids.size), "null_rows_isolated": null,
                   "padded_rows": padded,
                   "largest": {"id": int(ids[sizes.argmax()]), "name": tax.nm(int(ids[sizes.argmax()])),
                               "size": int(sizes.max())},
                   "degree_hist_logbinned": bc.log_binned([d for d, _ in st["degree_hist"]],
                                                          [c for _, c in st["degree_hist"]])})
        del st["degree_hist"]
        out[lev] = st
        chart.append({"level": lev, "edges": st["edges"], "density": st["density"]})
        projs.append({"name": f"compound x compound, same {lev}", "mode": "compound",
                      "nodes": st["nodes"], "edges": st["edges"], "density": st["density"],
                      "clustering": st["avg_clustering"], "components": st["components"],
                      "note": f"disjoint union of {ids.size} cliques; {null:,} null-slot rows "
                              "isolated; transitivity 1, assortativity 1 (within-clique)"})
    out["assortativity_note"] = ("Every edge joins two members of one clique, so both ends "
                                 "have degree n_c - 1: Pearson degree assortativity is exactly 1 "
                                 "whenever at least two clique sizes occur.")
    out["padding_note"] = ("At the subclass level, rows whose direct parent is at depth 3 repeat "
                           "the class id in the subclass slot (repeat-padding, how_this_was_built "
                           "B.2); they are kept as their own clique keyed by that id "
                           "(padded_rows). Null slots (119,795 rows at subclass) are isolated.")
    return out, chart, projs


def pair_rank(Bs, I, J, rank1):
    """g = 1 + max hub-rank of labels shared by rows I, J (0 if none). An edge survives
    removal of the top-k hubs iff g > k."""
    X = Bs[I].multiply(Bs[J]).tocsr().astype(np.int32)
    X.data = rank1[X.indices]
    return np.asarray(X.max(axis=1).todense()).ravel()


def sample_projection(deg, tax, rng):
    t = time.time()
    Bs, Ss, dps, cid = bc.load_sample()
    n = Bs.shape[0]
    out = {"sample_rows": n, "load_s": round(time.time() - t, 1)}

    # --- representativeness: random sample or the head of the file?
    slots, cnt = bc.path_counts()
    full = np.bincount(slots[:, 4][slots[:, 4] >= 0], weights=cnt[slots[:, 4] >= 0], minlength=bc.N)
    pf = full / full.sum()
    obs = np.bincount(dps[dps >= 0], minlength=bc.N)
    tvd = 0.5 * np.abs(obs / obs.sum() - pf).sum()
    sims = [0.5 * np.abs(rng.multinomial(obs.sum(), pf) / obs.sum() - pf).sum() for _ in range(200)]
    kM_s = np.bincount(Bs.indices, minlength=bc.N)
    c = cid[cid > 0]
    out["representativeness"] = {
        "tvd_direct_parent_vs_full": tvd, "tvd_null_mean": float(np.mean(sims)),
        "tvd_null_sd": float(np.std(sims)), "tvd_null_p": float(np.mean(np.array(sims) >= tvd)),
        "mean_M_sample": float(Bs.nnz / n), "mean_M_full": float(deg["kM"].sum() / bc.ROWS),
        "label_share_corr": float(np.corrcoef(kM_s / n, deg["kM"] / bc.ROWS)[0, 1]),
        "cid_min": int(c.min()), "cid_median": int(np.median(c)), "cid_max": int(c.max()),
        "cid_sorted": bool(np.all(np.diff(c) >= 0)),
        "verdict": "random" if np.mean(np.array(sims) >= tvd) > 0.01 else "not a uniform random sample"}

    # --- simple-graph density and transitivity of the multi-label projection
    order = np.argsort(-deg["kM"], kind="stable")
    rank1 = np.empty(bc.N, np.int32); rank1[order] = np.arange(1, bc.N + 1)
    I = rng.integers(0, n, PAIRS); J = rng.integers(0, n, PAIRS)
    keep = I != J; I, J = I[keep], J[keep]
    g = pair_rank(Bs, I, J, rank1)
    T = rng.integers(0, n, (TRIPLES, 3))
    T = T[(T[:, 0] != T[:, 1]) & (T[:, 0] != T[:, 2]) & (T[:, 1] != T[:, 2])]
    g01, g02, g12 = (pair_rank(Bs, T[:, a], T[:, b], rank1) for a, b in ((0, 1), (0, 2), (1, 2)))
    hub = []
    for k in HUBS:
        e = int((g > k).sum())
        open_ = (g01 > k) & (g02 > k)
        no = int(open_.sum())
        tr = float((open_ & (g12 > k)).sum() / no) if no >= 30 else None
        hub.append({"removed": k, "density": e / len(g), "density_ci99": bc.wilson(e, len(g)),
                    "transitivity": tr, "connected_triples_sampled": no,
                    "removed_last": None if k == 0 else tax.nm(int(order[k - 1]))})
    out["pairs"] = int(len(g)); out["triples"] = int(len(T))
    out["hub_removal"] = hub
    out["top_hubs"] = [{"rank": r + 1, "id": int(i), "name": tax.nm(int(i)),
                        "share": deg["kM"][i] / bc.ROWS} for r, i in enumerate(order[:10])]
    # pairs that share only the hubs: what fraction of edges hang on the top-1 / top-3 labels
    out["edges_only_via_top"] = {str(k): float(((g > 0) & (g <= k)).sum() / max((g > 0).sum(), 1))
                                 for k in (1, 2, 3, 5, 10)}

    # --- projected degree, as a share of the sample, for random focal compounds
    F = rng.choice(n, FOCAL, replace=False)
    frac = []
    for s in range(0, FOCAL, 100):
        X = (Bs @ Bs[F[s:s + 100]].T).tocsc()
        frac.extend(np.diff(X.indptr) / n)
    frac = np.array(frac)
    out["projected_degree_share"] = {"mean": frac.mean(), "median": float(np.median(frac)),
                                     "p01": float(np.quantile(frac, 0.01)), "min": frac.min(),
                                     "focal": FOCAL}

    # --- same for the substituent projection of compounds
    rS = np.empty(bc.N, np.int32); rS[np.argsort(-deg["kS"], kind="stable")] = np.arange(1, bc.N + 1)
    gS = pair_rank(Ss, I, J, rS)
    eS = int((gS > 0).sum())
    out["substituent_projection_density"] = {"density": eS / len(gS), "ci99": bc.wilson(eS, len(gS))}
    return out


def signatures(deg, scan):
    sz, ct = deg["sig_sizes"], deg["sig_count"]
    rows = int((sz * ct).sum())
    assert rows == scan["rows"]
    single = int(ct[sz == 1].sum())
    return {"distinct_label_sets": int(ct.sum()), "singleton_sets": single,
            "compounds_sharing_set": rows - single, "share_indistinguishable": (rows - single) / rows,
            "largest_set": int(sz.max()), "sets_of_size_ge_1000": int(ct[sz >= 1000].sum()),
            "compounds_in_sets_ge_1000": int((sz * ct)[sz >= 1000].sum()),
            "compression": rows / int(ct.sum()),
            "note": "label sets compared by a 64-bit hash of the sorted set (09_scan_bipartite); "
                    "collision probability over 2.4e7 sets ~ 1.6e-5."}


def substituent_projection(deg):
    i, j, v = bc.upper_pairs("substituent_projection")
    nodes = int((deg["kS"] > 0).sum())
    avg, tr = bc.clustering(bc.adjacency(i, j))
    return {"nodes": nodes, "edges": int(len(i)), "density": len(i) / (nodes * (nodes - 1) / 2),
            "avg_clustering": avg, "transitivity": tr}


def main():
    t0 = time.time()
    rng = np.random.default_rng(42)
    scan, deg, tax = bc.scan_summary(), bc.degrees(), bc.taxonomy()
    A, modes = two_mode(deg, scan, tax)
    levels, chart, projs = level_projections(tax)
    k = deg["kM"].astype(np.float64)
    kS = deg["kS"].astype(np.float64)
    multi = {"multigraph_edges": float((k * (k - 1) / 2).sum()),
             "multigraph_edges_exact": str(int(sum(int(x) * (int(x) - 1) // 2 for x in deg["kM"]))),
             "possible_pairs": bc.ROWS * (bc.ROWS - 1) / 2,
             "substituent_multigraph_edges": float((kS * (kS - 1) / 2).sum())}
    smp = sample_projection(deg, tax, rng)
    sig = signatures(deg, scan)
    subp = substituent_projection(deg)
    pairs = bc.ROWS * (bc.ROWS - 1) / 2
    d0, h = smp["hub_removal"][0], smp["hub_removal"]
    projs += [
        {"name": "compound x compound, share any label in M", "mode": "compound", "nodes": bc.ROWS,
         "edges": d0["density"] * pairs, "density": d0["density"], "clustering": d0["transitivity"],
         "components": None,
         "note": f"multigraph edges = sum_c C(k_c,2) = {multi['multigraph_edges']:.3e} (exact); simple "
                 f"density estimated on {smp['pairs']:,} sample pairs, 99% CI {d0['density_ci99']}; "
                 "clustering = transitivity estimated on sample triples"},
        {"name": "compound x compound, share any substituent", "mode": "compound", "nodes": bc.ROWS,
         "edges": smp["substituent_projection_density"]["density"] * pairs,
         "density": smp["substituent_projection_density"]["density"], "clustering": None,
         "components": None, "note": "density estimated on the 200k sample"},
        {"name": "substituent x substituent (F = S^T S)", "mode": "substituent", "nodes": subp["nodes"],
         "edges": subp["edges"], "density": subp["density"], "clustering": subp["avg_clustering"],
         "components": None, "note": f"transitivity {subp['transitivity']:.3f}"}]
    hk = {x["removed"]: x for x in h}
    headline = [
        {"label": "Compound projection density (shared label)", "value": f"{d0['density']:.4f}",
         "note": f"99% CI {d0['density_ci99'][0]:.4f}-{d0['density_ci99'][1]:.4f}; "
                 f"{smp['pairs']:,} random pairs on the 200k sample"},
        {"label": "... after removing top-3 hub labels", "value": f"{hk[3]['density']:.3f}",
         "note": ", ".join(x["name"] for x in smp["top_hubs"][:3])},
        {"label": "Distinct label sets", "value": f"{sig['distinct_label_sets']:,}",
         "note": f"{sig['share_indistinguishable']:.1%} of compounds share their whole label set "
                 "with another compound"}]
    ui = {"headline": headline, "modes": modes, "projections": projs,
          "charts": {"density_by_level": chart,
                     "hub_removal": [{"removed": x["removed"], "density": x["density"],
                                      "transitivity": x["transitivity"]} for x in h],
                     "signature_sizes": [[int(a), int(b)] for a, b in zip(deg["sig_sizes"],
                                                                          deg["sig_count"])]}}
    bc.save_section("A_B_two_mode_and_compound_projection",
                    {"two_mode": A, "compound_projection_by_level": levels,
                     "compound_projection_multilabel": {**multi, **smp},
                     "structural_equivalence": sig, "substituent_projection": subp,
                     "elapsed_s": round(time.time() - t0, 1)}, ui)
    pl = A["powerlaw"]["classes"]
    print(f"A  classes {A['M']['other_nodes']}  <k_c> {A['M']['other_mean_degree']:.0f}  "
          f"<|M|> {A['M']['compound_mean_degree']:.3f}  power law alpha {pl['alpha']:.2f}"
          f"+-{pl['sigma']:.2f} xmin {pl['xmin']:.0f} (tail {pl['n_tail']}) "
          f"vs lognormal R={pl['compare']['lognormal']['R']:.2f} p={pl['compare']['lognormal']['p']:.3f}")
    for lev in bc.LEVELS:
        s = levels[lev]
        print(f"B  {lev:13s} cliques {s['labels']:5d} edges {s['edges']:.3e} density {s['density']:.4f} "
              f"giant {s['giant_share']:.3f} null {s['null_rows_isolated']}")
    print(f"B  multi-label density {d0['density']:.5f} {d0['density_ci99']}  transitivity "
          f"{d0['transitivity']:.4f}  sample {smp['representativeness']['verdict']}")
    for x in h:
        print(f"   hubs removed {x['removed']:4d}  density {x['density']:.4f}  trans {x['transitivity']}"
              f"  (triples {x['connected_triples_sampled']})")
    print(f"B  distinct label sets {sig['distinct_label_sets']:,}  indistinguishable "
          f"{sig['share_indistinguishable']:.3f}   ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
