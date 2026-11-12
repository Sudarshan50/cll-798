#!/usr/bin/env python3
"""
Phase 2, step 5 -- network measures on the competing overlay, with nulls (plan §5, §6.1).

Only the measures the plan's triage calls informative are computed, each against the null
that can actually falsify it: degree-preserving rewiring (weights permuted over the rewired
edges) is the decisive null for H1b; Erdos-Renyi is shown once to fail; a Chung-Lu null on
label frequency tests "the overlay is just class size".

One structural fact governs the reading of everything here: criterion A(theta) keeps only
pairs whose LCA is at depth >= theta, so for every theta >= 2 no edge crosses a SuperClass
(and at theta = 4 none crosses a SubClass). Communities on such a graph are nested in the
tree by construction. Two extra comparisons are therefore added: a tree-constrained
rewiring null (swaps only within a SubClass block, so the theta-constraint is preserved)
and the tree-free B-only overlay, the only variant on which H1c can fail or pass on merit.
"""
import csv, json, sys, time
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase2_common as C  # noqa: E402

t0 = time.time()
R_DP, R_TREE, R_ER, R_CL, R_SWEEP, R_VAR = 200, 100, 50, 50, 5, 20
tax = C.taxonomy()
A = C.anc_matrix(tax, C.K)
sc = C.superclass_label(tax)
meta = json.loads((C.OUT / "overlay_meta.json").read_text())
TAU, NMIN = meta["separation"]["tau"], meta["separation"]["nmin"]
nodes = {int(r["id"]): r for r in csv.DictReader((C.OUT / "overlay_nodes.tsv").open(), delimiter="\t")}
depth = {i: int(r["depth"]) for i, r in nodes.items()}
size = {i: int(r["size"]) for i, r in nodes.items()}
lfreq = {i: int(r["label_freq"]) for i, r in nodes.items()}
dtot = {i: int(r["dp_alt_total"]) for i, r in nodes.items()}
nm = lambda i: tax.name[i]

# ---- read the overlay ----------------------------------------------------------
s1, s2 = [], []
with (C.OUT / "overlay_edges.tsv").open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["src"] == "1":
            s1.append((int(r["u"]), int(r["v"]), float(r["w_pmi"] or 0), int(r["w_raw"]),
                       int(r["lca_depth"]), int(r["n_uv"]), int(r["n_vu"]), r["primary"] == "1"))
        elif r["src"] == "2" and r["primary"] == "1":
            s2.append((int(r["u"]), int(r["v"]), float(r["w_pmi"]), int(r["w_raw"])))
E = np.array([e[:7] for e in s1], dtype=np.float64)
U_, V_, WP, WR, LCA = E[:, 0].astype(int), E[:, 1].astype(int), E[:, 2], E[:, 3], E[:, 4].astype(int)
PRIM = np.array([e[7] for e in s1])
# thresholds use the exact PMI (the file stores w_pmi to 4 dp, which moves ~3 pairs across tau)
_d, _a, _n = C.load_s1()
_U = C.Undirected(_d, _a)
assert np.array_equal(_U.u, U_) and np.array_equal(_U.v, V_), "overlay row order changed"
PSYM = np.nan_to_num(_U.collapse(_n, C.pmi_directed(_d, _a, _n)[0])[2], nan=-np.inf)
assert np.array_equal((LCA >= C.THETA0) & C.criterion_b(PSYM, WR, TAU, NMIN), PRIM)
print(f"read overlay: {time.time()-t0:.0f}s")


def graph(mask, extra=()):
    G = nx.Graph()
    for i in np.flatnonzero(mask):
        G.add_edge(int(U_[i]), int(V_[i]), weight=float(WP[i]), raw=float(WR[i]))
    for u, v, w, r in extra:
        if G.has_edge(u, v):
            G[u][v]["weight"] = max(G[u][v]["weight"], w)
        else:
            G.add_edge(u, v, weight=w, raw=float(r))
    return G


def louvain(G, seed=0):
    comms = nx.community.louvain_communities(G, weight="weight", resolution=1, seed=seed)
    return comms, nx.community.modularity(G, comms, weight="weight")


def sym_corr(G, val):
    x = np.array([(val[u], val[v]) for u, v in G.edges()], dtype=float)
    if len(x) < 3:
        return None
    a, b = np.concatenate([x[:, 0], x[:, 1]]), np.concatenate([x[:, 1], x[:, 0]])
    return float(np.corrcoef(a, b)[0, 1])


logsz = {i: np.log10(lfreq[i]) if lfreq[i] > 0 else 0.0 for i in nodes}


def measures(G, comm_seed=0, full=True):
    out = {"transitivity": nx.transitivity(G), "avg_clustering": nx.average_clustering(G),
           "weighted_clustering": nx.average_clustering(G, weight="weight")}
    comms, q = louvain(G, comm_seed)
    out["modularity"] = q
    out["communities"] = len(comms)
    if full:
        out["degree_assortativity"] = nx.degree_assortativity_coefficient(G)
        out["max_degree"] = max(d for _, d in G.degree())
        degs = np.array([d for _, d in G.degree()])
        out["degree_cv"] = float(degs.std() / degs.mean())
        out["size_assortativity"] = sym_corr(G, logsz)
        out["within_superclass_share"] = float(np.mean([sc[u] == sc[v] for u, v in G.edges()]))
    return out, comms


def permute_weights(G, rng, src):
    """Weights of the source graph `src`, shuffled over the (rewired) edges of G."""
    w = rng.permutation([d["weight"] for _, _, d in src.edges(data=True)])
    for (u, v), x in zip(G.edges(), w):
        G[u][v]["weight"] = float(x)
    return G


def tree_swap(G, rng, nswap):
    """Double-edge swaps restricted to edge pairs inside one SubClass block, so every
    rewired edge still has LCA depth >= 4."""
    E = [tuple(e) for e in G.edges()]
    S = set(E) | {(v, u) for u, v in E}
    grp = defaultdict(list)
    for i, (u, v) in enumerate(E):
        grp[int(A[u, 4])].append(i)
    blocks = [g for g in grp.values() if len(g) >= 2]
    if not blocks:
        return G.copy()
    p = np.array([len(g) for g in blocks], float); p /= p.sum()
    bi = rng.choice(len(blocks), nswap, p=p)
    r1, r2, flip = rng.random(nswap), rng.random(nswap), rng.random(nswap) < 0.5
    for k in range(nswap):
        g = blocks[bi[k]]
        i, j = g[int(r1[k] * len(g))], g[int(r2[k] * len(g))]
        if i == j:
            continue
        (a, b), (c, d) = E[i], E[j]
        if flip[k]:
            c, d = d, c
        if a == d or c == b or (a, d) in S or (c, b) in S:
            continue
        for x, y in ((a, b), (c, d)):
            S.discard((x, y)); S.discard((y, x))
        for x, y in ((a, d), (c, b)):
            S.add((x, y)); S.add((y, x))
        E[i], E[j] = (a, d), (c, b)
    H = nx.Graph(); H.add_nodes_from(G); H.add_edges_from(E, weight=1.0)
    return H


def null_stats(obs, reps, keys):
    out = {}
    for k in keys:
        x = np.array([r[k] for r in reps if r.get(k) is not None], float)
        if not len(x) or obs.get(k) is None:
            continue
        sd = float(x.std(ddof=1)) if len(x) > 1 else 0.0
        out[k] = {"observed": float(obs[k]), "null_mean": float(x.mean()), "null_sd": sd,
                  "z": float((obs[k] - x.mean()) / sd) if sd > 0 else None,
                  "p_upper": (1 + int((x >= obs[k]).sum())) / (1 + len(x)),
                  "central95": [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))],
                  "R": int(len(x))}
    return out


KEYS = ["transitivity", "avg_clustering", "weighted_clustering", "modularity",
        "degree_assortativity", "max_degree", "degree_cv", "size_assortativity"]

# ---- primary G_c -----------------------------------------------------------------
G = graph(PRIM)
n_, m_ = G.number_of_nodes(), G.number_of_edges()
obs, comms = measures(G)
comp = sorted(nx.connected_components(G), key=len, reverse=True)
LCC = G.subgraph(comp[0]).copy()
basic = {"nodes_nonisolated": n_, "edges": m_, "density": 2 * m_ / (n_ * (n_ - 1)),
         "components": len(comp), "lcc_nodes": len(comp[0]), "lcc_edges": LCC.number_of_edges(),
         "lcc_name_examples": sorted(nm(i) for i in list(comp[0]))[:8]}
print(f"G_c: {n_} nodes {m_} edges; {time.time()-t0:.0f}s")

# centralities
strength = dict(G.degree(weight="weight"))
deg = dict(G.degree())
eig = nx.eigenvector_centrality_numpy(LCC, weight="weight")
DG = nx.DiGraph()
for i in np.flatnonzero(PRIM):
    u, v, nuv, nvu = int(U_[i]), int(V_[i]), E[i, 5], E[i, 6]
    if nuv > 0:
        DG.add_edge(u, v, weight=nuv / dtot[u])
    if nvu > 0:
        DG.add_edge(v, u, weight=nvu / dtot[v])
pr = nx.pagerank(DG, weight="weight")
btw = nx.betweenness_centrality(G) if basic["density"] <= 0.05 else None
from scipy.stats import spearmanr  # noqa: E402


def top(c, k=10):
    return [{"id": i, "name": nm(i), "value": round(float(v), 5), "superclass": nm(sc[i])}
            for i, v in sorted(c.items(), key=lambda x: -x[1])[:k]]


def rho(c):
    ids = list(c)
    r, p = spearmanr([c[i] for i in ids], [logsz[i] for i in ids])
    return {"spearman_vs_log10_label_freq": round(float(r), 4), "p": float(p)}


cent = {"strength": {"top": top(strength), **rho(strength)},
        "degree": {"top": top(deg), **rho(deg)},
        "eigenvector_lcc": {"top": top(eig), **rho(eig),
                            "note": "computed on the largest component; the graph is disconnected"},
        "pagerank_directed_wcond": {"top": top(pr), **rho(pr)},
        "betweenness": ({"top": top(btw), **rho(btw)} if btw else
                        "omitted: density > 0.05 (plan §8 gate)")}
assort = {"degree": obs["degree_assortativity"],
          "depth_numeric": nx.numeric_assortativity_coefficient(
              nx.set_node_attributes(G, depth, "depth") or G, "depth"),
          "log10_label_freq_pearson": obs["size_assortativity"],
          "superclass": nx.attribute_assortativity_coefficient(
              nx.set_node_attributes(G, {i: sc[i] for i in G}, "sc") or G, "sc"),
          "within_superclass_share": obs["within_superclass_share"],
          "within_subclass_share": float(np.mean([A[u, 4] == A[v, 4] and A[u, 4] >= 0 for u, v in G.edges()])),
          "mean_abs_depth_difference": float(np.mean([abs(depth[u] - depth[v]) for u, v in G.edges()])),
          "same_depth_share": float(np.mean([depth[u] == depth[v] for u, v in G.edges()])),
          "parent_child_share": float(np.mean([tax.parent[u] == v or tax.parent[v] == u for u, v in G.edges()])),
          "sibling_share": float(np.mean([tax.parent[u] == tax.parent[v] for u, v in G.edges()]))}

# ---- nulls -------------------------------------------------------------------------
rng = np.random.default_rng(C.SEED)
dp_reps, tree_reps, er_reps, cl_reps = [], [], [], []
for r in range(R_DP):
    H = G.copy()
    nx.double_edge_swap(H, nswap=10 * m_, max_tries=1000 * m_, seed=int(rng.integers(1 << 30)))
    dp_reps.append(measures(permute_weights(H, rng, G), comm_seed=r)[0])
print(f"degree-preserving null x{R_DP}: {time.time()-t0:.0f}s")
for r in range(R_TREE):
    H = tree_swap(G, rng, 10 * m_)
    tree_reps.append(measures(permute_weights(H, rng, G), comm_seed=r)[0])
print(f"tree-constrained null x{R_TREE}: {time.time()-t0:.0f}s")
for r in range(R_ER):
    H = nx.gnm_random_graph(n_, m_, seed=int(rng.integers(1 << 30)))
    for (u, v), w in zip(H.edges(), rng.permutation([d["weight"] for *_, d in G.edges(data=True)])):
        H[u][v]["weight"] = float(w)
    er_reps.append(measures(H, comm_seed=r, full=False)[0]
                   | {"max_degree": max(d for _, d in H.degree()),
                      "degree_cv": float(np.std([d for _, d in H.degree()]) / np.mean([d for _, d in H.degree()]))})
# Chung-Lu on label frequency over every class that carries a label at all
cl_nodes = np.array(sorted(i for i in nodes if lfreq[i] > 0))
f = np.array([lfreq[i] for i in cl_nodes], float)
iu, iv = np.triu_indices(len(cl_nodes), 1)
ff = f[iu] * f[iv]
lo, hi = 0.0, 1.0 / ff.min()
for _ in range(60):
    c = (lo + hi) / 2
    lo, hi = (c, hi) if np.minimum(1, c * ff).sum() < m_ else (lo, c)
pcl = np.minimum(1, c * ff)
wts = np.array([d["weight"] for *_, d in G.edges(data=True)])
for r in range(R_CL):
    k = np.flatnonzero(rng.random(len(pcl)) < pcl)
    H = nx.Graph()
    H.add_weighted_edges_from(zip(cl_nodes[iu[k]].tolist(), cl_nodes[iv[k]].tolist(),
                                  rng.choice(wts, len(k)).tolist()))
    cl_reps.append(measures(H, comm_seed=r)[0])
del iu, iv, ff, pcl
print(f"ER x{R_ER}, Chung-Lu x{R_CL}: {time.time()-t0:.0f}s")
nulls = {"degree_preserving": null_stats(obs, dp_reps, KEYS),
         "tree_constrained": null_stats(obs, tree_reps, KEYS),
         "erdos_renyi": null_stats(obs, er_reps, KEYS),
         "chung_lu_label_freq": null_stats(obs, cl_reps, KEYS)}

# ---- H1b / H1c -----------------------------------------------------------------
dpn = nulls["degree_preserving"]
zc, zq = dpn["weighted_clustering"]["z"], dpn["modularity"]["z"]
inside = lambda k: dpn[k]["central95"][0] <= obs[k] <= dpn[k]["central95"][1]
h1b = ("supported" if (zc or 0) > 3 and (zq or 0) > 3 else
       "not supported" if inside("weighted_clustering") and inside("modularity") else "inconclusive")
cid = {i: k for k, cm in enumerate(sorted(comms, key=len, reverse=True)) for i in cm}
ids = sorted(G)
lab_sc = [sc[i] for i in ids]
lab_sub = [int(A[i, 4]) if A[i, 4] >= 0 else i for i in ids]
part = [cid[i] for i in ids]
NMI, ARI = C.nmi(part, lab_sc), C.ari(part, lab_sc)
cross = sum(1 for cm in comms if len({sc[i] for i in cm}) > 1)
h1c_rule = "supported" if NMI < 0.7 else "not supported" if NMI >= 0.9 else "inconclusive"
h1c = {"nmi_superclass": NMI, "ari_superclass": ARI, "nmi_subclass": C.nmi(part, lab_sub),
       "communities": len(comms), "communities_crossing_superclass": cross,
       "prereg_rule_verdict": h1c_rule,
       "verdict": "inconclusive",
       "reason": "every A(theta>=2) edge joins two classes of one SuperClass (theta=4: one "
                 "SubClass), so the Louvain partition is a refinement of the tree by "
                 "construction; NMI < 0.7 then measures granularity, not re-cutting. The "
                 "test is informative only on the tree-free B-only overlay (see variants)."}

# spectral bisection of the largest component (deterministic second opinion)
fv = nx.fiedler_vector(LCC, weight="weight", normalized=True, seed=0, method="tracemin_lu")
side = {n: int(x >= 0) for n, x in zip(LCC.nodes(), fv)}
spec = {"lcc_nodes": LCC.number_of_nodes(),
        "side_sizes": [sum(1 for s in side.values() if s == k) for k in (0, 1)],
        "side_top_superclass": [Counter(nm(sc[n]) for n, s in side.items() if s == k).most_common(3)
                                for k in (0, 1)],
        "side_top_subclass": [Counter(nm(int(A[n, 4])) if A[n, 4] >= 0 else nm(n)
                                      for n, s in side.items() if s == k).most_common(3) for k in (0, 1)],
        "ari_vs_louvain_on_lcc": C.ari([side[n] for n in LCC], [cid[n] for n in LCC])}

# small-world on the LCC only (the full graph is fragmented)
Lobs, Cobs = nx.average_shortest_path_length(LCC), nx.average_clustering(LCC)
Ls, Cs = [], []
for r in range(10):
    H = LCC.copy()
    nx.double_edge_swap(H, nswap=10 * H.number_of_edges(), max_tries=1000 * H.number_of_edges(),
                        seed=int(rng.integers(1 << 30)))
    if nx.is_connected(H):
        Ls.append(nx.average_shortest_path_length(H)); Cs.append(nx.average_clustering(H))
small = {"lcc_L": Lobs, "lcc_C": Cobs, "rewired_L": float(np.mean(Ls)) if Ls else None,
         "rewired_C": float(np.mean(Cs)) if Cs else None, "connected_rewires": len(Ls),
         "sigma": float((Cobs / np.mean(Cs)) / (Lobs / np.mean(Ls))) if Ls and np.mean(Cs) > 0 else None}

# power law (Clauset-Shalizi-Newman) on strength and degree
import powerlaw  # noqa: E402
import warnings  # noqa: E402
pl = {}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for k, data, disc in (("strength", np.array(list(strength.values())), False),
                          ("degree", np.array(list(deg.values())), True)):
        fit = powerlaw.Fit(data[data > 0], discrete=disc, verbose=False)
        Rl, pl_ = fit.distribution_compare("power_law", "lognormal", normalized_ratio=True)
        Re, pe = fit.distribution_compare("power_law", "exponential", normalized_ratio=True)
        pl[k] = {"alpha": fit.power_law.alpha, "sigma": fit.power_law.sigma, "xmin": float(fit.xmin),
                 "n_tail": int((data >= fit.xmin).sum()), "n": int((data > 0).sum()),
                 "LR_vs_lognormal": float(Rl), "p_vs_lognormal": float(pl_),
                 "LR_vs_exponential": float(Re), "p_vs_exponential": float(pe),
                 "reading": ("power law favoured over lognormal" if Rl > 0 and pl_ < 0.1 else
                             "lognormal favoured over power law" if Rl < 0 and pl_ < 0.1 else
                             "power law and lognormal not distinguishable")}
print(f"spectral/small-world/powerlaw: {time.time()-t0:.0f}s")


# ---- variants: tree-free B-only, and S1 ∪ S2 ------------------------------------
def variant(Gv, R):
    o, cm = measures(Gv)
    reps = []
    for r in range(R):
        H = Gv.copy()
        mm = H.number_of_edges()
        nx.double_edge_swap(H, nswap=10 * mm, max_tries=1000 * mm, seed=int(rng.integers(1 << 30)))
        reps.append(measures(permute_weights(H, rng, Gv), comm_seed=r, full=False)[0])
    ii = sorted(Gv)
    cc = {i: k for k, c_ in enumerate(cm) for i in c_}
    pt = [cc[i] for i in ii]
    return {"nodes": Gv.number_of_nodes(), "edges": Gv.number_of_edges(), "observed": o,
            "null_degree_preserving": null_stats(o, reps, ["weighted_clustering", "modularity", "transitivity"]),
            "nmi_superclass": C.nmi(pt, [sc[i] for i in ii]), "ari_superclass": C.ari(pt, [sc[i] for i in ii]),
            "communities_crossing_superclass": sum(1 for c_ in cm if len({sc[i] for i in c_}) > 1),
            "largest_communities": [{"size": len(c_), "top_superclasses": Counter(nm(sc[i]) for i in c_).most_common(3)}
                                    for c_ in sorted(cm, key=len, reverse=True)[:8]]}, R


Bm = C.criterion_b(PSYM, WR, TAU, NMIN)
vB, _ = variant(graph(Bm), R_VAR // 2)          # 73k edges: each rewire is ~0.7M swaps
vS2, _ = variant(graph(PRIM, s2), R_VAR)
print(f"variants: {time.time()-t0:.0f}s")

# ---- 36-cell surface: H1b/H1c columns (sensitivity) ------------------------------
sweep = []
for th in C.THETAS:
    for t in C.TAUS:
        for mn in C.NMINS:
            Gs = graph((LCA >= th) & C.criterion_b(PSYM, WR, t, mn))
            o, cm = measures(Gs, full=False)
            reps = []
            for r in range(R_SWEEP):
                H = Gs.copy(); mm = H.number_of_edges()
                nx.double_edge_swap(H, nswap=10 * mm, max_tries=1000 * mm, seed=int(rng.integers(1 << 30)))
                reps.append(measures(permute_weights(H, rng, Gs), comm_seed=r, full=False)[0])
            ns = null_stats(o, reps, ["weighted_clustering", "modularity"])
            ii = sorted(Gs); cc = {i: k for k, c_ in enumerate(cm) for i in c_}
            sweep.append({"theta": th, "tau": t, "nmin": mn, "nodes": Gs.number_of_nodes(),
                          "edges": Gs.number_of_edges(), "modularity": round(o["modularity"], 4),
                          "weighted_clustering": round(o["weighted_clustering"], 4),
                          "z_modularity": round(ns["modularity"]["z"], 2) if ns["modularity"]["z"] else None,
                          "z_wclust": round(ns["weighted_clustering"]["z"], 2) if ns["weighted_clustering"]["z"] else None,
                          "nmi_superclass": round(C.nmi([cc[i] for i in ii], [sc[i] for i in ii]), 4)})
print(f"sweep: {time.time()-t0:.0f}s")

# ---- write -----------------------------------------------------------------------
order = sorted(comms, key=len, reverse=True)
with (C.OUT / "communities.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["id", "name", "community", "superclass", "superclass_name", "subclass_name",
                "depth", "degree", "strength", "size", "spectral_side"])
    for k, cm in enumerate(order):
        for i in sorted(cm):
            w.writerow([i, nm(i), k, sc[i], nm(sc[i]), nm(int(A[i, 4])) if A[i, 4] >= 0 else "",
                        depth[i], deg[i], round(strength[i], 4), size[i], side.get(i, "")])
comm_summary = [{"id": k, "size": len(cm), "label": Counter(nm(sc[i]) for i in cm).most_common(1)[0][0],
                 "purity": round(Counter(sc[i] for i in cm).most_common(1)[0][1] / len(cm), 4),
                 "subclass": Counter(nm(int(A[i, 4])) if A[i, 4] >= 0 else nm(i) for i in cm).most_common(1)[0][0],
                 "examples": sorted(nm(i) for i in cm)[:5]}
                for k, cm in enumerate(order)]
out = {"graph": "primary G_c = S1 A(4) ∩ B(%d,%d), weight w_pmi" % (TAU, NMIN),
       "basic": basic, "observed": obs, "centrality": cent, "assortativity": assort,
       "nulls": nulls, "replicates": {"degree_preserving": R_DP, "tree_constrained": R_TREE,
                                      "erdos_renyi": R_ER, "chung_lu": R_CL, "sweep": R_SWEEP,
                                      "variants": R_VAR},
       "H1b": {"z_weighted_clustering": zc, "z_modularity": zq, "verdict": h1b,
               "tree_constrained_z_weighted_clustering": nulls["tree_constrained"]["weighted_clustering"]["z"],
               "tree_constrained_z_modularity": nulls["tree_constrained"]["modularity"]["z"]},
       "H1c": h1c, "spectral_bisection": spec, "small_world_lcc": small, "power_law": pl,
       "variants": {"B_only_tree_free": vB, "S1_union_S2": vS2},
       "sweep_h1bc": sweep, "communities": comm_summary[:60], "n_communities": len(order),
       "runtime_s": round(time.time() - t0, 1)}
C.dump(out, "network_measures.json")
print(json.dumps({"basic": basic, "observed": obs, "H1b": out["H1b"],
                  "H1c": {k: v for k, v in h1c.items() if k != "reason"},
                  "tree_null": {k: nulls["tree_constrained"][k]["z"] for k in nulls["tree_constrained"]},
                  "assort": assort, "power_law": pl, "small": small,
                  "B_only": {k: vB[k] for k in ("nodes", "edges", "nmi_superclass", "ari_superclass",
                                                 "communities_crossing_superclass")},
                  "B_only_z": {k: v["z"] for k, v in vB["null_degree_preserving"].items()},
                  "S1S2": {k: vS2[k] for k in ("nodes", "edges", "nmi_superclass")},
                  "runtime": out["runtime_s"]}, indent=1, default=float))
