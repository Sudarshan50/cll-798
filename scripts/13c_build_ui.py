#!/usr/bin/env python3
"""
Phase 2, last step -- the compact frontend payload results/phase2/ui.json, plus
results/phase2/sweep_full.tsv (the 36-cell surface with the separation, H1a and H1b/H1c
columns of 11_, 12_ and 13_ merged into one table).

Pure assembly: every number is read from a JSON/TSV written by an earlier Phase-2 script,
none is computed here except formatting, so the payload cannot drift from the ledger.
"""
import csv, json, math, sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase2_common as C  # noqa: E402

J = lambda f: json.loads((C.OUT / f).read_text())
meta, val, net, ep, sep = (J(f) for f in ("overlay_meta.json", "validation.json",
                                          "network_measures.json", "error_propagation.json",
                                          "separation_summary.json"))
tax = C.taxonomy()
TAU, NMIN = meta["separation"]["tau"], meta["separation"]["nmin"]
r3 = lambda x: None if x is None else round(float(x), 3)

# ---- nodes ---------------------------------------------------------------------
comm = {int(r["id"]): int(r["community"])
        for r in csv.DictReader((C.OUT / "communities.tsv").open(), delimiter="\t")}
flag = {r["id"] for r in val["join_audit"] if r["is_dp"]}
nodes, scatter, scatter_f = {}, [], []
for r in csv.DictReader((C.OUT / "overlay_nodes.tsv").open(), delimiter="\t"):
    i, dp = int(r["id"]), r["is_dp"] == "1"
    nodes[str(i)] = {"H": r3(r["H_filtered"]) if dp else 0, "Hu": r3(r["H_unfiltered"]) if dp else 0,
                     "comm": comm.get(i, -1), "str": r3(r["strength"]), "deg": int(r["degree"]),
                     "flag": int(i in flag), "dp": int(dp)}
    if dp:
        ls = round(math.log10(int(r["size"])), 3)
        scatter.append([ls, r3(r["H_unfiltered"]), int(i in flag)])
        scatter_f.append([ls, r3(r["H_filtered"]), int(i in flag)])

# ---- competing edges: S1 and S2 primary, S3 tier --------------------------------
edges = []
with (C.OUT / "overlay_edges.tsv").open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["primary"] == "1" or r["src"] == "3":
            edges.append([int(r["u"]), int(r["v"]), r3(r["w_pmi"] or 0), int(r["w_raw"]), int(r["src"])])
edges.sort(key=lambda e: (-e[2], -e[3]))
edges = edges[:6000]

# ---- communities ----------------------------------------------------------------
members = defaultdict(list)
for r in csv.DictReader((C.OUT / "communities.tsv").open(), delimiter="\t"):
    members[int(r["community"])].append(r)
communities = []
for k in sorted(members):
    sc = Counter(r["superclass_name"] for r in members[k])
    lab, c = sc.most_common(1)[0]
    sub = Counter(r["subclass_name"] or r["name"] for r in members[k]).most_common(1)[0][0]
    communities.append({"id": k, "size": len(members[k]), "label": lab,
                        "purity": round(c / len(members[k]), 3), "subclass": sub})

# ---- sweep ----------------------------------------------------------------------
key = lambda r: (int(r["theta"]), int(r["tau"]), int(r["nmin"]))
s11 = {key(r): r for r in csv.DictReader((C.OUT / "separation_sweep.tsv").open(), delimiter="\t")}
s12 = {key(r): r for r in val["sensitivity_sweep_h1a"]}
s13 = {key(r): r for r in net["sweep_h1bc"]}
sweep, full = [], []
for k in sorted(s11):
    a, b, c = s11[k], s12.get(k, {}), s13.get(k, {})
    sweep.append({"theta": k[0], "tau": k[1], "nmin": k[2], "edges": int(a["edges"]),
                  "jaccard": float(a["jaccard"]), "p_h1a": b.get("p_h1a"),
                  "stable95": float(a["stable95"]) if a["stable95"] else None,
                  "nmi": c.get("nmi_superclass"), "modularity": c.get("modularity"),
                  "z_modularity": c.get("z_modularity"), "z_wclust": c.get("z_wclust"),
                  "face": [int(a["aryl_thioethers"]), int(a["oxosteroids"]), int(a["depsipeptides_axis"])],
                  "primary": int(a["primary"])})
    full.append({**a, **{f"h1a_{x}": b.get(x) for x in ("T_obs", "null_mean", "p_h1a", "flagged_nonzero")},
                 **{f"net_{x}": c.get(x) for x in ("nodes", "modularity", "weighted_clustering",
                                                   "z_modularity", "z_wclust", "nmi_superclass")}})
with (C.OUT / "sweep_full.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(full[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows({k: ("" if v is None else v) for k, v in r.items()} for r in full)

# ---- hypotheses -------------------------------------------------------------------
P, S, D = val["primary"], val["secondary"], val["descriptive"]
pilot = S["exact20_unfiltered_matched_labeltotal_deciles (pilot definition)"]
n_sweep_sig = sum(1 for s in sweep if s["p_h1a"] is not None and s["p_h1a"] < 0.05)
H1b, H1c = net["H1b"], net["H1c"]
dpn, trn = net["nulls"]["degree_preserving"], net["nulls"]["tree_constrained"]
vB = net["variants"]["B_only_tree_free"]
hyp = [
    {"id": "H1a", "title": "Expert-flagged classes are more ambiguous than size-matched controls",
     "statistic": f"mean filtered entropy H_f over the {P['n']} flagged direct-parent classes",
     "observed": round(P["T_obs"], 4),
     "null": f"size-decile stratified permutation, R={P['R']}; null mean {P['null_mean']:.3f} ± {P['null_sd']:.3f}",
     "p": round(P["p"], 4), "verdict": P["verdict"],
     "detail": (f"Pre-registered, run once. Effect {D['effect_bits_filtered']:+.3f} bits, 95% CI "
                f"[{D['effect_ci95_filtered'][0]:.3f}, {D['effect_ci95_filtered'][1]:.3f}]. Only "
                f"{D['flagged_nonzero_H_filtered']}/{P['n']} flagged classes have any retained competitor; the "
                f"positive controls Sesquiterpenoids/Diterpenoids score 0 on G_c. Unfiltered pilot reproduced: "
                f"p={pilot['p']} (n=20, label-total deciles). Significant in {n_sweep_sig}/36 sweep cells. "
                f"Exploratory: tree-free B-only score p={S['B_only_matched']['p']}.")},
    {"id": "H1b", "title": "The overlay has more clustering and modularity than a degree-preserving null",
     "statistic": "Onnela weighted clustering (w_pmi) and Louvain modularity",
     "observed": round(net["observed"]["weighted_clustering"], 4),
     "null": f"degree-preserving rewiring, weights permuted, R={net['replicates']['degree_preserving']}",
     "p": round(dpn["weighted_clustering"]["p_upper"], 4), "verdict": H1b["verdict"],
     "detail": (f"z(weighted clustering)={H1b['z_weighted_clustering']:.1f}, z(modularity)="
                f"{H1b['z_modularity']:.1f} (Q={net['observed']['modularity']:.3f}). Against a tree-constrained "
                f"rewiring that keeps every edge inside its SubClass: z={H1b['tree_constrained_z_weighted_clustering']:.1f} "
                f"and {H1b['tree_constrained_z_modularity']:.1f}. p is the empirical minimum 1/(R+1).")},
    {"id": "H1c", "title": "Overlay communities are not a relabelling of the SuperClass partition",
     "statistic": "NMI(Louvain communities, SuperClass)", "observed": round(H1c["nmi_superclass"], 4),
     "null": "pre-registered threshold rule: supported iff NMI < 0.7 (no null distribution)",
     "p": None, "verdict": H1c["verdict"],
     "detail": (f"Rule returns '{H1c['prereg_rule_verdict']}', but {H1c['communities_crossing_superclass']} of "
                f"{H1c['communities']} communities cross a SuperClass and NMI vs SubClass is "
                f"{H1c['nmi_subclass']:.3f}: with theta>=2 every edge stays inside one SuperClass, so the "
                f"partition refines the tree by construction. On the tree-free B-only overlay NMI="
                f"{vB['nmi_superclass']:.3f} and {vB['communities_crossing_superclass']} of "
                f"{vB['observed']['communities']} communities cross a SuperClass.")},
]

# ---- headline -------------------------------------------------------------------
pb = ep["headline"]
lk = val["leakage"]
headline = [
    {"label": "Competing edges (primary A∩B)", "value": f"{meta['separation']['A_and_B']:,}",
     "note": f"of {meta['s1']['undirected_pairs']:,} dp–alt pairs; θ=4, τ={TAU}, n_min={NMIN} fixed on the calibration set"},
    {"label": "H1a primary p", "value": f"{P['p']:.3f}",
     "note": f"n={P['n']} flagged classes vs size-matched controls; {P['verdict']}"},
    {"label": "Pilot reproduced (unfiltered)", "value": f"p = {pilot['p']:.4f}",
     "note": "exact-join n=20, label-total deciles; the brief reported 0.0101"},
    {"label": "Tree-free score (exploratory)", "value": f"p = {S['B_only_matched']['p']:.4f}",
     "note": "B-only overlay, same matched design; not pre-registered"},
    {"label": "H1b clustering z", "value": f"{H1b['z_weighted_clustering']:.0f}",
     "note": f"vs degree-preserving null; {H1b['tree_constrained_z_weighted_clustering']:.1f} vs tree-constrained"},
    {"label": "NMI vs SuperClass", "value": f"{H1c['nmi_superclass']:.2f}",
     "note": f"nested in the tree by construction; tree-free B-only {vB['nmi_superclass']:.2f}"},
    {"label": "Swaps reaching SuperClass", "value": f"{pb['B_only_superclass_or_above']*100:.1f}%",
     "note": f"competitor swap (B-only) vs {pb['uniform_superclass_or_above']*100:.1f}% uniform, "
             f"{pb['size_superclass_or_above']*100:.1f}% size-proportional"},
    {"label": "Leakage", "value": f"|ΔH| ≤ {lk['max_abs_delta_H_filtered_flagged']:.4f}",
     "note": f"{lk['found_in_release_by_cid']}/800 test compounds removed; unchanged to 3 dp"},
]

# ---- charts ---------------------------------------------------------------------
null_compare = []
for null, label, keys in (("degree_preserving", "degree-preserving",
                           ("weighted_clustering", "modularity", "transitivity", "avg_clustering")),
                          ("tree_constrained", "tree-constrained", ("weighted_clustering", "modularity", "transitivity")),
                          ("erdos_renyi", "Erdős–Rényi", ("max_degree", "degree_cv", "transitivity")),
                          ("chung_lu_label_freq", "Chung–Lu (label freq.)", ("transitivity", "weighted_clustering", "modularity"))):
    for k in keys:
        x = net["nulls"][null].get(k)
        if x:
            null_compare.append({"measure": f"{k.replace('_', ' ')} · {label}", "observed": r3(x["observed"]),
                                 "null_mean": r3(x["null_mean"]), "null_sd": r3(x["null_sd"]),
                                 "z": r3(x["z"]) if x["z"] is not None else None})
prof = ep["profiles"]
ui = {
    "meta": {"generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "primary": f"S1 dp–alt pairs with LCA depth ≥ 4 AND PMI_sym ≥ {TAU} AND rows ≥ {NMIN} "
                        f"(A∩B; τ, n_min set on the R5/PATH calibration set), weight w_pmi",
             "nodes": len(nodes), "edges_all": meta["s1"]["undirected_pairs"],
             "edges_competing": meta["separation"]["A_and_B"],
             "edges_shipped": len(edges), "edges_s2_competing": meta["s2"].get("primary"),
             "prereg_sha256": val["prereg"]["sha256"], "prereg_frozen_at": val["prereg"]["frozen_at"]},
    "headline": headline, "hypotheses": hyp, "edges": edges, "nodes": nodes,
    "communities": communities, "sweep": sweep,
    "errorProp": {"levels": ep["levels"], "overlay": prof["B_only"]["mean"],
                  "uniform": prof["uniform"]["mean"], "size": prof["size"]["mean"],
                  "overlay_ci": [prof["B_only"]["ci95_lo"], prof["B_only"]["ci95_hi"]],
                  "overlay_primary": prof["primary_A_and_B"]["mean"],
                  "overlay_unfiltered": prof["unfiltered_S1"]["mean"],
                  "note": "overlay = tree-free B-only competitors; on the primary A∩B every swap is "
                          "'none' by construction (θ=4)"},
    "charts": {"lca_depth": [{"depth": r["depth"], "pairs": r["pairs"], "mass": r["mass"]}
                             for r in meta["lca_depth_directed"]],
               "entropy_vs_size": scatter, "entropy_filtered_vs_size": scatter_f,
               "perm_null": val["perm_null_hist"], "null_compare": null_compare},
}
out = C.OUT / "ui.json"
out.write_text(json.dumps(ui, separators=(",", ":"), ensure_ascii=False))
print(f"wrote {out} ({out.stat().st_size/1e6:.2f} MB): {len(nodes)} nodes, {len(edges)} edges, "
      f"{len(communities)} communities, {len(sweep)} sweep cells")
