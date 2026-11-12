#!/usr/bin/env python3
"""
Phase 2, step 4 -- H1a: do the expert-flagged classes score higher on the overlay's
ambiguity statistic than size-matched controls? (plan §6.2, docs/phase2_prereg.md §4)

The primary test is run exactly once. Its result is written to primary_test.json with the
sha256 of the score vector it used; any later run reloads that file instead of
recomputing, and refuses to proceed if the scores have changed since. The prereg's hash and
the freeze timestamp are taken before the test, and copied into validation.json.

Everything after the primary is secondary and labelled so: the uncontrolled p, the
unfiltered 2x2, the exact-join n=20 set, the alternative score, the gate-fallback A-only
overlay, the bootstrap CI, positive controls, the leakage check (the 800 test compounds,
pulled from the release by CID, subtracted from the counts), the ChEBI axis, and the
36-cell sensitivity surface.

Run `python scripts/12_validate_flagged.py --extract` once first if
results/phase2/evaluation_compounds_release.tsv is missing (one 2 GB pass, ~1 min).
"""
import csv, hashlib, json, sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase2_common as C  # noqa: E402

EVAL = C.OUT / "evaluation_compounds_release.tsv"
if "--extract" in sys.argv:
    test, chebi = C.evaluation_cids()
    rows = C.extract_release_rows(test + [r[0] for r in chebi], EVAL)
    sys.exit(print(f"extracted {len(rows)} release rows -> {EVAL}"))

now = lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
init = dict(l.split(" ", 1) for l in (C.OUT / "prereg_hash_initial.txt").read_text().splitlines())
prereg = {"file": "docs/phase2_prereg.md", "sha256": C.sha256(C.PREREG), "frozen_at": now(),
          "sha256_sections_1_8": init["sha256"], "sections_1_8_written_at": init["written_at"],
          "note": "frozen by hash, not by git commit (commits not authorised); §9 appended "
                  "after calibration, before this script first ran. On re-runs sha256/frozen_at "
                  "are those recorded with the primary test, and sha256_now must equal them"}

tax = C.taxonomy()
nodes = {int(r["id"]): r for r in csv.DictReader((C.OUT / "overlay_nodes.tsv").open(), delimiter="\t")}
meta = json.loads((C.OUT / "overlay_meta.json").read_text())
TAU, NMIN = meta["separation"]["tau"], meta["separation"]["nmin"]
col = lambda k: np.array([float(nodes[i][k] or 0) if i in nodes else 0.0 for i in range(C.K)])
H_f, H_u, MS_f = col("H_filtered"), col("H_unfiltered"), col("maxshare_filtered")
size, lab = col("size"), col("dp_alt_total")
pop = np.array(sorted(i for i, r in nodes.items() if r["is_dp"] == "1"))      # 3,604 classes
pos = {c: j for j, c in enumerate(pop)}
strata = C.deciles(size[pop])
strata_lab = C.deciles(lab[pop])
flat = np.zeros(len(pop), int)

audit = C.flagged_classes(tax)
for r in audit:
    r["is_dp"] = r["id"] in pos
fl = [pos[r["id"]] for r in audit if r["is_dp"]]
fl20 = [pos[r["id"]] for r in audit if r["is_dp"] and r["resolved_by"] == "exact"]
chain = {"fp_fn_markings": sum(r["fp"] + r["fn"] for r in audit), "distinct_names": len(audit),
         "resolved_exact": sum(r["resolved_by"] == "exact" for r in audit),
         "resolved_total": sum(r["id"] is not None for r in audit),
         "direct_parent_classes_primary": len(fl), "direct_parent_classes_exact_join": len(fl20)}
assert chain["resolved_total"] == chain["distinct_names"], "unresolved names: fix before the test"


# ---- PRIMARY (once) -----------------------------------------------------------
score = H_f[pop]
digest = hashlib.sha256(np.round(score, 4).tobytes()).hexdigest()
pfile = C.OUT / "primary_test.json"
if pfile.exists():
    primary = json.loads(pfile.read_text())
    assert primary["score_sha256"] == digest, "scores changed since the primary test was run"
    print(f"primary test already run at {primary['run_at']}; reloaded, not recomputed")
    null = np.array(primary["null"])
    prereg.update(sha256_now=prereg["sha256"], sha256=primary["prereg_sha256"],
                  frozen_at=primary.get("prereg_frozen_at", primary["run_at"]),
                  unchanged_since_primary=prereg["sha256"] == primary["prereg_sha256"])
    if "prereg_frozen_at" not in primary:
        prereg["frozen_at_note"] = ("first run computed the hash in-process immediately before "
                                    "the test; its timestamp was not persisted, run_at is the bound")
else:
    T, null, p = C.perm_test(score, fl, strata, R=10000, seed=C.SEED)
    primary = {"run_at": now(), "prereg_sha256": prereg["sha256"], "prereg_frozen_at": prereg["frozen_at"],
               "score_sha256": digest,
               "statistic": "mean H_filtered over flagged direct-parent classes",
               "n": len(fl), "R": 10000, "seed": C.SEED, "T_obs": float(T),
               "null_mean": float(null.mean()), "null_sd": float(null.std()), "p": p,
               "alpha": 0.05, "verdict": "supported" if p < 0.05 else "not supported",
               "null": np.round(null, 5).tolist()}
    pfile.write_text(json.dumps(primary))
    print(f"PRIMARY  T={T:.4f}  null={null.mean():.4f}±{null.std():.4f}  p={p:.4f}")


# ---- secondaries ----------------------------------------------------------------
def test(sc, idx, st, R=10000):
    T, nl, p = C.perm_test(sc, idx, st, R=R, seed=C.SEED)
    return {"T_obs": round(float(T), 4), "null_mean": round(float(nl.mean()), 4),
            "null_sd": round(float(nl.std()), 4), "p": round(p, 5), "n": len(idx)}


def pct(sc, x):
    return round(float((sc <= x).mean() * 100), 1)


sec = {
    "filtered_uncontrolled": test(score, fl, flat),
    "unfiltered_matched": test(H_u[pop], fl, strata),
    "unfiltered_uncontrolled": test(H_u[pop], fl, flat),
    "filtered_matched_labeltotal_deciles": test(score, fl, strata_lab),
    "unfiltered_matched_labeltotal_deciles": test(H_u[pop], fl, strata_lab),
    "exact20_filtered_matched": test(score, fl20, strata),
    "exact20_filtered_uncontrolled": test(score, fl20, flat),
    "exact20_unfiltered_matched_labeltotal_deciles (pilot definition)": test(H_u[pop], fl20, strata_lab),
    "exact20_unfiltered_uncontrolled (pilot definition)": test(H_u[pop], fl20, flat),
    "altscore_1_minus_maxshare_matched": test(MS_f[pop], fl, strata),
    "altscore_1_minus_maxshare_uncontrolled": test(MS_f[pop], fl, flat),
}
eff, ci = C.bootstrap_effect(score, fl, strata)
eff_u, ci_u = C.bootstrap_effect(H_u[pop], fl, strata)
descr = {"flagged_mean_H_filtered": round(float(score[fl].mean()), 4),
         "overall_mean_H_filtered": round(float(score.mean()), 4),
         "flagged_mean_H_filtered_percentile": pct(score, score[fl].mean()),
         "flagged_nonzero_H_filtered": int((score[fl] > 0).sum()),
         "overall_nonzero_share_H_filtered": round(float((score > 0).mean()), 4),
         "flagged_mean_H_unfiltered": round(float(H_u[pop][fl].mean()), 4),
         "overall_mean_H_unfiltered": round(float(H_u[pop].mean()), 4),
         "flagged_mean_H_unfiltered_percentile": pct(H_u[pop], H_u[pop][fl].mean()),
         "exact20_mean_H_unfiltered": round(float(H_u[pop][fl20].mean()), 4),
         "exact20_mean_H_unfiltered_percentile": pct(H_u[pop], H_u[pop][fl20].mean()),
         "effect_bits_filtered": round(eff, 4), "effect_ci95_filtered": [round(x, 4) for x in ci],
         "effect_bits_unfiltered": round(eff_u, 4), "effect_ci95_unfiltered": [round(x, 4) for x in ci_u]}
by = {v: k for k, v in tax.name.items()}
posctl = {nm: {"H_filtered": round(float(H_f[by[nm]]), 4), "pct_filtered": pct(score, H_f[by[nm]]),
               "H_unfiltered": round(float(H_u[by[nm]]), 4), "pct_unfiltered": pct(H_u[pop], H_u[by[nm]])}
          for nm in ("Sesquiterpenoids", "Diterpenoids")}
# residual size imbalance: flagged vs the pooled matched-control distribution
lg = np.log10(size[pop])
ctrl = np.concatenate([lg[(strata == strata[i]) & ~np.isin(np.arange(len(pop)), fl)] for i in fl])
size_bal = {"flagged_log10size_quartiles": np.round(np.percentile(lg[fl], [25, 50, 75]), 3).tolist(),
            "matched_pool_log10size_quartiles": np.round(np.percentile(ctrl, [25, 50, 75]), 3).tolist(),
            "all_log10size_quartiles": np.round(np.percentile(lg, [25, 50, 75]), 3).tolist(),
            "flagged_deciles": sorted(Counter(strata[fl].tolist()).items())}

# ---- recompute scores under other edge sets ------------------------------------
d, a, n = C.load_s1()
A_ = C.anc_matrix(tax, C.K)
U = C.Undirected(d, a)
lca = C.lca_depth(A_, U.u, U.v)


def scores(nn, pair_rule):
    pm, _, _ = C.pmi_directed(d, a, nn)
    nu, nv, ps = U.collapse(nn, pm)
    keep = pair_rule(ps, nu + nv)
    return C.entropy_by_dp(d, nn, C.directed_mask(U, keep)), C.entropy_by_dp(d, nn, nn > 0)


Hp, _ = scores(n, lambda ps, wr: (lca >= C.THETA0) & C.criterion_b(ps, wr, TAU, NMIN))
assert np.allclose(np.round(Hp[pop], 4), np.round(score, 4), atol=1e-4), "rebuild disagrees with 10_"
HA, _ = scores(n, lambda ps, wr: (lca >= C.THETA0) & (wr > 0))
HB, _ = scores(n, lambda ps, wr: C.criterion_b(ps, wr, TAU, NMIN))
sec["gate_fallback_A4_only_matched"] = test(HA[pop], fl, strata)
sec["gate_fallback_A4_only_uncontrolled"] = test(HA[pop], fl, flat)
sec["B_only_matched"] = test(HB[pop], fl, strata)

# ---- leakage (C5): subtract the 800 test compounds' own rows -------------------
test_cids, chebi_rows = C.evaluation_cids()
rel = {}
for p in csv.reader(EVAL.open(), delimiter="\t"):
    if len(p) >= 6 and p[1].isdigit():
        rel[int(p[1])] = (json.loads(p[4]), json.loads(p[5]))
sub = Counter()
found = [c for c in test_cids if c in rel]
for c in found:
    tree, other = rel[c]
    if tree[4] is not None:
        for x in other.get("alternative_parents") or []:
            sub[(tree[4], x)] += 1
key = d * C.K + a
delta = np.zeros(len(n))
kidx = {k: i for i, k in enumerate(key.tolist())}
missing = 0
for (x, y), c in sub.items():
    i = kidx.get(x * C.K + y)
    if i is None:
        missing += 1
    else:
        delta[i] = c
Hl, Hlu = scores(n - delta, lambda ps, wr: (lca >= C.THETA0) & C.criterion_b(ps, wr, TAU, NMIN))
fl_ids = pop[fl]
leak = {"test_compounds": len(test_cids), "found_in_release_by_cid": len(found),
        "label_instances_removed": int(delta.sum()), "pairs_touched": int((delta > 0).sum()),
        "pairs_not_in_counts": missing,
        "flagged_with_test_compound_as_dp": int(sum(1 for c in found if rel[c][0][4] in set(fl_ids.tolist()))),
        "max_abs_delta_H_filtered_flagged": float(np.abs(Hl[fl_ids] - H_f[fl_ids]).max()),
        "max_abs_delta_H_unfiltered_flagged": float(np.abs(Hlu[fl_ids] - H_u[fl_ids]).max()),
        "delta_T_obs": float(Hl[fl_ids].mean() - H_f[fl_ids].mean()),
        "unchanged_to_3dp": bool(np.all(np.round(Hl[fl_ids], 3) == np.round(H_f[fl_ids], 3)))}

# ---- ChEBI axis (Additional file 5; 20 compounds, not the 20/21 classes) -------
ch = []
for r in chebi_rows:
    cid, chebi_n, missed, sugg = r[0], r[4], r[5], r[6]
    t = rel.get(cid)
    dp = t[0][4] if t else None
    ch.append({"cid": cid, "disagreement": (missed + sugg) / (chebi_n + sugg),
               "direct_parent": tax.name.get(dp) if dp is not None else None,
               "H_filtered": float(H_f[dp]) if dp is not None else None,
               "H_unfiltered": float(H_u[dp]) if dp is not None else None})
ok = [x for x in ch if x["H_filtered"] is not None]
chebi = {"n_compounds": len(ch), "found_in_release": len(ok), "rows": ch,
         "definition": "(missed + suggestions) / (ChEBI terms + suggestions); exploratory"}
if len(ok) >= 5:
    chebi["spearman_H_filtered"] = C.spearman([x["disagreement"] for x in ok], [x["H_filtered"] for x in ok])
    chebi["spearman_H_unfiltered"] = C.spearman([x["disagreement"] for x in ok], [x["H_unfiltered"] for x in ok])

# ---- 36-cell sensitivity surface (after the primary) ---------------------------
pm, _, _ = C.pmi_directed(d, a, n)
nu, nv, ps = U.collapse(n, pm)
sweep = []
for th in C.THETAS:
    for t in C.TAUS:
        for m in C.NMINS:
            keep = (lca >= th) & C.criterion_b(ps, nu + nv, t, m)
            Hc = C.entropy_by_dp(d, n, C.directed_mask(U, keep))[pop]
            T, nl, p = C.perm_test(Hc, fl, strata, R=2000, seed=C.SEED)
            sweep.append({"theta": th, "tau": t, "nmin": m, "T_obs": round(float(T), 4),
                          "null_mean": round(float(nl.mean()), 4), "p_h1a": round(p, 4),
                          "flagged_nonzero": int((Hc[fl] > 0).sum())})

# ---- audit table + histogram ----------------------------------------------------
for r in audit:
    i = r["id"]
    r.update({"size": int(size[i]), "decile": int(strata[pos[i]]) if r["is_dp"] else None,
              "H_filtered": round(float(H_f[i]), 4) if r["is_dp"] else None,
              "H_unfiltered": round(float(H_u[i]), 4) if r["is_dp"] else None,
              "maxshare_filtered": round(float(MS_f[i]), 4) if r["is_dp"] else None,
              "pct_H_filtered": pct(score, H_f[i]) if r["is_dp"] else None,
              "depth": int(tax.depth[i])})
cnt, edges = np.histogram(null, bins=40)
out = {"prereg": prereg, "join_chain": chain, "join_audit": audit,
       "primary": {k: v for k, v in primary.items() if k != "null"},
       "secondary": sec, "descriptive": descr, "positive_controls": posctl,
       "size_balance": size_bal, "leakage": leak, "chebi_secondary": chebi,
       "sensitivity_sweep_h1a": sweep,
       "perm_null_hist": {"bins": np.round(edges, 4).tolist(), "counts": cnt.tolist(),
                          "observed": primary["T_obs"]},
       "config": {"theta": C.THETA0, "tau": TAU, "nmin": NMIN}}
C.dump(out, "validation.json")
print(json.dumps({"chain": chain, "primary": out["primary"], "descriptive": descr,
                  "secondary": sec, "positive_controls": posctl, "leakage": leak,
                  "chebi": {k: v for k, v in chebi.items() if k != "rows"}}, indent=1))
print("sweep p_h1a:", [(s["theta"], s["tau"], s["nmin"], s["p_h1a"]) for s in sweep])
