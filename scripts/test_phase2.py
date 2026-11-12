#!/usr/bin/env python3
"""
Tests for the Phase-2 machinery (phase2_common.py), same plain-assert style as
test_checks.py.

Run:  .venv/bin/python scripts/test_phase2.py           (needs numpy; scipy for spearman)

  L1  toy taxonomy: PMI against hand-computed values, LCA depth, padded label slots,
      undirected collapse, filtered entropy.
  L2  statistics: the stratified permutation test recovers a planted effect, does not
      invent one under the null, and its matching really is within-stratum.
  L3  non-modification gate (plan C4/R8): the Phase-1 outputs are byte-identical to the
      hashes snapshotted before Phase 2 began (results/phase2/phase1_hashes_before.txt).
"""
import math, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont            # noqa: E402
import phase2_common as C  # noqa: E402

PASS, FAIL = [], []


def ok(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  -- ' + detail) if detail and not cond else ''}")


def toy():
    """          root(9999999)
                      |
                   0 kingdom
                      |
                   1 superclass
                 /           \\
             2 classA       3 classB
             /     \\
         4 subA    5 subB
           |
         6 d5
    """
    return chemont.Taxonomy({9999999: -1, 0: 9999999, 1: 0, 2: 1, 3: 1, 4: 2, 5: 2, 6: 4},
                            {9999999: "root", 0: "k", 1: "s", 2: "cA", 3: "cB", 4: "sA", 5: "sB", 6: "d5"})


def l1():
    print("\nL1 UNIT -- toy taxonomy")
    t = toy()
    A, dep = C.anc_matrix(t, 7), C.depth_array(t, 7)
    u, v = np.array([4, 4, 6, 2, 6, 0]), np.array([5, 3, 5, 4, 4, 0])
    ok("LCA depth: siblings sA/sB meet at classA (3)", C.lca_depth(A, u[:1], v[:1])[0] == 3)
    ok("LCA depth: sA vs classB meet at superclass (2)", C.lca_depth(A, u[1:2], v[1:2])[0] == 2)
    ok("LCA depth: d5 vs sB meet at classA (3)", C.lca_depth(A, u[2:3], v[2:3])[0] == 3)
    ok("LCA depth: ancestor pair classA/sA has LCA = classA (3)", C.lca_depth(A, u[3:4], v[3:4])[0] == 3)
    ok("LCA depth: parent-child sA/d5 = depth of sA (4)", C.lca_depth(A, u[4:5], v[4:5])[0] == 4)
    ok("LCA depth: a node with itself = its depth", C.lca_depth(A, u[5:], v[5:])[0] == 1)
    ok("padded slots: classB (depth 3) repeats itself in the subclass slot",
       C.label_slots(A, dep, [3]).tolist() == [[0, 1, 3, 3]])
    ok("first changed level: sA vs classB differ first at class",
       C.first_changed_level(A, dep, np.array([4]), np.array([3]))[0] == 2)
    ok("first changed level: d5 vs sA share all four slots -> none",
       C.first_changed_level(A, dep, np.array([6]), np.array([4]))[0] == 4)

    # PMI by hand: rows (d, a, n)
    d = np.array([4, 4, 5, 5, 3]); a = np.array([5, 0, 4, 0, 0]); n = np.array([30., 10., 5., 15., 40.])
    pmi, dtot, amarg = C.pmi_directed(d, a, n, k=7)
    N = 100.0
    exp = math.log2((30 / 40) / (30 / N))          # P(5|4)=30/40, P(5)=30/100
    ok("PMI(sA, sB) = log2[(30/40)/(30/100)]", abs(pmi[0] - exp) < 1e-12, f"{pmi[0]} vs {exp}")
    ok("PMI of a ubiquitous descriptor is low: PMI(sA, kingdom) < 0",
       pmi[1] < 0 and abs(pmi[1] - math.log2((10 / 40) / (65 / N))) < 1e-12)
    U = C.Undirected(d, a, k=7)
    nu, nv, ps = U.collapse(n, pmi)
    i = int(np.flatnonzero((U.u == 4) & (U.v == 5))[0])
    ok("undirected collapse sums both directions", nu[i] + nv[i] == 35)
    ok("PMI_sym = min of the two directions", abs(ps[i] - min(pmi[0], pmi[2])) < 1e-12)
    j = int(np.flatnonzero((U.u == 0) & (U.v == 3))[0])
    ok("one-directional pair keeps its single PMI", abs(ps[j] - pmi[4]) < 1e-12)
    H = C.entropy_by_dp(d, n, np.ones(5, bool), k=7)
    ok("unfiltered entropy of sA = H(3/4, 1/4)",
       abs(H[4] - (-(0.75 * math.log2(0.75) + 0.25 * math.log2(0.25)))) < 1e-12)
    H1 = C.entropy_by_dp(d, n, np.array([True, False, True, True, True]), k=7)
    ok("filtered entropy with one retained alternative is 0", H1[4] == 0 and H1[3] == 0)
    ok("NMI of identical partitions is 1, ARI 1", abs(C.nmi([0, 0, 1, 1], [5, 5, 7, 7]) - 1) < 1e-12
       and abs(C.ari([0, 0, 1, 1], [5, 5, 7, 7]) - 1) < 1e-12)


def l2():
    print("\nL2 STATISTICS -- planted effect")
    rng = np.random.default_rng(0)
    m = 3000
    size = rng.lognormal(5, 2, m)
    strata = C.deciles(size)
    base = 0.4 * np.log(size) + rng.normal(0, 1, m)        # score confounded by size
    flagged = rng.choice(np.flatnonzero(strata >= 7), 20, replace=False)   # size-biased set
    planted = base.copy(); planted[flagged] += 1.5
    _, _, p_eff = C.perm_test(planted, flagged, strata, R=4000, seed=1)
    ok("planted +1.5 effect is recovered (p < 0.01)", p_eff < 0.01, f"p={p_eff}")
    _, _, p_null = C.perm_test(base, flagged, strata, R=4000, seed=1)
    ok("no planted effect -> size-matched p is not small (p > 0.05)", p_null > 0.05, f"p={p_null}")
    _, _, p_unc = C.perm_test(base, flagged, np.zeros(m, int), R=4000, seed=1)
    ok("...while the uncontrolled test is fooled by the size confound (p < 0.01)", p_unc < 0.01, f"p={p_unc}")
    ps = [C.perm_test(rng.normal(0, 1, 500), rng.choice(500, 20, replace=False),
                      C.deciles(rng.random(500)), R=500, seed=s)[2] for s in range(200)]
    ok("p-values are roughly uniform under H0 (5 % level holds within 2x)",
       np.mean(np.array(ps) < 0.05) < 0.10, f"{np.mean(np.array(ps) < 0.05)}")
    lo, hi = C.bootstrap_effect(planted, flagged, strata)[1]
    ok("bootstrap CI of the planted effect covers 1.5", lo < 1.5 < hi, f"[{lo:.2f}, {hi:.2f}]")


def l3():
    print("\nL3 NON-MODIFICATION -- Phase-1 outputs unchanged")
    snap = C.OUT / "phase1_hashes_before.txt"
    if not snap.exists():
        ok("Phase-1 hash snapshot exists", False, str(snap)); return
    for line in snap.read_text().splitlines():
        h, f = line.split()
        ok(f"unchanged: {f}", C.sha256(C.ROOT / f) == h)


if __name__ == "__main__":
    l1(); l2(); l3()
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("failed: " + ", ".join(FAIL))
    sys.exit(1 if FAIL else 0)
