#!/usr/bin/env python3
"""
Toy tests for the bipartite extension (scripts/bipartite_common.py).

Each numeric shortcut the analysis takes at 73M-row scale is checked against the slow,
obvious computation on a matrix small enough to brute-force:

    BiCM        grouped-by-degree fit == per-row fit; constraints hold; expected degrees
                and V-motifs equal exhaustive enumeration of all 2^12 matrices
    tails       our hypergeometric / Poisson log-tails vs scipy; Poisson >= exact
                Poisson-binomial tail (why the Poisson approximation is conservative)
    cliques     closed-form projection statistics == networkx on the explicit graph
    SVN         a planted over-expressed pair is the one Bonferroni flags
    NODF        a perfectly nested matrix scores 100; curveball preserves margins
    partitions  NMI / ARI identities; sparse clustering == networkx
"""

import itertools
import sys
import unittest
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bipartite_common as bc  # noqa: E402


class TestBiCM(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        x, y = rng.lognormal(0, 1, 40), rng.lognormal(-1, 1, 15)
        B = rng.random((40, 15)) < (x[:, None] * y) / (1 + x[:, None] * y)
        B = B[B.sum(1) > 0][:, B.sum(0) > 0]
        self.B = B.astype(int)

    def test_grouped_equals_rowwise_and_constraints(self):
        d, k = self.B.sum(1), self.B.sum(0)
        L, n = np.unique(d, return_counts=True)
        g = bc.bicm_fit(L, n, k)
        self.assertLess(g["err_rows"], 1e-8); self.assertLess(g["err_cols"], 1e-8)
        r = bc.bicm_fit(d, np.ones(len(d)), k)
        self.assertLess(r["err_rows"], 1e-8)
        Pg = g["P"][np.searchsorted(L, d)]
        self.assertTrue(np.allclose(Pg, r["P"], atol=1e-8))

    def test_brute_force_enumeration(self):
        B = np.array([[1, 1, 0], [1, 0, 1], [1, 1, 1], [0, 1, 0]])
        d, k = B.sum(1), B.sum(0)
        f = bc.bicm_fit(d, np.ones(4), k)
        P = f["P"]
        Ed, Ek, EV = np.zeros(4), np.zeros(3), np.zeros((3, 3))
        for bits in itertools.product((0, 1), repeat=12):
            M = np.array(bits).reshape(4, 3)
            w = np.prod(np.where(M == 1, P, 1 - P))
            Ed += w * M.sum(1); Ek += w * M.sum(0); EV += w * (M.T @ M)
        self.assertTrue(np.allclose(Ed, d, atol=1e-7)); self.assertTrue(np.allclose(Ek, k, atol=1e-7))
        i, j = np.array([0, 0, 1]), np.array([1, 2, 2])
        mom = bc.vmotif_moments(P, np.ones(4), i, j)
        self.assertTrue(np.allclose(mom[0], EV[i, j], atol=1e-10))

    def test_poisson_conservative_vs_exact_poisson_binomial(self):
        q = np.random.default_rng(3).uniform(0, 0.6, 50)
        pmf = np.array([1.0])
        for p in q:
            pmf = np.convolve(pmf, [1 - p, p])
        for c in range(int(q.sum()) + 2, 40, 3):       # the upper tail, where tests live
            exact = pmf[c:].sum()
            self.assertGreaterEqual(np.exp(bc.poisson_logsf(np.array([c]), np.array([q.sum()]))[0]),
                                    exact)


class TestTails(unittest.TestCase):
    def test_hypergeom_matches_scipy(self):
        rng = np.random.default_rng(0)
        Nt = 50000
        ka, kb = rng.integers(50, 20000, 300), rng.integers(50, 20000, 300)
        mu = ka * kb / Nt
        c = np.maximum(1, (mu * rng.uniform(1.0, 2.5, 300)).astype(int))
        c = np.minimum(c, np.minimum(ka, kb))
        ours = bc.hypergeom_logsf(c, Nt, ka, kb)
        ref = stats.hypergeom.logsf(c - 1, Nt, ka, kb)
        sig = ref < np.log(1e-3)
        self.assertTrue(np.allclose(ours[sig], ref[sig], rtol=1e-6))
        self.assertTrue(np.all(ours >= ref - 1e-9))          # never anti-conservative

    def test_poisson_matches_scipy_and_extends(self):
        lam = np.array([3.0, 50.0, 1e4, 2.0])
        c = np.array([10, 120, 11000, 4000])
        ours = bc.poisson_logsf(c, lam)
        ref = stats.poisson.logsf(c - 1, lam)
        fin = np.isfinite(ref)
        self.assertTrue(np.allclose(ours[fin], ref[fin], rtol=1e-6))
        self.assertTrue(np.isfinite(ours).all() and ours[3] < -1e4)


class TestCliques(unittest.TestCase):
    def test_closed_form_equals_networkx(self):
        sizes, iso = [5, 3, 3, 2, 1, 7], 4
        G = nx.Graph()
        nid = 0
        for s in sizes:
            G.add_nodes_from(range(nid, nid + s))
            G.add_edges_from(itertools.combinations(range(nid, nid + s), 2)); nid += s
        G.add_nodes_from(range(nid, nid + iso))
        cf = bc.clique_stats(sizes, iso)
        n = G.number_of_nodes()
        self.assertEqual(cf["nodes"], n); self.assertEqual(cf["edges"], G.number_of_edges())
        self.assertAlmostEqual(cf["density"], nx.density(G))
        self.assertAlmostEqual(cf["avg_clustering"], nx.average_clustering(G))
        self.assertAlmostEqual(cf["transitivity"], nx.transitivity(G))
        self.assertEqual(cf["components"], nx.number_connected_components(G))
        self.assertAlmostEqual(cf["giant_share"], max(map(len, nx.connected_components(G))) / n)
        self.assertAlmostEqual(cf["assortativity"], nx.degree_assortativity_coefficient(G))
        hist = dict(cf["degree_hist"])
        self.assertEqual(hist, dict(zip(*np.unique([d for _, d in G.degree()], return_counts=True))))


class TestSVN(unittest.TestCase):
    def test_planted_pair_flagged(self):
        rng = np.random.default_rng(7)
        R, K = 5000, 30
        B = rng.random((R, K)) < rng.uniform(0.02, 0.3, K)
        planted = rng.random(R) < 0.15
        B[planted, 3] = True; B[planted, 11] = True
        C = B.T.astype(np.int64) @ B
        i, j = np.triu_indices(K, 1)
        k = B.sum(0)
        lp = bc.hypergeom_logsf(C[i, j], R, k[i], k[j])
        hits = {(a, b) for a, b in zip(i[lp < np.log(0.01 / len(i))], j[lp < np.log(0.01 / len(i))])}
        self.assertEqual(hits, {(3, 11)})

    def test_bh(self):
        lp = np.log(np.array([1e-6, 0.004, 0.03, 0.5]))
        self.assertEqual(bc.bh(lp, 0.05).tolist(), [True, True, True, False])


class TestNested(unittest.TestCase):
    def test_perfectly_nested_is_100(self):
        A = np.triu(np.ones((8, 8)))[:, ::-1]
        self.assertAlmostEqual(bc.nodf(A)[0], 100.0)
        self.assertLess(bc.nodf(np.eye(6))[0], 1e-9)

    def test_curveball_preserves_margins(self):
        rng = np.random.default_rng(0)
        A = (rng.random((40, 25)) < 0.3).astype(int)
        R = bc.curveball(A, rng)
        self.assertTrue((R.sum(0) == A.sum(0)).all() and (R.sum(1) == A.sum(1)).all())
        self.assertFalse((R == A).all())

    def test_fitness_complexity_ranks_diversified_first(self):
        A = np.triu(np.ones((6, 6)))[:, ::-1]
        F, Q = bc.fitness_complexity(A, it=50)
        self.assertTrue(np.all(np.diff(F) <= 1e-12))


class TestBipartiteModularity(unittest.TestCase):
    def test_block_diagonal_true_partition_beats_random(self):
        rng = np.random.default_rng(4)
        blocks, per, rows = 3, 6, 300
        B = np.zeros((rows, blocks * per), int)
        for r in range(rows):
            b = r % blocks
            B[r, b * per + rng.choice(per, rng.integers(2, per + 1), replace=False)] = 1
        C = B.T @ B
        k = B.sum(0)
        i, j = np.triu_indices(B.shape[1], 1)
        keep = C[i, j] > 0
        i, j, c = i[keep], j[keep], C[i, j][keep]
        dh = np.bincount(B.sum(1))
        true = np.repeat(np.arange(blocks), per)
        mt = bc.bimod(i, j, c, k, true, dh)
        null = [bc.bimod(i, j, c, k, rng.permutation(true), dh) for _ in range(200)]
        self.assertGreater(mt, max(null))
        # one module: first term is exactly 1 (every ordered co-assignment pair is internal)
        one = bc.bimod(i, j, c, k, np.zeros(len(k), int), dh)
        self.assertAlmostEqual(one, 1 - (k.sum() ** 2 - (k ** 2).sum()) / k.sum() ** 2)
        self.assertAlmostEqual(bc.bimod(i, j, c, k, -np.ones(len(k), int), dh), 0.0)
        # perfectly separated blocks: M_B = 1 - sum_s (K_s^2 - sum k^2) / (sum k)^2
        K = np.bincount(true, weights=k); K2 = np.bincount(true, weights=k ** 2)
        self.assertAlmostEqual(mt, 1 - (K ** 2 - K2).sum() / k.sum() ** 2)


class TestPartitions(unittest.TestCase):
    def test_nmi_ari(self):
        a = np.array([0, 0, 1, 1, 2, 2])
        self.assertAlmostEqual(bc.nmi(a, a[::-1] * 7), 1.0)
        self.assertAlmostEqual(bc.ari(a, a + 3), 1.0)
        b = np.array([0, 1, 0, 1, 0, 1])
        self.assertLess(abs(bc.ari(a, b)), 0.5)

    def test_sparse_clustering_equals_networkx(self):
        G = nx.gnm_random_graph(200, 900, seed=2)
        G.remove_nodes_from([n for n, d in dict(G.degree()).items() if d == 0])
        u, v = map(np.array, zip(*G.edges()))
        avg, tr = bc.clustering(bc.adjacency(u, v, 200))
        self.assertAlmostEqual(avg, nx.average_clustering(G))
        self.assertAlmostEqual(tr, nx.transitivity(G))

    def test_double_edge_swap_preserves_degrees(self):
        G = nx.barabasi_albert_graph(300, 4, seed=1)
        u, v = map(np.array, zip(*G.edges()))
        ru, rv, done = bc.double_edge_swap(u, v, len(u), np.random.default_rng(0))
        deg = lambda a, b: np.bincount(np.r_[a, b], minlength=300)
        self.assertTrue((deg(u, v) == deg(ru, rv)).all())
        keys = bc.edge_keys(ru, rv)
        self.assertEqual(len(set(keys)), len(keys)); self.assertTrue((ru != rv).all())
        self.assertEqual(done, len(u))
        self.assertLess(bc.jaccard(bc.edge_keys(u, v), keys), 0.5)

    def test_weighted_clustering_equals_networkx(self):
        G = nx.gnm_random_graph(40, 200, seed=5)
        rng = np.random.default_rng(5)
        for a, b in G.edges():
            G[a][b]["weight"] = float(rng.integers(1, 50))
        W = nx.to_numpy_array(G, weight="weight")
        self.assertAlmostEqual(bc.weighted_clustering_dense(W),
                               np.mean(list(nx.clustering(G, weight="weight").values())))


if __name__ == "__main__":
    unittest.main(verbosity=1)
