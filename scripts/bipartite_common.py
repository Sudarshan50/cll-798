"""
Shared machinery for the bipartite extension (scripts 14-16).

Everything here is small, pure and unit-tested on toys (scripts/test_bipartite.py), so the
three analysis scripts read as a sequence of statements about the data rather than as
numerics. The scan (09_scan_bipartite.py) reduced 73.1M rows to category-sized matrices;
nothing in this module touches the 2 GB file.

    loading      scan outputs, taxonomy, names, the 200k sample as sparse incidence
    tails        log upper tails of the hypergeometric and Poisson laws, fast at N = 7e7
    BiCM         bipartite configuration model fitted on degree classes (Saracco 2015)
    graphs       clustering / transitivity from sparse adjacency, partition agreement
    nested       NODF, curveball null, Fitness-Complexity
    output       one results JSON shared by the three scripts, assembled into ui.json
"""

import csv
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
from scipy import optimize, sparse, special

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont  # noqa: E402
import paths  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "results" / "bipartite"
RESULTS = OUT / "bipartite_results.json"
N = 4832
ROWS = 73105281
LEVELS = ("kingdom", "superclass", "class", "subclass", "direct_parent")


# ---- loading ------------------------------------------------------------------------
def scan_summary(require_complete=True):
    s = json.loads((OUT / "bipartite_scan.json").read_text())
    if require_complete and not (s.get("complete") and s["rows"] == ROWS):
        raise SystemExit("bipartite_scan.json is not a complete 73,105,281-row scan")
    return s


def degrees():
    d = np.load(OUT / "degrees.npz")
    return {k: d[k] for k in d.files}


def sym(name):
    """Full symmetric matrix from an upper-triangle CSR (diagonal kept once)."""
    m = sparse.load_npz(OUT / f"{name}.npz").tocsr()
    return (m + sparse.triu(m, 1).T).tocsr()


def upper_pairs(name, k=1):
    """(i, j, v) of the strictly-upper triangle, i < j."""
    m = sparse.triu(sparse.load_npz(OUT / f"{name}.npz"), k).tocoo()
    return m.row.astype(np.int64), m.col.astype(np.int64), m.data


def taxonomy():
    return chemont.from_dictionary_tsv(ROOT / "data" / "chemont_dictionary.tsv")


def level_map(tax, ids, level):
    """Ancestor at `level` (2 = SuperClass, 3 = Class); shallower nodes map to themselves."""
    return np.array([tax.anc_at[int(i)].get(level, int(i)) for i in ids])


def path_counts():
    with (ROOT / "results" / "path_counts.tsv").open() as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)
        rows = [x for x in r]
    slots = np.array([[int(v) if v else -1 for v in x[:5]] for x in rows])
    return slots, np.array([int(x[5]) for x in rows], dtype=np.int64)


def load_sample(limit=None):
    """The 200k sample as sparse incidence: B (M = dp U alts), S (substituents), dp, cid."""
    p = paths.require("sample_200k.tsv")
    mi, mp, si, sp, dps, cid = [], [0], [], [0], [], []
    with p.open() as fh:
        next(fh)
        for n, line in enumerate(fh):
            if limit and n >= limit:
                break
            c = line.rstrip("\n").split("\t")
            t, o = json.loads(c[4]), json.loads(c[5])
            m = set(o.get("alternative_parents") or ())
            if t[4] is not None:
                m.add(t[4])
            mi.extend(m); mp.append(len(mi))
            s = set(o.get("substituents") or ())
            si.extend(s); sp.append(len(si))
            dps.append(-1 if t[4] is None else t[4])
            cid.append(int(c[1]) if c[1] else -1)
    rows = len(dps)

    def csr(i, p):
        return sparse.csr_matrix((np.ones(len(i), np.int8), np.asarray(i, np.int32),
                                  np.asarray(p, np.int64)), shape=(rows, N))
    return csr(mi, mp), csr(si, sp), np.asarray(dps), np.asarray(cid)


# ---- statistics ---------------------------------------------------------------------
def log_binned(values, counts=None, per_decade=8):
    """[[k, count]] histogram; exact for k <= 20, log-binned (geometric centre) above."""
    v = np.asarray(values, dtype=np.int64)
    c = np.ones_like(v) if counts is None else np.asarray(counts, dtype=np.int64)
    keep = c > 0
    v, c = v[keep], c[keep]
    out = [[int(k), int(c[v == k].sum())] for k in np.unique(v[v <= 20])]
    big = v > 20
    if big.any():
        edges = np.unique(np.floor(np.logspace(np.log10(21), np.log10(v.max() + 1),
                                               max(2, int(np.log10(v.max() / 21 + 1)
                                                          * per_decade) + 2))).astype(int))
        idx = np.digitize(v[big], edges)
        for b in np.unique(idx):
            lo = edges[b - 1] if b > 0 else 21
            hi = edges[b] if b < len(edges) else v.max() + 1
            out.append([int(round(np.sqrt(lo * max(hi - 1, lo)))), int(c[big][idx == b].sum())])
    return out


def wilson(k, n, z=2.5758):
    """99 % Wilson interval for a binomial proportion."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(c - h), float(c + h)]


def bh(logp, alpha, m=None):
    """Benjamini-Hochberg on natural-log p-values; m = number of tests (>= len)."""
    m = len(logp) if m is None else m
    o = np.argsort(logp)
    thr = np.log(alpha * np.arange(1, len(logp) + 1) / m)
    ok = np.nonzero(logp[o] <= thr)[0]
    keep = np.zeros(len(logp), bool)
    if ok.size:
        keep[o[:ok[-1] + 1]] = True
    return keep


def _tail(logpmf_c, ratio, c, terms):
    """log P(X >= c) = log pmf(c) + log sum_t prod ratio, remainder bounded geometrically.

    Valid for any log-concave pmf: once the step ratio r(x) = pmf(x+1)/pmf(x) is < 1 it
    keeps falling, so the remainder after `terms` steps is <= last / (1 - r). The result
    is therefore an upper bound on p (conservative) and exact to ~1e-12 relative for the
    significant pairs, where r is far below 1. Where c + terms is still below the mode
    the tail is O(1), and log p = 0 is returned (also conservative)."""
    acc = np.ones_like(logpmf_c)
    cur = np.ones_like(logpmf_c)
    x = c.astype(np.float64)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        for _ in range(terms):
            r = ratio(x)
            cur = cur * r
            acc += cur
            x = x + 1
        r = ratio(x)
        rem = np.where(r < 1, cur * r / np.maximum(1 - r, 1e-300), np.inf)
        out = logpmf_c + np.log(acc + rem)
    out = np.where(np.isnan(out), 0.0, out)
    return np.minimum(out, 0.0)


def hypergeom_logsf(c, Ntot, ka, kb, terms=64, chunk=1 << 18):
    """log P(X >= c), X ~ Hypergeom(Ntot, ka, kb). scipy's is ~1 ms/pair at Ntot = 7e7."""
    c, ka, kb = (np.asarray(a, dtype=np.float64) for a in (c, ka, kb))
    out = np.empty(len(c))
    for s in range(0, len(c), chunk):
        cc, a, b = c[s:s + chunk], ka[s:s + chunk], kb[s:s + chunk]
        lp = (special.gammaln(a + 1) - special.gammaln(cc + 1) - special.gammaln(a - cc + 1)
              + special.gammaln(Ntot - a + 1) - special.gammaln(b - cc + 1)
              - special.gammaln(Ntot - a - b + cc + 1)
              - special.gammaln(Ntot + 1) + special.gammaln(b + 1)
              + special.gammaln(Ntot - b + 1))
        ratio = lambda x: np.clip((a - x) * (b - x) / ((x + 1) * (Ntot - a - b + x + 1)), 0, None)
        out[s:s + chunk] = _tail(lp, ratio, cc, terms)
    return out


def poisson_logsf(c, lam, terms=64):
    """log P(X >= c), X ~ Poisson(lam); exact (scipy) where it is representable."""
    c, lam = np.asarray(c, np.float64), np.asarray(lam, np.float64)
    with np.errstate(divide="ignore"):
        ex = np.log(special.pdtrc(c - 1, lam))
    bad = ~np.isfinite(ex)
    if bad.any():
        cb, lb = c[bad], lam[bad]
        lp = cb * np.log(lb) - lb - special.gammaln(cb + 1)
        ex[bad] = _tail(lp, lambda x: lb / (x + 1), cb, terms)
    return np.where(c <= 0, 0.0, ex)


# ---- BiCM ---------------------------------------------------------------------------
def _cols_given_rows(a, n, k, b0=None, it=200):
    """For fixed row log-fitnesses a (one per degree class, multiplicity n), the column
    equations sum_L n_L expit(a_L + b_c) = k_c are independent and monotone in b_c:
    safeguarded Newton, vectorised over columns."""
    lo = np.full(len(k), -80.0); hi = np.full(len(k), 80.0)
    b = np.zeros(len(k)) if b0 is None else b0.copy()
    for _ in range(it):
        p = special.expit(a[:, None] + b[None, :])
        g = n @ p - k
        hi = np.where(g > 0, b, hi); lo = np.where(g <= 0, b, lo)
        dg = n @ (p * (1 - p))
        nb = b - g / np.maximum(dg, 1e-300)
        out = (nb <= lo) | (nb >= hi) | ~np.isfinite(nb)
        nb = np.where(out, (lo + hi) / 2, nb)
        if np.max(np.abs(nb - b)) < 1e-13:
            return nb
        b = nb
    return b


def bicm_fit(row_deg, row_mult, col_deg, tol=1e-10):
    """Bipartite configuration model with rows grouped by degree.

    p_rc = x_r y_c / (1 + x_r y_c),  sum_c p_rc = d_r,  sum_r p_rc = k_c.
    Rows with equal degree have equal x (the solution is unique), so only the
    degree histogram is needed: G degree classes instead of 73M rows. The column
    log-fitnesses are solved exactly given the rows; the G-1 free row parameters
    (one is a gauge) by scipy root. Returns log x (per class), log y, and diagnostics."""
    d = np.asarray(row_deg, np.float64); n = np.asarray(row_mult, np.float64)
    k = np.asarray(col_deg, np.float64)
    assert np.isclose((d * n).sum(), k.sum()), "row and column degree sums differ"
    E = k.sum()
    # initial guess: Chung-Lu, x_r y_c ~ d_r k_c / E
    a0 = np.log(d / np.sqrt(E)); b = np.log(k / np.sqrt(E))
    g0 = int(np.argmax(n))                         # fix the gauge on the largest class
    state = {"b": b}

    def resid(af):
        a = np.insert(af, g0, a0[g0])
        state["b"] = _cols_given_rows(a, n, k, state["b"])
        p = special.expit(a[:, None] + state["b"][None, :])
        r = (p.sum(1) - d) / d
        return np.delete(r, g0)

    sol = optimize.root(resid, np.delete(a0, g0), method="hybr", tol=tol)
    a = np.insert(sol.x, g0, a0[g0])
    b = _cols_given_rows(a, n, k, state["b"])
    p = special.expit(a[:, None] + b[None, :])
    err_r = np.max(np.abs(p.sum(1) - d) / d)
    err_c = np.max(np.abs(n @ p - k) / k)
    return {"a": a, "b": b, "P": p, "err_rows": float(err_r), "err_cols": float(err_c),
            "converged": bool(sol.success), "nfev": int(sol.nfev)}


def vmotif_moments(P, n, i, j, chunk=1 << 17):
    """Mean, variance and third central moment of V_ij = sum_r B_ri B_rj under BiCM, for
    the listed pairs. V is Poisson-binomial with q_r = p_ri p_rj, n_L rows per class."""
    out = np.zeros((3, len(i)))
    for s in range(0, len(i), chunk):
        q = P[:, i[s:s + chunk]] * P[:, j[s:s + chunk]]
        out[0, s:s + chunk] = n @ q
        out[1, s:s + chunk] = n @ (q * (1 - q))
        out[2, s:s + chunk] = n @ (q * (1 - q) * (1 - 2 * q))
    return out


def rna_logsf(c, mu, var, m3):
    """Refined normal approximation (Volkova 1996) to the Poisson-binomial upper tail."""
    sd = np.sqrt(var)
    x = (c - 0.5 - mu) / sd
    g = m3 / sd ** 3
    sf = special.ndtr(-x) - g * (1 - x * x) * np.exp(-x * x / 2) / (6 * np.sqrt(2 * np.pi))
    with np.errstate(divide="ignore"):
        return np.log(np.clip(sf, 1e-320, 1.0))


# ---- compound-mode projection of a partition ----------------------------------------
def clique_stats(sizes, isolated=0):
    """Exact statistics of the one-mode projection of a single-label affiliation: one
    clique per label of size n_c, plus `isolated` unlabelled nodes. No graph is built."""
    s = np.asarray(sizes, dtype=np.float64)
    s = s[s > 0]
    n = s.sum() + isolated
    e = (s * (s - 1) / 2).sum()
    deg, cnt = np.unique(np.r_[s - 1, np.zeros(1 if isolated else 0)], return_inverse=True)
    w = np.bincount(cnt, weights=np.r_[s, [isolated] if isolated else []])
    ends = s[s >= 2]                       # both endpoints of every edge have degree n_c - 1
    assort = 1.0 if len(np.unique(ends)) > 1 else float("nan")
    return {"nodes": int(n), "edges": float(e), "density": float(e / (n * (n - 1) / 2)),
            "mean_degree": float(2 * e / n),
            "degree_hist": [[int(d), int(c)] for d, c in zip(deg, w)],
            "avg_clustering": float(s[s >= 3].sum() / n),
            "transitivity": 1.0 if (s >= 3).any() else 0.0,
            "components": int(len(s) + isolated), "giant_share": float(s.max() / n),
            "isolated": int(isolated), "assortativity": assort}


# ---- graphs -------------------------------------------------------------------------
def adjacency(u, v, n=N):
    A = sparse.coo_matrix((np.ones(len(u)), (u, v)), shape=(n, n)).tocsr()
    A = ((A + A.T) > 0).astype(np.float64)
    A.setdiag(0); A.eliminate_zeros()
    return A.tocsr()


def clustering(A):
    """Average clustering (networkx convention: degree < 2 counts as 0, over nodes with
    degree >= 1) and transitivity, from a symmetric 0/1 sparse adjacency."""
    k = np.asarray(A.sum(1)).ravel()
    tri = np.asarray((A @ A).multiply(A).sum(1)).ravel() / 2
    node = k > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.where(k > 1, 2 * tri / (k * (k - 1)), 0.0)
    trip = (k * (k - 1) / 2).sum()
    return float(c[node].mean()) if node.any() else 0.0, float(tri.sum() / trip) if trip else 0.0


def weighted_clustering_dense(Wd):
    """Onnela et al. (2005) weighted clustering, networkx's definition, on a dense matrix:
    c_i = sum_jk (w_ij w_jk w_ki)^(1/3) / (k_i (k_i - 1)), weights scaled by max(w)."""
    X = np.cbrt(Wd / Wd.max())
    k = (Wd > 0).sum(1).astype(np.float64)
    num = ((X @ X) * X).sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.where(k > 1, num / (k * (k - 1)), 0.0)
    return float(c[k > 0].mean())


def bimod(i, j, c, k, labels, deg_hist):
    """Bipartite modularity of a partition of the category layer (Guimera, Sales-Pardo &
    Amaral 2007, PRE 76:036102, Eq. 5), from the projection and degrees alone:

        M_B = sum_s [ sum_{a!=b in s} C_ab / sum_d n_d d(d-1)  -  sum_{a!=b in s} k_a k_b / (sum_d n_d d)^2 ]

    sums over ordered pairs; (i, j, c) the strictly-upper off-diagonal C; labels < 0 are
    singleton modules (contribute 0). The null is bipartite: no projection clique is
    counted as 'within-module' structure, unlike Newman-Girvan Q on C."""
    L = np.arange(len(deg_hist), dtype=np.float64)
    n = np.asarray(deg_hist, np.float64)
    den1, den2 = (n * L * (L - 1)).sum(), (n * L).sum() ** 2
    lab = np.asarray(labels)
    same = (lab[i] == lab[j]) & (lab[i] >= 0)
    t1 = 2 * np.asarray(c, np.float64)[same].sum() / den1
    kk = np.asarray(k, np.float64)
    m = lab >= 0
    _, inv = np.unique(lab[m], return_inverse=True)
    K = np.bincount(inv, weights=kk[m]); K2 = np.bincount(inv, weights=kk[m] ** 2)
    return float(t1 - (K ** 2 - K2).sum() / den2)


def double_edge_swap(u, v, nswap, rng):
    """Degree-preserving rewiring, the algorithm of nx.double_edge_swap ((a,b),(c,d) ->
    (a,d),(c,b), rejecting self-loops and multi-edges), on int arrays and one key set:
    ~10x faster than networkx at 1e5-1e6 edges. Returns rewired (u, v), accepted swaps."""
    u, v = u.astype(np.int64).copy(), v.astype(np.int64).copy()
    key = lambda a, b: a * N + b if a < b else b * N + a
    S = {key(a, b) for a, b in zip(u.tolist(), v.tolist())}
    E = len(u)
    tries = 0; done = 0
    while done < nswap and tries < 20 * nswap:
        i, j, flip = rng.integers(0, E, 4096), rng.integers(0, E, 4096), rng.random(4096) < 0.5
        for x, y, f in zip(i.tolist(), j.tolist(), flip.tolist()):
            tries += 1
            a, b = int(u[x]), int(v[x])
            c, d = (int(v[y]), int(u[y])) if f else (int(u[y]), int(v[y]))
            if x == y or a == d or c == b or a == c or b == d:
                continue
            k1, k2 = key(a, d), key(c, b)
            if k1 in S or k2 in S:
                continue
            S.discard(key(a, b)); S.discard(key(c, d)); S.add(k1); S.add(k2)
            v[x] = d; u[y], v[y] = c, b
            done += 1
            if done >= nswap:
                break
    return u, v, done


def _contingency(a, b):
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    M = sparse.coo_matrix((np.ones(len(a)), (ai, bi))).tocsr()
    return M


def nmi(a, b):
    """Normalised mutual information, arithmetic-mean normalisation (sklearn default)."""
    M = _contingency(a, b); n = M.sum()
    pa = np.asarray(M.sum(1)).ravel() / n; pb = np.asarray(M.sum(0)).ravel() / n
    m = M.tocoo(); pij = m.data / n
    I = float((pij * np.log(pij / (pa[m.row] * pb[m.col]))).sum())
    ha = -float((pa * np.log(pa)).sum()); hb = -float((pb * np.log(pb)).sum())
    return 1.0 if ha + hb == 0 else 2 * I / (ha + hb)


def ari(a, b):
    """Adjusted Rand index (Hubert & Arabie 1985)."""
    M = _contingency(a, b)
    c2 = lambda x: x * (x - 1) / 2
    sij = c2(M.data).sum()
    sa = c2(np.asarray(M.sum(1)).ravel()).sum(); sb = c2(np.asarray(M.sum(0)).ravel()).sum()
    exp = sa * sb / c2(M.sum())
    den = (sa + sb) / 2 - exp
    return 1.0 if den == 0 else float((sij - exp) / den)


def jaccard(e1, e2):
    s1, s2 = set(e1), set(e2)
    return len(s1 & s2) / len(s1 | s2) if s1 | s2 else float("nan")


def edge_keys(u, v):
    u, v = np.minimum(u, v), np.maximum(u, v)
    return (u.astype(np.int64) * N + v).tolist()


# ---- nestedness and economic complexity --------------------------------------------
def nodf(A):
    """NODF (Almeida-Neto et al. 2008), 0-100, order-free form: for every pair of rows
    (and of columns) with k_i > k_j, the share of j's ones that i also has; equal-degree
    pairs score 0; averaged over all C(n,2) + C(m,2) pairs. Empty rows/cols dropped."""
    A = np.asarray(A, dtype=np.float32)
    A = A[A.sum(1) > 0][:, A.sum(0) > 0]

    def side(X):
        k = X.sum(1)
        O = X @ X.T
        num = np.where(k[:, None] > k[None, :], O / np.maximum(k[None, :], 1), 0.0)
        n = len(k)
        return float(num.sum()), n * (n - 1) / 2
    r, nr = side(A); c, nc = side(A.T)
    return 100 * (r + c) / (nr + nc), 100 * r / nr, 100 * c / nc


def curveball(A, rng, trades=None):
    """Fixed-fixed null (Strona et al. 2014): random pairwise trades of non-shared
    presences between rows preserve every row and column total."""
    rows = [set(np.nonzero(r)[0]) for r in np.asarray(A)]
    nr = len(rows)
    for _ in range(trades or 5 * nr):
        i, j = rng.choice(nr, 2, replace=False)
        a, b = rows[i] - rows[j], rows[j] - rows[i]
        if not a or not b:
            continue
        shared = rows[i] & rows[j]
        pool = list(a | b)
        rng.shuffle(pool)
        rows[i] = shared | set(pool[:len(a)])
        rows[j] = shared | set(pool[len(a):])
    out = np.zeros_like(np.asarray(A))
    for r, s in enumerate(rows):
        out[r, list(s)] = 1
    return out


def fitness_complexity(M, it=200):
    """Tacchella et al. (2012): F_c = sum_p M_cp Q_p,  Q_p = 1 / sum_c M_cp / F_c,
    each normalised to mean 1 per iteration. Returns F (rows), Q (cols)."""
    M = np.asarray(M, np.float64)
    F = np.ones(M.shape[0]); Q = np.ones(M.shape[1])
    for _ in range(it):
        Fn = M @ Q
        Qn = 1.0 / np.maximum(M.T @ (1.0 / np.maximum(F, 1e-300)), 1e-300)
        F, Q = Fn / Fn.mean(), Qn / Qn.mean()
    return F, Q


# ---- output -------------------------------------------------------------------------
def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return None if not np.isfinite(f) else f
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def save_section(key, payload, ui=None):
    """Merge one script's section (and its UI fragment) into bipartite_results.json,
    then rebuild ui.json from whatever sections exist."""
    res = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    res[key] = _clean(payload)
    res.setdefault("ui_parts", {})[key] = _clean(ui or {})
    RESULTS.write_text(json.dumps(res, indent=1))
    build_ui(res)


def build_ui(res):
    s = scan_summary(require_complete=False)
    parts = res.get("ui_parts", {})
    ui = {"meta": {"generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                   "rows": s["rows"],
                   "classes": s["label_bipartite"]["class_nodes_nonzero"],
                   "substituents": s["substituent_bipartite"]["substituent_nodes_nonzero"],
                   "incidences": s["label_bipartite"]["edges"],
                   "substituent_incidences": s["substituent_bipartite"]["edges"]},
          "headline": [], "modes": [], "projections": [], "backbones": [], "edges": [],
          "nodes": {}, "communities": [], "fg": {},
          "charts": {"density_by_level": [], "hub_removal": [], "signature_sizes": []}}
    for key in sorted(parts):
        p = parts[key]
        for f in ("headline", "modes", "projections", "backbones"):
            ui[f] += p.get(f, [])
        if p.get("edges"):
            ui["edges"] = p["edges"]
        if p.get("communities"):
            ui["communities"] = p["communities"]
        if p.get("fg"):
            ui["fg"] = p["fg"]
        for nid, v in p.get("nodes", {}).items():
            ui["nodes"].setdefault(nid, {"k": 0, "comm": -1, "fit": None}).update(v)
        for c, v in p.get("charts", {}).items():
            ui["charts"][c] = v
        for f, v in p.get("extra", {}).items():
            ui[f] = v
    (OUT / "ui.json").write_text(json.dumps(_clean(ui), separators=(",", ":")))
    return ui
