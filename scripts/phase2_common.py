"""
Shared machinery for Phase 2 (scripts 10-13c).

Everything here works off the Phase-1 sufficient statistics (dp_alt_pairs.tsv,
path_counts.tsv, inode_sets.tsv) plus the frozen tree; nothing rescans the 2 GB file
except `extract_release_rows`, which pulls the ~820 evaluation compounds by CID for the
leakage check. Kept separate so the toy-taxonomy tests in test_phase2.py exercise the
same PMI / LCA / permutation code the analysis uses.
"""
import csv, hashlib, json, math, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont  # noqa: E402
import paths    # noqa: E402

ROOT = HERE.parent
RES = ROOT / "results"
OUT = RES / "phase2"
DICT = ROOT / "data" / "chemont_dictionary.tsv"
OBO = ROOT / "data" / "ChemOnt_2_1.obo"
SUP = ROOT / "data" / "paper_supplementary"
M4 = SUP / "13321_2016_174_MOESM4_ESM.xlsx"
M5 = SUP / "13321_2016_174_MOESM5_ESM.xlsx"
PREREG = ROOT / "docs" / "phase2_prereg.md"
ROOT_N = chemont.ROOT_NUM
K = 4832                                   # dense release id space 0..4831
N_COMPOUNDS = 73105281
LEVELS = ("kingdom", "superclass", "class", "subclass")
THETAS, TAUS, NMINS = (2, 3, 4, 5), (1, 2, 3), (10, 100, 1000)
THETA0 = 4
SEED = 798


def sha256(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def dump(obj, name):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=1, default=float))


# ---------------------------------------------------------------- the tree
def taxonomy(path=DICT):
    return chemont.from_dictionary_tsv(path)


def anc_matrix(tax, k=None):
    """A[n, L] = ancestor-or-self of n at level L (root = level 0), -1 below n's depth."""
    k = k or (max(x for x in tax.parent if x != ROOT_N) + 1)
    A = np.full((k, 12), -1, dtype=np.int64)
    for n, chain in tax.anc_at.items():
        if n != ROOT_N:
            for L, x in chain.items():
                A[n, L] = x
    return A


def depth_array(tax, k=None):
    k = k or (max(x for x in tax.parent if x != ROOT_N) + 1)
    d = np.full(k, -1, dtype=np.int64)
    for n, v in tax.depth.items():
        if n != ROOT_N:
            d[n] = v
    return d


def lca_depth(A, u, v):
    """Ancestor chains are prefixes, so the LCA depth is the length of the shared prefix."""
    au, av = A[u], A[v]
    return ((au == av) & (au >= 0)).sum(1) - 1


def parent_array(tax, k=None):
    k = k or (max(x for x in tax.parent if x != ROOT_N) + 1)
    p = np.full(k, -1, dtype=np.int64)
    for n, x in tax.parent.items():
        if n != ROOT_N:
            p[n] = x
    return p


def label_slots(A, depth, ids):
    """Padded 4-slot label (kingdom..subclass) as the release writes it: slots deeper than
    the class repeat the class itself (chemont.Taxonomy.expected_path)."""
    ids = np.asarray(ids)
    L = np.arange(1, 5)
    return np.where(L[None, :] <= depth[ids][:, None], A[ids][:, 1:5], ids[:, None])


def first_changed_level(A, depth, d, a):
    """0..3 = shallowest of kingdom/superclass/class/subclass whose slot differs, 4 = none."""
    diff = label_slots(A, depth, d) != label_slots(A, depth, a)
    return np.where(diff.any(1), diff.argmax(1), 4)


def superclass_label(tax):
    """SuperClass of each node; Kingdom-level nodes get their own kingdom."""
    out = {}
    for n, chain in tax.anc_at.items():
        if n != ROOT_N:
            out[n] = chain.get(2, chain.get(1))
    return out


# ---------------------------------------------------------------- S1
def load_s1(path=RES / "dp_alt_pairs.tsv"):
    """Directed (d, a, n). The 5 empty-direct-parent rows (35 labels) are dropped."""
    d, a, n = [], [], []
    with Path(path).open() as fh:
        next(fh)
        for line in fh:
            x, y, c = line.rstrip("\n").split("\t")
            if x:
                d.append(int(x)); a.append(int(y)); n.append(int(c))
    return np.array(d), np.array(a), np.array(n, dtype=np.float64)


def pmi_directed(d, a, n, k=K):
    """PMI(d,a) = log2[ P(a|d) / P(a) ], marginals from the rows themselves."""
    dtot = np.bincount(d, n, k)
    amarg = np.bincount(a, n, k)
    N = n.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        pmi = np.log2(n / dtot[d]) - np.log2(amarg[a] / N)
    return pmi, dtot, amarg


class Undirected:
    """Collapse directed (d, a) rows to undirected pairs {u < v}. `inv` maps each directed
    row to its pair so the bootstrap can re-collapse resampled counts cheaply."""

    def __init__(self, d, a, k=K):
        self.k = k
        key = np.minimum(d, a) * k + np.maximum(d, a)
        uniq, self.inv = np.unique(key, return_inverse=True)
        self.u, self.v = uniq // k, uniq % k
        self.fwd = d < a                     # row runs u -> v
        self.m = len(uniq)

    def collapse(self, n, pmi):
        f, b, inv = self.fwd, ~self.fwd, self.inv
        n_uv = np.bincount(inv[f], n[f], self.m)
        n_vu = np.bincount(inv[b], n[b], self.m)
        p_uv = np.full(self.m, np.nan); p_uv[inv[f]] = pmi[f]
        p_vu = np.full(self.m, np.nan); p_vu[inv[b]] = pmi[b]
        p_uv[n_uv == 0] = np.nan; p_vu[n_vu == 0] = np.nan
        return n_uv, n_vu, np.fmin(p_uv, p_vu)


def criterion_b(pmi_sym, w_raw, tau, nmin):
    return (pmi_sym >= tau) & (w_raw >= nmin)


def directed_mask(U, pair_mask):
    """Directed rows whose undirected pair is retained."""
    return pair_mask[U.inv]


def entropy_by_dp(d, n, mask, k=K):
    """H(d) over p(a|d) renormalised to retained rows. 0 with < 2 retained alternatives."""
    dm, nm = d[mask], n[mask]
    tot = np.bincount(dm, nm, k)
    p = nm / tot[dm]
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -np.bincount(dm, np.where(p > 0, p * np.log2(p), 0.0), k)
    return np.where(tot > 0, h, 0.0) + 0.0


def max_share_score(d, n, mask, k=K):
    """Alternative ambiguity score 1 - max_a p(a|d) over retained rows (0 if none)."""
    dm, nm = d[mask], n[mask]
    tot = np.bincount(dm, nm, k)
    mx = np.zeros(k)
    np.maximum.at(mx, dm, nm)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(tot > 0, 1 - mx / tot, 0.0)


# ---------------------------------------------------------------- sizes
def compounds_as_dp(path=RES / "path_counts.tsv", k=K):
    s = np.zeros(k)
    with Path(path).open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["direct_parent"]:
                s[int(r["direct_parent"])] += int(r["rows"])
    return s


# ---------------------------------------------------------------- S3 / calibration
def r5_pairs(tax, d, a):
    """Directed S1 rows where one end is a strict ancestor of the other (R5 violations)."""
    out = []
    for i, (x, y) in enumerate(zip(d.tolist(), a.tolist())):
        if x in tax.parent and y in tax.parent and (y in tax.anc[x] or x in tax.anc[y]):
            out.append(i)
    return np.array(out, dtype=np.int64)


def path_shape_pairs(tax, path=RES / "path_counts.tsv"):
    """PATH-violating label shapes and the (reported, expected) node pairs at their
    non-null mismatched slots. Returns (n_shapes, {pair: rows})."""
    shapes, pairs = set(), Counter()
    with Path(path).open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = [int(r[s]) if r[s] else None
                 for s in ("kingdom", "superclass", "class", "subclass", "direct_parent")]
            if tax.check_path(t) is False:
                shapes.add(tuple(t))
                exp = tax.expected_path(t[4])
                for rep, e in zip(t[:4], exp[:4]):
                    if rep is not None and rep != e:
                        pairs[(min(rep, e), max(rep, e))] += int(r["rows"])
    return len(shapes), pairs


def r6_pairs(tax, path=RES / "inode_sets.tsv"):
    """(intermediate node, direct parent) pairs violating R6 (inclusive reading)."""
    pairs = Counter()
    with Path(path).open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["direct_parent"]:
                continue
            sub = int(r["subclass"]) if r["subclass"] else None
            dp = int(r["direct_parent"])
            inodes = [int(x) for x in r["intermediate_nodes"].split(",") if x]
            for x in tax.check_r6(sub, dp, inodes) or []:
                if x != dp:
                    pairs[(min(x, dp), max(x, dp))] += int(r["rows"])
    return pairs


# ---------------------------------------------------------------- S2
def load_s2(d, a, n, k=K):
    """Alt-alt co-occurrence (a < b, rows). Uses results/alt_cooccurrence_pairs.tsv when
    present, else derives it from the label projection C[a,b] - n(a,b) - n(b,a).
    Returns (a, b, c, source) or None if neither input is complete."""
    tsv = RES / "alt_cooccurrence_pairs.tsv"
    if tsv.exists() and tsv.stat().st_size > 0:
        x, y, c = [], [], []
        with tsv.open() as fh:
            next(fh)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                x.append(int(p[0])); y.append(int(p[1])); c.append(int(p[2]))
        x, y, c = np.array(x), np.array(y), np.array(c, dtype=np.float64)
        keep = (x != y) & (c > 0)
        return x[keep], y[keep], c[keep], str(tsv.relative_to(ROOT))
    scan = RES / "bipartite" / "bipartite_scan.json"
    npz = RES / "bipartite" / "class_projection.npz"
    if not (scan.exists() and npz.exists() and json.loads(scan.read_text()).get("complete")):
        return None
    import scipy.sparse as sp
    C = sp.load_npz(npz).tocsr()
    C = sp.triu(C, k=1).tocoo()
    dp = Counter()
    for x, y, c in zip(d.tolist(), a.tolist(), n.tolist()):
        dp[(min(x, y), max(x, y))] += c
    cc = C.data.astype(np.float64) - np.array([dp.get((x, y), 0.0)
                                               for x, y in zip(C.row.tolist(), C.col.tolist())])
    keep = cc > 0
    return C.row[keep], C.col[keep], cc[keep], "derived from results/bipartite/class_projection.npz"


# ---------------------------------------------------------------- validation inputs
def flagged_classes(tax):
    """Every distinct FP/FN category name of Additional file 4 sheet 2, resolved
    exact -> case-insensitive -> OBO synonym. Returns one audit row per name."""
    import warnings
    import openpyxl
    warnings.filterwarnings("ignore", category=UserWarning)
    wb = openpyxl.load_workbook(M4, read_only=True)
    rows = list(wb.worksheets[1].iter_rows(min_row=2, values_only=True))
    wb.close()
    fp = Counter(r[3] for r in rows if r[3] not in (None, ""))
    fn = Counter(r[4] for r in rows if r[4] not in (None, ""))
    exact = {v: k for k, v in tax.name.items()}
    lower = {v.lower(): k for k, v in tax.name.items()}
    syn = obo_synonyms()
    out = []
    for nm in sorted(set(fp) | set(fn)):
        i, how = exact.get(nm), "exact"
        if i is None:
            i, how = lower.get(nm.lower()), "case-insensitive"
        if i is None:
            i, how = syn.get(nm.lower()), "obo-synonym"
        if i is None:
            how = "unresolved"
        out.append({"name": nm, "fp": fp.get(nm, 0), "fn": fn.get(nm, 0), "id": i,
                    "resolved_by": how, "dictionary_name": tax.name.get(i) if i is not None else None})
    return out


def obo_synonyms():
    num = {}
    with DICT.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            num[r["chemont_id"]] = int(r["numeric_id"])
    syn, cur = {}, None
    for line in OBO.open(encoding="utf-8"):
        if line.startswith("id: "):
            cur = num.get(line[4:].strip())
        elif line.startswith("synonym: ") and cur is not None:
            s = line.split('"')[1].lower()
            syn.setdefault(s, cur)
    return syn


def evaluation_cids():
    """CIDs of the 800 test compounds (Additional file 4 sheet 1) and of the 20 ChEBI
    comparison compounds (Additional file 5) -- two different sets of 'twenty'."""
    import warnings
    import openpyxl
    warnings.filterwarnings("ignore", category=UserWarning)
    wb = openpyxl.load_workbook(M4, read_only=True)
    test = [int(r[0]) for r in wb.worksheets[0].iter_rows(min_row=2, values_only=True) if r[0]]
    wb.close()
    wb = openpyxl.load_workbook(M5, read_only=True)
    chebi = [r for r in wb.worksheets[0].iter_rows(min_row=2, values_only=True)
             if isinstance(r[0], int)]
    wb.close()
    return test, chebi


def extract_release_rows(cids, out):
    """One streaming pass: the release rows whose cid field is in `cids`. The pattern
    '\\t<cid>\\t' can only hit the cid column (inchikey and zinc ids are not numeric)."""
    src = paths.require("classyfire_main.tsv.zst")
    pat = OUT / "_cid_patterns.txt"
    pat.write_text("".join(f"\t{c}\t\n" for c in sorted(set(cids))))
    cmd = f"zstd -dcq '{src}' | rg -F -f '{pat}' > '{out}'"
    subprocess.run(cmd, shell=True, check=False)
    pat.unlink()
    rows = []
    want = set(cids)
    for line in Path(out).open():
        p = line.rstrip("\n").split("\t")
        if len(p) >= 6 and p[1].isdigit() and int(p[1]) in want:
            rows.append(p)
    return rows


# ---------------------------------------------------------------- statistics
def deciles(x):
    """Rank-based deciles (ties broken by order, so bins are equal-sized)."""
    r = np.argsort(np.argsort(x, kind="stable"), kind="stable")
    return (10 * r) // len(x)


def perm_test(score, flagged, strata, R=10000, seed=SEED):
    """Stratified permutation: each control set takes, per stratum, as many non-flagged
    members as there are flagged ones, without replacement. Returns (T_obs, null, p)."""
    score, flagged, strata = map(np.asarray, (score, flagged, strata))
    rng = np.random.default_rng(seed)
    T = score[flagged].mean()
    fl = np.zeros(len(score), bool); fl[flagged] = True
    tot = np.zeros(R)
    for s, k in Counter(strata[flagged].tolist()).items():
        pool = np.flatnonzero((strata == s) & ~fl)
        pick = np.argsort(rng.random((R, len(pool))), axis=1)[:, :k]
        tot += score[pool][pick].sum(1)
    null = tot / len(flagged)
    return T, null, (1 + int((null >= T).sum())) / (1 + R)


def bootstrap_effect(score, flagged, strata, R=10000, seed=SEED):
    """95 % CI of mean(flagged) - mean(matched expectation), resampling the flagged."""
    score, flagged, strata = map(np.asarray, (score, flagged, strata))
    fl = np.zeros(len(score), bool); fl[flagged] = True
    cm = {s: score[(strata == s) & ~fl].mean() for s in set(strata[flagged].tolist())}
    diff = score[flagged] - np.array([cm[s] for s in strata[flagged]])
    rng = np.random.default_rng(seed)
    bs = diff[rng.integers(0, len(diff), (R, len(diff)))].mean(1)
    return float(diff.mean()), [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def nmi(x, y):
    """Normalised mutual information, arithmetic-mean normalisation."""
    x, y = np.unique(x, return_inverse=True)[1], np.unique(y, return_inverse=True)[1]
    n = len(x)
    c = np.zeros((x.max() + 1, y.max() + 1)); np.add.at(c, (x, y), 1)
    px, py, pxy = c.sum(1) / n, c.sum(0) / n, c / n
    nz = pxy > 0
    mi = (pxy[nz] * np.log(pxy[nz] / np.outer(px, py)[nz])).sum()
    hx, hy = -(px * np.log(px)).sum(), -(py * np.log(py)).sum()
    return float(2 * mi / (hx + hy)) if hx + hy > 0 else 1.0


def ari(x, y):
    x, y = np.unique(x, return_inverse=True)[1], np.unique(y, return_inverse=True)[1]
    c = np.zeros((x.max() + 1, y.max() + 1)); np.add.at(c, (x, y), 1)
    comb = lambda v: (v * (v - 1) / 2).sum()
    s, a, b, t = comb(c), comb(c.sum(1)), comb(c.sum(0)), comb(np.array([len(x)]))
    e = a * b / t
    return float((s - e) / ((a + b) / 2 - e)) if (a + b) / 2 != e else 1.0


def spearman(x, y):
    from scipy.stats import spearmanr
    r, p = spearmanr(x, y)
    return float(r), float(p)
