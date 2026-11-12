#!/usr/bin/env python3
"""
Extension — one parallel pass that turns the release into bipartite incidence statistics.

The release is an affiliation (two-mode) network: 73.1M compounds on one side, ChemOnt
categories on the other. Two incidence relations are read from every row:

    B   compound -> label set   M(r) = {direct parent} U {alternative parents}
    S   compound -> substituent set (the release's functional-group annotations)

M(r) is the paper's non-redundant label set: by R5 no member is an ancestor of another,
so every ancestor is implied and counting the closure would only re-count the tree.

What is accumulated, all bounded by the category count (N = 4832 ids) not the row count:

    k_M[c]          compounds carrying category c in M          (class-mode degree)
    k_S[s]          compounds carrying substituent s
    deg_M[L]        compounds with |M| = L                      (compound-mode degree)
    deg_S[L]        compounds with |S| = L
    C = B^T B       class-mode one-mode projection, upper triangle incl. diagonal
    W = B^T diag(1/(|M|-1)) B   Newman's (2001) collaboration-weighted projection
    D = P^T S       direct parent x substituent (P = one-hot direct parent)
    F = S^T S       substituent-mode one-mode projection, upper triangle
    sig             64-bit hash of each compound's sorted M; compounds with equal M are
                    structurally equivalent in the compound-mode projection

The compound-mode projection is never materialised (Sum_c C(k_c, 2) ~ 1e15 edges); its
statistics are closed-form in k (scripts/14_bipartite_projections.py).

Rows are decompressed by the `zstd` CLI in this process and parsed by a worker pool in
byte chunks with bounded look-ahead, so memory stays at a few chunks whatever the input.
"""

import argparse
import json
import subprocess
import sys
import time
from collections import deque
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import paths  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "results" / "bipartite"
N = 4832                       # dense id space: release numeric ids are 0..4831
EXPECTED_ROWS = 73105281
CHUNK = 48 << 20               # bytes of decompressed text per task
HMAX = 1024                    # label-count histogram length; some compounds exceed 64
AP, SB = b'"alternative_parents"', b'"substituents"'


def field(blob, key):
    i = blob.find(key)
    if i < 0:
        return ()
    s = blob.find(b"[", i)
    e = blob.find(b"]", s)
    body = blob[s + 1:e]
    if not body:
        return ()
    return [int(t) for t in body.split(b",") if t != b"null"]


def upper(mat):
    """COO triplets of the upper triangle (diagonal included) of a symmetric product."""
    m = sparse.triu(mat, format="coo")
    return m.row.astype(np.int32), m.col.astype(np.int32), m.data.astype(np.int64)


def work(chunk):
    m_idx, m_ptr, s_idx, s_ptr, dps, sig = [], [0], [], [0], [], []
    rows = n_alt = dp_in_alts = no_dp = 0
    for line in chunk.split(b"\n"):
        if not line or line.startswith(b"inchikey\t"):
            continue
        c = line.split(b"\t")
        if len(c) < 6:
            continue
        rows += 1
        tree = c[4][1:-1].split(b",")
        dp = None if len(tree) != 5 or tree[4] == b"null" else int(tree[4])
        alts = field(c[5], AP)
        n_alt += len(alts)
        m = set(alts)
        if dp is None:
            no_dp += 1
        else:
            if dp in m:
                dp_in_alts += 1
            m.add(dp)
        m_idx.extend(m); m_ptr.append(len(m_idx))
        sig.append(hash(tuple(sorted(m))))      # int tuples hash deterministically
        subs = set(field(c[5], SB))
        s_idx.extend(subs); s_ptr.append(len(s_idx))
        dps.append(-1 if dp is None else dp)

    def csr(idx, ptr):
        idx = np.asarray(idx, dtype=np.int32)
        return sparse.csr_matrix((np.ones(len(idx), dtype=np.int32), idx,
                                  np.asarray(ptr, dtype=np.int64)), shape=(rows, N))

    B, S = csr(m_idx, m_ptr), csr(s_idx, s_ptr)
    dps = np.asarray(dps, dtype=np.int64)
    has = dps >= 0
    P = sparse.csr_matrix((np.ones(int(has.sum()), dtype=np.int32),
                           (np.nonzero(has)[0], dps[has])), shape=(rows, N))
    lm, ls = np.diff(B.indptr), np.diff(S.indptr)
    inv = np.zeros(rows); np.divide(1.0, lm - 1, out=inv, where=lm > 1)
    Bw = sparse.diags(inv) @ B
    return {
        "sig": np.asarray(sig, dtype=np.int64),
        "W": (lambda m: (m.row.astype(np.int32), m.col.astype(np.int32),
                         m.data.astype(np.float64)))(sparse.triu(Bw.T @ B, k=1,
                                                                 format="coo")),
        "rows": rows, "n_alt": n_alt, "dp_in_alts": dp_in_alts, "no_dp": no_dp,
        "kM": np.bincount(B.indices, minlength=N), "kS": np.bincount(S.indices, minlength=N),
        "degM": np.bincount(lm, minlength=HMAX), "degS": np.bincount(ls, minlength=HMAX),
        "C": upper(B.T @ B), "F": upper(S.T @ S),
        "D": (lambda m: (m.row.astype(np.int32), m.col.astype(np.int32),
                         m.data.astype(np.int64)))((P.T @ S).tocoo()),
    }


def chunks(stream):
    tail = b""
    while True:
        buf = stream.read(CHUNK)
        if not buf:
            if tail:
                yield tail
            return
        buf = tail + buf
        cut = buf.rfind(b"\n")
        if cut < 0:
            tail = buf
            continue
        yield buf[:cut]
        tail = buf[cut + 1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(paths.find("classyfire_main.tsv.zst") or ""))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit-chunks", type=int, default=0, help="smoke test: stop early")
    a = ap.parse_args()
    if not a.input:
        paths.require("classyfire_main.tsv.zst")

    C = np.zeros((N, N), dtype=np.int64)
    D = np.zeros((N, N), dtype=np.int64)
    F = np.zeros((N, N), dtype=np.int64)
    W = np.zeros((N, N), dtype=np.float64)
    sigs = []
    kM = np.zeros(N, dtype=np.int64); kS = np.zeros(N, dtype=np.int64)
    degM = np.zeros(HMAX, dtype=np.int64); degS = np.zeros(HMAX, dtype=np.int64)
    tot = {"rows": 0, "n_alt": 0, "dp_in_alts": 0, "no_dp": 0}

    def absorb(r):
        for k in tot:
            tot[k] += r[k]
        kM[:] += r["kM"]; kS[:] += r["kS"]
        if len(r["degM"]) > HMAX or len(r["degS"]) > HMAX:
            raise SystemExit(f"a compound has >= {HMAX} labels; raise HMAX")
        degM[:len(r["degM"])] += r["degM"]; degS[:len(r["degS"])] += r["degS"]
        sigs.append(r["sig"])
        for mat, key in ((C, "C"), (D, "D"), (F, "F"), (W, "W")):
            i, j, v = r[key]
            mat[i, j] += v                      # triplets are unique within one product

    proc = subprocess.Popen(["zstd", "-dcq", a.input], stdout=subprocess.PIPE,
                            bufsize=1 << 22)
    t0 = time.time()
    pending = deque()
    with Pool(a.workers) as pool:
        for k, ch in enumerate(chunks(proc.stdout)):
            if a.limit_chunks and k >= a.limit_chunks:
                break
            pending.append(pool.apply_async(work, (ch,)))
            while len(pending) > 2 * a.workers:
                absorb(pending.popleft().get())
            if k % 20 == 0:
                el = time.time() - t0
                sys.stderr.write(f"\r{tot['rows']/1e6:6.2f}M rows  {el:5.0f}s  "
                                 f"{tot['rows']/max(el,1e-9)/1e3:6.0f}k rows/s")
                sys.stderr.flush()
        while pending:
            absorb(pending.popleft().get())
    proc.stdout.close(); proc.kill() if a.limit_chunks else proc.wait()
    el = time.time() - t0
    sys.stderr.write("\n")

    OUT.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(OUT / "class_projection.npz", sparse.csr_matrix(C))
    sparse.save_npz(OUT / "dp_substituent.npz", sparse.csr_matrix(D))
    sparse.save_npz(OUT / "substituent_projection.npz", sparse.csr_matrix(F))
    sparse.save_npz(OUT / "class_projection_newman.npz", sparse.csr_matrix(W))
    _, mult = np.unique(np.concatenate(sigs), return_counts=True)
    del sigs
    sig_sizes, sig_count = np.unique(mult, return_counts=True)
    np.savez(OUT / "degrees.npz", kM=kM, kS=kS, degM=degM, degS=degS,
             sig_sizes=sig_sizes, sig_count=sig_count)

    edges_M = int(kM.sum()); edges_S = int(kS.sum())
    L = np.arange(HMAX)
    assert int(degM.sum()) == tot["rows"] == int(degS.sum()), "degree histogram lost rows"
    assert int((L * degM).sum()) == edges_M and int((L * degS).sum()) == edges_S
    assert int(np.triu(C, 1).sum()) == int((L * (L - 1) // 2 * degM).sum())
    degM, degS = np.trim_zeros(degM, "b"), np.trim_zeros(degS, "b")
    summary = {
        "input": a.input, "elapsed_s": round(el, 1), "rows": tot["rows"],
        "complete": tot["rows"] == EXPECTED_ROWS,
        "rows_without_direct_parent": tot["no_dp"],
        "alternative_parent_labels": tot["n_alt"],
        "rows_where_direct_parent_is_also_alt": tot["dp_in_alts"],
        "label_bipartite": {"definition": "M = {direct parent} U alternative parents",
                            "edges": edges_M, "class_nodes_nonzero": int((kM > 0).sum()),
                            "mean_compound_degree": edges_M / tot["rows"],
                            "projection_pairs_nonzero_offdiag":
                                int(np.count_nonzero(np.triu(C, 1)))},
        "substituent_bipartite": {"edges": edges_S,
                                  "substituent_nodes_nonzero": int((kS > 0).sum()),
                                  "mean_compound_degree": edges_S / tot["rows"],
                                  "projection_pairs_nonzero_offdiag":
                                      int(np.count_nonzero(np.triu(F, 1)))},
        "label_signatures": {"distinct": int(mult.size),
                             "singletons": int((mult == 1).sum()),
                             "largest": int(mult.max())},
        "dp_substituent_pairs_nonzero": int(np.count_nonzero(D)),
        "compound_degree_hist_M": degM.tolist(),
        "compound_degree_hist_S": degS.tolist(),
    }
    (OUT / "bipartite_scan.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if not k.startswith("compound_")},
                     indent=2))


if __name__ == "__main__":
    main()
