#!/usr/bin/env python3
"""
Derive the Phase-2 S2 channel -- alternative-parent co-occurrence -- without a rescan.

scripts/09_scan_bipartite.py accumulates C[a,b] = #compounds whose M = {dp} U alts holds
both a and b. A compound contributes to C[a,b] in exactly one of three disjoint ways:
both are alternative parents, or a is its direct parent and b an alternative, or the
reverse. The direct parent is never also listed as an alternative (the scan asserts the
count is 0), so

    altalt[a,b] = C[a,b] - n(dp=a, alt=b) - n(dp=b, alt=a)

with n from results/dp_alt_pairs.tsv. Writes the file `02_scan_dataset.py --pairs`
would have written: results/alt_cooccurrence_pairs.tsv  (a, b, rows; a < b).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
BIP = RES / "bipartite"

scan = json.loads((BIP / "bipartite_scan.json").read_text())
if not scan["complete"] or scan["rows_where_direct_parent_is_also_alt"] != 0:
    sys.exit("bipartite scan incomplete, or dp-in-alts != 0: the subtraction is invalid")

C = sparse.load_npz(BIP / "class_projection.npz").tocsr().astype(np.int64)
N = C.shape[0]
d, a, n = [], [], []
with (RES / "dp_alt_pairs.tsv").open() as fh:
    next(fh)
    for line in fh:
        x, y, c = line.rstrip("\n").split("\t")
        if x:
            d.append(int(x)); a.append(int(y)); n.append(int(c))
d, a, n = np.array(d), np.array(a), np.array(n, dtype=np.int64)
lo, hi = np.minimum(d, a), np.maximum(d, a)
DA = sparse.csr_matrix((n, (lo, hi)), shape=(N, N))       # duplicates summed: both directions
A = sparse.triu(C, k=1).tocsr() - sparse.triu(DA, k=1).tocsr()
A.eliminate_zeros()
A = A.tocoo()
if (A.data < 0).any():
    sys.exit(f"{(A.data < 0).sum()} negative alt-alt counts: inputs disagree")

order = np.argsort(-A.data, kind="stable")
with (RES / "alt_cooccurrence_pairs.tsv").open("w") as fh:
    fh.write("a\tb\trows\n")
    for i in order:
        fh.write(f"{A.row[i]}\t{A.col[i]}\t{A.data[i]}\n")

scan_path = RES / "scan.json"
print(f"alt-alt pairs: {A.nnz:,}  label-pair instances: {int(A.data.sum()):,}")
print(f"wrote {RES / 'alt_cooccurrence_pairs.tsv'}  (scan.json left untouched: {scan_path})")
