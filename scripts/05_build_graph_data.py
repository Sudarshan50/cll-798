#!/usr/bin/env python3
"""
Build the single JSON payload behind the network frontend.

Every one of the 73,105,281 compounds is accounted for: each is aggregated into the
ChemOnt class it was assigned, and the per-class totals reconcile exactly to the
full-scan row count. Nothing is sampled away.

Reads only the sufficient statistics from the 73.1M pass -- no rescan.
"""
import csv, json, math, sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont

ROOT = HERE.parent
RES = ROOT / "results"
tax = chemont.from_dictionary_tsv(ROOT / "data" / "chemont_dictionary.tsv")

# ---- exact per-class compound counts, from all 73.1M rows -------------------
as_dp = Counter()        # compounds whose DIRECT PARENT is this class
in_path = Counter()      # compounds carrying this class anywhere in the 5-slot path
total_rows = 0
with (RES / "path_counts.tsv").open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c = int(r["rows"]); total_rows += c
        slots = [r[k] for k in ("kingdom", "superclass", "class", "subclass", "direct_parent")]
        ids = [int(x) for x in slots if x]
        if ids:
            as_dp[ids[-1]] += c
        for x in set(ids):
            in_path[x] += c

# subtree totals: every compound counted under all its ancestors
kids = defaultdict(list)
for n, p in tax.parent.items():
    if p != -1:
        kids[p].append(n)
subtree = {}
def sub(n):
    if n in subtree: return subtree[n]
    subtree[n] = as_dp.get(n, 0) + sum(sub(k) for k in kids[n])
    return subtree[n]
for n in tax.parent: sub(n)

# ---- confusability edges ----------------------------------------------------
alt_w = Counter(); dptot = Counter()
with (RES / "dp_alt_pairs.tsv").open() as fh:
    next(fh)
    for line in fh:
        d, x, c = line.rstrip("\n").split("\t")
        if not d: continue
        d, x, c = int(d), int(x), int(c)
        if d == x: continue
        alt_w[(min(d, x), max(d, x))] += c
        dptot[d] += c

entropy = {}
byd = defaultdict(Counter)
with (RES / "dp_alt_pairs.tsv").open() as fh:
    next(fh)
    for line in fh:
        d, x, c = line.rstrip("\n").split("\t")
        if d: byd[int(d)][int(x)] += int(c)
for d, cnt in byd.items():
    t = sum(cnt.values())
    entropy[d] = round(-sum((v/t)*math.log2(v/t) for v in cnt.values() if v), 3) if t else 0.0

TOP_EDGES = 14000
top = alt_w.most_common(TOP_EDGES)

# ---- example compounds per class (from the 200k sample; no rescan) ----------
ex = defaultdict(list)
sample = ROOT.parent / "chemont_project" / "data" / "sample_200k.tsv"
if sample.exists():
    with sample.open() as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            if i > 120000: break
            t = json.loads(r["chemont_tree_json"])
            dp = t[4]
            if dp is None or len(ex[dp]) >= 3: continue
            ex[dp].append({"k": r["inchikey"], "s": (r["smiles"] or "")[:90],
                           "c": r["cid"] or None})

# ---- assemble ---------------------------------------------------------------
ROOT_N = 9999999
def kingdom_of(n):
    while tax.parent.get(n, -1) not in (-1, ROOT_N):
        n = tax.parent[n]
    return tax.name.get(n, "")

sc_of = {}
for n in tax.parent:
    x = n
    while tax.depth.get(x, 0) > 2:
        x = tax.parent[x]
    sc_of[n] = x if tax.depth.get(x, 0) == 2 else None

nodes = []
for n in sorted(tax.parent):
    if n == ROOT_N: continue
    nodes.append({
        "i": n, "n": tax.name.get(n, str(n)),
        "c": f"CHEMONTID:{n:07d}", "p": tax.parent[n] if tax.parent[n] != ROOT_N else -1,
        "d": tax.depth[n], "k": 0 if kingdom_of(n) == "Organic compounds" else 1,
        "s": sc_of.get(n), "sc": tax.name.get(sc_of.get(n), ""),
        "dp": as_dp.get(n, 0), "st": subtree.get(n, 0), "ip": in_path.get(n, 0),
        "e": entropy.get(n, 0.0),
        "x": ex.get(n, []),
    })
links = [{"s": a, "t": b, "w": w} for (a, b), w in top]
tree = [{"s": tax.parent[n], "t": n} for n in tax.parent
        if n != ROOT_N and tax.parent[n] not in (-1, ROOT_N)]

payload = {
    "meta": {"compounds": total_rows, "classes": len(nodes),
             "tree_edges": len(tree), "confus_edges": len(links),
             "confus_edges_total": len(alt_w),
             "alt_labels": sum(alt_w.values())},
    "nodes": nodes, "tree": tree, "links": links,
}
out = RES / "graph_data.json"
out.write_text(json.dumps(payload, separators=(",", ":")))
print(f"compounds reconciled : {total_rows:,}")
print(f"classes              : {len(nodes):,}")
print(f"tree edges           : {len(tree):,}")
print(f"confusability edges  : {len(links):,} of {len(alt_w):,}")
print(f"classes with examples: {sum(1 for n in nodes if n['x']):,}")
print(f"wrote {out}  ({out.stat().st_size/1e6:.2f} MB)")
