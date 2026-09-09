#!/usr/bin/env python3
"""
Phase 1b + 1c — one streaming pass over the release.

DESIGN NOTE (why this is structured the way it is)
--------------------------------------------------
The first version applied every taxonomy rule inside the row loop, using 4825-bit
integer ancestor masks. That ran at ~0.011 M rows/s (109 CPU-minutes without finishing).
Raw parsing alone runs at ~0.64 M rows/s, so essentially all of it was rule evaluation.

This version instead accumulates BOUNDED SUFFICIENT STATISTICS in the hot loop and
evaluates the paper's rules afterwards, over distinct keys only:

    label paths          ~2.9k distinct  (vs 73.1M rows)
    intermediate-node sets ~1k distinct
    (direct_parent, alternative_parent) pairs   sparse, <= 4824^2

Two consequences, both deliberate:
  * It is fast (minutes, not hours).
  * Re-analysis needs no second pass for the PATH, R6 and dp-vs-alt rules: those are
    recomputable from the dumped statistics under any rule reading or taxonomy version.
    R5's alt-vs-alt clause is the exception -- it needs per-row alternative-parent sets,
    which are not stored. Measured over 200,000 rows, zero rows violate that clause
    without also violating dp-vs-alt, so no reported figure depends on it.

Alternative-parent CO-OCCURRENCE pairs (the Phase-2 class-class network input) are the
one statistic that is not free: ~78 pairs/row x 73.1M rows. They are behind --pairs so a
Phase-1 run is not taxed by Phase-2 work.

1b  PROVENANCE (not paper replication): the Zenodo release-note figures. These describe a
    2026 InChIKey-deduplicated third-party aggregation, NOT the 77M population the 2016
    paper reports.
1c  INVARIANTS (paper replication):
      PATH  the 5-slot label must be a genuine root-path of its direct parent
            (Additional file 1: the tree is "is_a complete")
      R5    Results: alternative parents have no "ancestor-descendant relationship" with
            each other or the direct parent. Table S1 says the weaker "parent-child".
            Both readings reported.
      R6    intermediate nodes are descendants of the subclass and ascendants of the
            direct parent. Strict and inclusive readings both reported.
"""

import argparse
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont  # noqa: E402
import paths    # noqa: E402

ROOT = HERE.parent
DEFAULT_INPUT = paths.find("classyfire_main.tsv.zst") or Path("classyfire_main.tsv.zst")
DICT = ROOT / "data" / "chemont_dictionary.tsv"
OBO = ROOT / "data" / "ChemOnt_2_1.obo"
RESULTS = ROOT / "results"

ZENODO = {"rows": 73105281, "with_cid": 68693298, "with_zinc": 30187323,
          "with_both": 25775344, "unresolved_slot_rows": 119795,
          "classes_seen": 4120, "classes_total": 4824}
SLOTS = ("kingdom", "superclass", "class", "subclass", "direct_parent")
EXPECTED_ROWS = ZENODO["rows"]


def progress(total, quiet):
    """tqdm when available, otherwise an equivalent stderr bar. No hard dependency."""
    if quiet:
        class _N:
            def update(self, n=1): pass
            def close(self): pass
            def set_postfix_str(self, s): pass
        return _N()
    try:
        from tqdm import tqdm
        return tqdm(total=total, unit="row", unit_scale=True, dynamic_ncols=True,
                    desc="scanning", smoothing=0.05, mininterval=0.5)
    except ImportError:
        class _Bar:
            def __init__(self):
                self.n, self.t0, self.last, self.post = 0, time.time(), 0.0, ""
            def set_postfix_str(self, s): self.post = s
            def update(self, k=1):
                self.n += k
                now = time.time()
                if now - self.last < 0.5:
                    return
                self.last = now
                el = now - self.t0
                frac = self.n / total if total else 0
                eta = (el / frac - el) if frac > 0.001 else 0
                bar = "#" * int(30 * frac) + "." * (30 - int(30 * frac))
                sys.stderr.write(
                    f"\rscanning [{bar}] {100*frac:5.1f}%  {self.n/1e6:6.2f}M/"
                    f"{total/1e6:.1f}M  {self.n/el/1e6:4.2f}M/s  "
                    f"eta {eta/60:4.1f}m  {self.post}")
                sys.stderr.flush()
            def close(self):
                sys.stderr.write("\n"); sys.stderr.flush()
        return _Bar()


def open_stream(path):
    if str(path).endswith(".zst"):
        p = subprocess.Popen(["zstd", "-dcq", str(path)], stdout=subprocess.PIPE)
        return p.stdout, p
    return open(path, "rb"), None


def ints(b):
    """b'0,12,null,3,9' -> (0,12,None,3,9)"""
    if not b:
        return ()
    return tuple(None if t == b"null" else int(t) for t in b.split(b","))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--json-only", action="store_true")
    ap.add_argument("--pairs", action="store_true",
                    help="also accumulate alternative-parent co-occurrence pairs "
                         "(Phase-2 class-class network); materially slower")
    ap.add_argument("--out", default=str(RESULTS / "scan.json"))
    a = ap.parse_args()

    # ---------------- hot loop: parse and tally only ----------------
    n = n_cid = n_zinc = n_both = 0
    paths = Counter()          # tree 5-tuple            -> rows
    isets = Counter()          # (subclass, dp, inodes)  -> rows
    dp_alt = Counter()         # (dp, alt)               -> rows
    alt_pairs = Counter()      # (a,b) a<b               -> rows   (--pairs only)
    n_alt_total = 0

    stream, proc = open_stream(a.input)
    total = a.limit or EXPECTED_ROWS
    bar = progress(total, a.quiet or a.json_only)
    t0 = time.time()
    stream.readline()
    AP, IN = b'"alternative_parents"', b'"intermediate_nodes"'
    for raw in stream:
        n += 1
        if a.limit and n > a.limit:
            n -= 1
            break
        c = raw.split(b"\t")
        if len(c) < 6:
            continue
        if c[1]:
            n_cid += 1
            if c[2]:
                n_both += 1
        if c[2]:
            n_zinc += 1

        t = ints(c[4][1:-1])
        paths[t] += 1
        blob = c[5]
        i = blob.find(AP)
        if i >= 0:
            s = blob.find(b"[", i); e = blob.find(b"]", s)
            alts = ints(blob[s + 1:e])
        else:
            alts = ()
        j = blob.find(IN)
        if j >= 0:
            s = blob.find(b"[", j); e = blob.find(b"]", s)
            inodes = ints(blob[s + 1:e])
        else:
            inodes = ()

        dp = t[4] if len(t) == 5 else None
        isets[(t[3] if len(t) == 5 else None, dp, inodes)] += 1
        n_alt_total += len(alts)
        for x in alts:
            dp_alt[(dp, x)] += 1
        if a.pairs and len(alts) > 1:
            sa = sorted(alts)
            for ii, x in enumerate(sa):
                for y in sa[ii + 1:]:
                    alt_pairs[(x, y)] += 1
        if not (n & 0xFFFF):
            bar.update(0x10000)
    bar.close()
    stream.close()
    if proc:
        proc.wait()
    el = time.time() - t0

    # ---------------- post-pass: rules over distinct keys only ----------------
    tax = chemont.from_dictionary_tsv(DICT)
    obo_tax = chemont.from_obo(OBO, DICT) if OBO.exists() else None
    drift = set(chemont.edge_diff(tax, obo_tax)) if obo_tax else set()

    null_by_slot = Counter(); dp_depth = Counter(); seen = set()
    n_null_slot = n_path_bad = n_path_bad_full = drift_touched = 0
    path_shapes = Counter()
    for t, cnt in paths.items():
        if len(t) != 5:
            continue
        if any(v is None for v in t):
            n_null_slot += cnt
            for i, v in enumerate(t):
                if v is None:
                    null_by_slot[SLOTS[i]] += cnt
        dp = t[4]
        seen.update(v for v in t if v is not None)
        if dp is None:
            continue
        dp_depth[tax.depth.get(dp)] += cnt
        if tax.check_path(list(t)) is False:
            n_path_bad += cnt
            # rows with a complete, non-null path: the population the release's
            # "correct for every row" claim is actually about. Subtracting
            # unresolved_slot_rows in prose mixed two different populations, because
            # the 7 rows with a null direct_parent never reach check_path at all.
            if all(v is not None for v in t):
                n_path_bad_full += cnt
            path_shapes[(t, tuple(tax.expected_path(dp)))] += cnt
            if any(x in drift for x in t if x is not None):
                drift_touched += cnt

    # R5 over distinct (dp, alt) pairs
    n_r5_pairs = n_r5_pairs_direct = 0
    r5_pair_rows = r5_pair_rows_direct = 0
    r5_pairs_detail = Counter()
    for (dp, x), cnt in dp_alt.items():
        if dp is None or dp not in tax.parent or x not in tax.parent:
            continue
        seen.add(x)
        rel = None
        if tax.parent.get(x) == dp or tax.parent.get(dp) == x:
            rel = "direct"
        elif x in tax.anc[dp] or dp in tax.anc[x]:
            rel = "deeper"
        if rel:
            n_r5_pairs += 1
            r5_pair_rows += cnt
            r5_pairs_detail[(tax.nm(dp), tax.nm(x), rel)] += cnt
            if rel == "direct":
                n_r5_pairs_direct += 1
                r5_pair_rows_direct += cnt

    # R6 over distinct (subclass, dp, inodes)
    n_r6 = n_r6_strict = 0
    r6_detail = Counter()
    for (sb, dp, inodes), cnt in isets.items():
        seen.update(inodes)
        if dp is None or dp not in tax.parent:
            continue
        bad = tax.check_r6(sb, dp, inodes) or []
        if bad:
            n_r6 += cnt
            for x in bad:
                r6_detail[(tax.nm(x), tax.depth[x], tax.nm(dp), tax.depth[dp])] += cnt
        if any(x in tax.parent and x not in tax.anc[dp] for x in inodes):
            n_r6_strict += cnt

    prov = {"rows": n, "with_cid": n_cid, "with_zinc": n_zinc, "with_both": n_both,
            "unresolved_slot_rows": n_null_slot,
            "classes_seen": len(seen & set(tax.parent)),
            "classes_total": len(tax.parent) - 1}
    out = {
        "input": a.input, "limit": a.limit, "elapsed_s": round(el, 1),
        "rows_per_s": round(n / el) if el else None,
        "sufficient_statistics": {
            "distinct_label_paths": len(paths),
            "distinct_intermediate_node_sets": len(isets),
            "distinct_dp_alt_pairs": len(dp_alt),
            "total_alternative_parent_labels": n_alt_total,
            "distinct_alt_cooccurrence_pairs": len(alt_pairs) if a.pairs else None,
        },
        "provenance_1b": {"observed": prov, "zenodo_stated": ZENODO,
                          "delta": {k: prov[k] - ZENODO[k] for k in ZENODO if k in prov}},
        "null_by_slot": dict(null_by_slot),
        "direct_parent_depth": {str(k): v for k, v in sorted(dp_depth.items())
                                if k is not None},
        "invariants": {
            "path_violation_rows": n_path_bad,
            "path_violation_rows_fully_populated": n_path_bad_full,
            "path_violation_shapes": len(path_shapes),
            "path_violation_rows_touching_drifted_edges": drift_touched,
            "r5_violating_dp_alt_pairs": n_r5_pairs,
            # summed over violating (dp,alt) pairs: an UPPER BOUND on distinct violating
            # rows, exact unless a row carries two violating pairs at once
            "r5_violating_label_instances": r5_pair_rows,
            "r5_violating_pairs_direct_parent_child": n_r5_pairs_direct,
            "r5_violating_label_instances_direct_parent_child": r5_pair_rows_direct,
            "r6_violation_rows_inclusive": n_r6,
            "r6_violation_rows_strict_tableS1": n_r6_strict,
        },
        "drifted_edges": sorted(tax.nm(x) for x in drift),
        "top_path_shapes": [{"reported": list(g), "expected": list(e), "rows": c}
                            for (g, e), c in path_shapes.most_common(25)],
        "top_r5_pairs": [{"direct_parent": x, "alternative_parent": y, "relation": r,
                          "rows": c} for (x, y, r), c in r5_pairs_detail.most_common(25)],
        "top_r6": [{"intermediate_node": k[0], "depth": k[1], "direct_parent": k[2],
                    "dp_depth": k[3], "rows": c} for k, c in r6_detail.most_common(25)],
    }
    if a.json_only:
        print(json.dumps(out))
        return

    RESULTS.mkdir(exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=2))
    # dump the sufficient statistics so no future change needs another pass
    with (RESULTS / "path_counts.tsv").open("w") as fh:
        fh.write("kingdom\tsuperclass\tclass\tsubclass\tdirect_parent\trows\n")
        for t, c in paths.most_common():
            fh.write("\t".join("" if v is None else str(v) for v in t) + f"\t{c}\n")
    with (RESULTS / "dp_alt_pairs.tsv").open("w") as fh:
        fh.write("direct_parent\talternative_parent\trows\n")
        for (d, x), c in dp_alt.most_common():
            fh.write(f"{'' if d is None else d}\t{x}\t{c}\n")
    with (RESULTS / "inode_sets.tsv").open("w") as fh:
        fh.write("subclass\tdirect_parent\tintermediate_nodes\trows\n")
        for (sb, dp, ino), c in isets.most_common():
            fh.write(f"{'' if sb is None else sb}\t{'' if dp is None else dp}\t"
                     f"{','.join(map(str, ino))}\t{c}\n")
    if a.pairs:
        with (RESULTS / "alt_cooccurrence_pairs.tsv").open("w") as fh:
            fh.write("a\tb\trows\n")
            for (x, y), c in alt_pairs.most_common():
                fh.write(f"{x}\t{y}\t{c}\n")

    print(f"\n{n:,} rows in {el:.0f}s  ({n/el/1e6:.2f} M rows/s)\n")
    print("1b PROVENANCE (Zenodo release notes, not paper claims)")
    for k in ZENODO:
        if k in prov:
            print(f"   {k:24} stated {ZENODO[k]:>12,}   observed {prov[k]:>12,}   "
                  f"delta {out['provenance_1b']['delta'][k]:+,}")
    print("\n1c INVARIANTS (paper claims)")
    for k, v in out["invariants"].items():
        print(f"   {k:48} {v:>12,}")
    print("\nsufficient statistics written (no rescan needed for re-analysis):")
    for f in ("path_counts.tsv", "dp_alt_pairs.tsv", "inode_sets.tsv"):
        print(f"   results/{f}")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
