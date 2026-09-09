#!/usr/bin/env python3
"""
Test suite for the Phase 1b/1c invariant checks.

Run:  python3 scripts/test_checks.py          (no dependencies, no dataset needed for L1)

Three levels, cheapest first:

  L1  UNIT      hand-built toy taxonomies with known answers. Catches logic errors in
                check_path / check_r5 / check_r6 without touching any real data. These
                are the tests that would have caught the two semantic bugs found while
                developing the scan (see docs/phase1a_audit.md and the notes below).
  L2  FIXTURE   the real ChemOnt taxonomy, hand-verified facts from the paper.
  L3  SAMPLE    the 200k-row sample: schema invariants, plus a determinism check that
                the streaming scan gives identical answers to a naive re-implementation.

Regression notes -- every one of these encodes a bug that was actually made:
  * test_path_padding      : shallow classifications repeat the deepest node in trailing
                             slots. Treating the slots as literal depths reported 7724
                             false violations per 50k rows.
  * test_r6_inclusive      : the release lists the depth-4 subclass itself among the
                             intermediate nodes. A strict depth>4 reading reported 2743
                             false violations per 100k rows.
  * test_r5_directions     : R5 is violated in three distinct directions; an early version
                             checked only alt-vs-alt and missed the dominant case.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chemont  # noqa: E402

# the dictionary ships in-repo (md5-verified against Zenodo); the sibling directory
# is not part of the repo and is absent on any other machine
DICT = HERE.parent / "data" / "chemont_dictionary.tsv"
OBO = HERE.parent / "data" / "ChemOnt_2_1.obo"
SAMPLE = HERE.parent.parent / "chemont_project" / "data" / "sample_200k.tsv"

PASS, FAIL, SKIPPED = [], [], []


def ok(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  -- ' + detail) if detail and not cond else ''}")


# ---------------------------------------------------------------- L1: unit tests
def toy():
    """        root(-1)
                 |
              1 kingdom
                 |
              2 superclass
                /       \
           3 classA    4 classB
              |
           5 subclass
              |
           6 depth5
              |
           7 depth6
    """
    return chemont.Taxonomy(
        {0: -1, 1: 0, 2: 1, 3: 2, 4: 2, 5: 3, 6: 5, 7: 6},
        {0: "root", 1: "kingdom", 2: "super", 3: "classA", 4: "classB",
         5: "subclass", 6: "d5", 7: "d6"})


def l1():
    print("\nL1 UNIT -- toy taxonomies, known answers")
    t = toy()
    ok("depth: root=0, kingdom=1, subclass=4",
       (t.depth[0], t.depth[1], t.depth[5]) == (0, 1, 4))
    ok("ancestors of d6 are exactly {d5, subclass, classA, super, kingdom, root}",
       t.anc[7] == {6, 5, 3, 2, 1, 0})

    # --- path
    ok("path: full-depth label validates", t.check_path([1, 2, 3, 5, 6]) is True)
    ok("test_path_padding: shallow dp repeats in trailing slots",
       t.check_path([1, 2, 4, 4, 4]) is True, "classB is depth 3; slots 4,5 must repeat it")
    ok("path: nulls in trailing slots instead of padding are rejected",
       t.check_path([1, 2, 4, None, 4]) is False)
    ok("path: a null direct parent is unclassifiable, not a violation",
       t.check_path([1, 2, 4, None, None]) is None)
    ok("path: wrong class is caught", t.check_path([1, 2, 4, 5, 6]) is False)
    ok("path: unknown direct parent returns None", t.check_path([1, 2, 3, 5, 999]) is None)

    # --- R5, all three directions
    ok("test_r5_directions/clean", t.check_r5(6, [4]) == set())
    ok("test_r5_directions/alt_vs_alt", t.check_r5(4, [3, 5]) == {"alt_vs_alt"},
       "classA is an ancestor of subclass, and neither relates to classB")
    ok("R5: sub-rules are reported independently, not collapsed",
       t.check_r5(4, [2, 3]) == {"alt_vs_alt", "alt_ancestor_of_dp"},
       "super is an ancestor of both classA and the direct parent classB")
    ok("test_r5_directions/dp_ancestor_of_alt", t.check_r5(3, [6]) == {"dp_ancestor_of_alt"},
       "classA is an ancestor of d5 -- the dominant real-world case")
    ok("test_r5_directions/alt_ancestor_of_dp", t.check_r5(6, [3]) == {"alt_ancestor_of_dp"})
    ok("R5: unknown ids are skipped, not crashed on", t.check_r5(6, [4, 12345]) == set())

    # --- R6
    ok("test_r6_inclusive: the subclass itself is an accepted intermediate node",
       t.check_r6(5, 7, [5]) == [], "release lists the depth-4 node; strict depth>4 is wrong")
    ok("R6: a genuine depth-5 intermediate node passes", t.check_r6(5, 7, [6]) == [])
    ok("R6: a node deeper than the direct parent is rejected", t.check_r6(5, 6, [7]) == [7])
    ok("R6: a node outside the subclass subtree is rejected", t.check_r6(5, 7, [4]) == [4])

    try:
        chemont.Taxonomy({0: 1, 1: 0})
        ok("cycles are detected", False, "no exception raised")
    except ValueError:
        ok("cycles are detected", True)


# ------------------------------------------------------------- L2: real taxonomy
def l2():
    print("\nL2 FIXTURE -- real ChemOnt, facts taken from the paper")
    if not DICT.exists():
        print(f"  SKIP  L2 needs {DICT}")
        SKIPPED.append("L2"); return None, None
    d = chemont.from_dictionary_tsv(DICT)
    o = chemont.from_obo(OBO, DICT) if OBO.exists() else None
    ok("dictionary holds 4825 terms (root + 4824 categories)", len(d.parent) == 4825,
       str(len(d.parent)))
    ok("exactly one root", sum(1 for v in d.parent.values() if v == -1) == 1)
    ok("two kingdoms at depth 1", sum(1 for n in d.parent if d.depth[n] == 1) == 2)
    ok("31 superclasses at depth 2", sum(1 for n in d.parent if d.depth[n] == 2) == 31)
    ok("1729 subclasses at depth 4", sum(1 for n in d.parent if d.depth[n] == 4) == 1729)
    ok("maximum depth 11", max(d.depth.values()) == 11)
    if o:
        diff = chemont.edge_diff(d, o)
        ok("release and OBO 2.1 differ on exactly 4 parent edges", len(diff) == 4,
           f"{len(diff)}: {[d.nm(n) for n in diff]}")
        print(f"        drifted: {sorted(d.nm(n) for n in diff)}")
    return d, o


# ------------------------------------------------------------------ L3: sample
def l3(d):
    print("\nL3 SAMPLE -- 200k rows, schema + cross-implementation agreement")
    if d is None or not SAMPLE.exists():
        # a skipped tier used to leave the suite exiting 0, so CI reported success
        # while only L1 had actually run
        print(f"  SKIP  L3 needs {SAMPLE}; run `make data` or fetch the sample")
        SKIPPED.append("L3"); return
    import csv as _csv
    n = bad_len = unknown = null_kingdom = 0
    naive = {"path": 0, "r5": 0, "r5pairs": 0, "r6": 0}
    with SAMPLE.open() as fh:
        for i, r in enumerate(_csv.DictReader(fh, delimiter="\t")):
            if i >= 20000:
                break
            n += 1
            t = json.loads(r["chemont_tree_json"])
            if len(t) != 5:
                bad_len += 1
            if t[0] is None:
                null_kingdom += 1
            ids = [x for x in t if x is not None]
            if any(x not in d.parent for x in ids):
                unknown += 1
            o = json.loads(r["chemont_other_json"])
            if d.check_path(t) is False:
                naive["path"] += 1
            alts = o.get("alternative_parents") or []
            if d.check_r5(t[4], alts):
                naive["r5"] += 1
            dp = t[4]
            if dp is not None and dp in d.parent:
                for a in alts:
                    if a in d.parent and a != dp and (a in d.anc[dp] or dp in d.anc[a]):
                        naive["r5pairs"] += 1
            if d.check_r6(t[3], t[4], o.get("intermediate_nodes") or []):
                naive["r6"] += 1
    ok("every label array has exactly 5 slots", bad_len == 0, str(bad_len))
    ok("every referenced id exists in the dictionary", unknown == 0, str(unknown))
    ok("kingdom slot is never null", null_kingdom == 0, str(null_kingdom))
    print(f"        naive re-implementation over {n} rows: {naive}")

    scan = HERE / "02_scan_dataset.py"
    if scan.exists():
        out = subprocess.run(
            [sys.executable, str(scan), "--input", str(SAMPLE), "--limit", "20000",
             "--quiet", "--json-only"],
            capture_output=True, text=True)
        if out.returncode != 0:
            ok("streaming scan runs on the sample", False, out.stderr.strip()[-300:])
            return
        got = json.loads(out.stdout)
        inv = got["invariants"]
        # r5_violating_label_instances counts (dp, alt) PAIR instances while naive["r5"]
        # counts ROWS, so a >= comparison passed for any inflated value. Compare against
        # a naive pair count, which actually bites.
        agree = (inv["path_violation_rows"] == naive["path"]
                 and inv["r5_violating_label_instances"] == naive["r5pairs"]
                 and inv["r6_violation_rows_inclusive"] == naive["r6"])
        ok("streaming scan agrees with the naive re-implementation "
           "(R5 as a documented upper bound)", agree,
           f"scan={inv} naive={naive}")
        rows = got["provenance_1b"]["observed"]["rows"]
        ok("streaming scan counted the right number of rows", rows == n, f"{rows} vs {n}")


if __name__ == "__main__":
    l1()
    d, _ = l2()
    l3(d)
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed"
          + (f", {len(SKIPPED)} tier(s) skipped: {', '.join(SKIPPED)}" if SKIPPED else ""))
    if FAIL:
        print("failed: " + ", ".join(FAIL))
    if SKIPPED and os.environ.get("ALLOW_SKIPPED_TIERS") != "1":
        print("a skipped tier is not a pass; set ALLOW_SKIPPED_TIERS=1 to accept")
        sys.exit(2)
    sys.exit(1 if FAIL else 0)
