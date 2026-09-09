"""
Shared ChemOnt taxonomy model for Phase 1b/1c.

Loads the taxonomy from either artefact and exposes the ancestor queries the paper's
invariants are written in terms of. Kept separate from the scan so it can be unit-tested
against hand-built toy taxonomies (see scripts/test_checks.py).
"""

import csv
from pathlib import Path

ROOT_NUM = 9999999


class Taxonomy:
    """A ChemOnt tree keyed by the dataset's compact numeric ids.

    depth(root) = 0, so Kingdom = 1, SuperClass = 2, Class = 3, SubClass = 4, exactly as
    the paper numbers its levels (Methods -> Component 1).
    """

    def __init__(self, parent, name=None):
        self.parent = dict(parent)              # node -> parent node, or -1 for the root
        self.name = dict(name or {})
        self.depth, self.anc_at, self.anc = {}, {}, {}
        for n in self.parent:
            self._depth(n)
        for n in self.parent:
            chain, x = {}, n
            while x != -1:
                chain[self.depth[x]] = x
                x = self.parent[x]
            self.anc_at[n] = chain                       # level -> ancestor-or-self
            self.anc[n] = set(chain.values()) - {n}      # strict ancestors

    def _depth(self, n, guard=None):
        if n in self.depth:
            return self.depth[n]
        guard = guard or set()
        if n in guard:
            raise ValueError(f"cycle through node {n}")
        p = self.parent[n]
        self.depth[n] = 0 if p == -1 else self._depth(p, guard | {n}) + 1
        return self.depth[n]

    def nm(self, n):
        return "NULL" if n is None else self.name.get(n, str(n))

    def expected_path(self, dp):
        """The 5-slot label a compound with this direct parent must carry.

        Slots 1-4 are the direct parent's ancestors at depths 1-4. When the direct parent
        is shallower than a slot, the release pads by repeating it -- a schema fact
        established empirically from the file, not an assumption from the paper.
        """
        d, a = self.depth[dp], self.anc_at[dp]
        return [a[L] if L <= d else dp for L in (1, 2, 3, 4)] + [dp]

    # ---- the paper's stated invariants ------------------------------------------
    def check_path(self, tree):
        """Is the reported 5-level label a genuine root-path of the direct parent?"""
        dp = tree[4]
        if dp is None or dp not in self.parent:
            return None
        return tree == self.expected_path(dp)

    def check_r5(self, dp, alts):
        """Paper: 'Alternative parents are categories that describe the compound but do
        not have an ancestor-descendant relationship with each other or with the Direct
        Parent.'  Returns a set of violated sub-rules."""
        if dp is None or dp not in self.parent:
            return None
        S = {a for a in alts if a in self.parent}
        anc_dp = self.anc[dp]
        bad = set()
        for a in S:
            anc_a = self.anc[a]
            if anc_a & S:
                bad.add("alt_vs_alt")           # an alternative parent is an ancestor of another
            if dp in anc_a:
                bad.add("dp_ancestor_of_alt")   # the direct parent is an ancestor of an alternative
        if anc_dp & S:
            bad.add("alt_ancestor_of_dp")       # an alternative parent is an ancestor of the direct
        return bad

    def check_r6(self, subclass, dp, inodes):
        """Paper (Results) and Additional file 1 Table S1: intermediate nodes are
        'descendants of the subclass and ascendants of the direct parent'.

        Ancestry, not depth. An earlier version tested `depth[n] <= depth[dp]`, which
        accepted 7,970 rows carrying a node at the direct parent's own depth that was
        not on its lineage at all, undercounting the violation total by 703 rows.
        Descendant-or-self of the subclass is accepted because the release lists the
        depth-4 subclass itself among the intermediate nodes.
        """
        if dp is None or dp not in self.parent:
            return None
        bad = []
        for n in inodes:
            if n not in self.parent:
                continue
            ok = (n in self.anc[dp]
                  and subclass is not None
                  and (n == subclass or subclass in self.anc[n]))
            if not ok:
                bad.append(n)
        return bad


def from_dictionary_tsv(path):
    parent, name = {}, {}
    with Path(path).open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n = int(r["numeric_id"])
            name[n] = r["name"]
            p = r["parent_numeric_id"]
            parent[n] = -1 if p in ("null", "", None) else int(p)
    return Taxonomy(parent, name)


def from_obo(path, dict_path):
    """The paper's own ChemOnt 2.1 OBO, re-keyed onto the release's numeric ids so the
    two taxonomies can be compared edge by edge."""
    num, name = {}, {}
    with Path(dict_path).open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            num[r["chemont_id"]] = int(r["numeric_id"])
            name[int(r["numeric_id"])] = r["name"]
    obo, cur = {}, None
    for line in Path(path).open(encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith("["):
            if cur and "id" in cur:
                obo[cur["id"]] = cur
            cur = {} if line == "[Term]" else None
            continue
        if not line.strip() or cur is None:
            continue
        k, _, v = line.partition(": ")
        if k == "id":
            cur["id"] = v
        elif k == "is_a":
            cur["parent"] = v.split("!")[0].strip()
    if cur and "id" in cur:
        obo[cur["id"]] = cur
    parent = {}
    for cid, t in obo.items():
        if cid not in num:
            continue
        p = t.get("parent")
        parent[num[cid]] = num[p] if p in num else -1
    return Taxonomy(parent, name)


def edge_diff(a, b):
    """Nodes whose parent differs between two taxonomies over their shared node set."""
    return {n: (a.parent[n], b.parent[n])
            for n in set(a.parent) & set(b.parent) if a.parent[n] != b.parent[n]}
