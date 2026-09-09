#!/usr/bin/env python3
"""
Phase 1a (v2) — reproduce the ChemOnt structural claims of the baseline paper.

Baseline paper:
  Djoumbou Feunang et al. (2016) ClassyFire: automated chemical classification with a
  comprehensive, computable taxonomy. J Cheminform 8:61. doi:10.1186/s13321-016-0174-y

Input   : data/ChemOnt_2_1.obo
          data-version 2.1, date 26:08:2016 16:28, OBO-Edit 2.3.1
          sha256 8616a6ecb96c8aeb204739a4de045cd290eab9a0e164d0ceaa5d888d0751fe22
Output  : results/taxonomy_claims.json, results/taxonomy_nodes.tsv, verdict table on stdout

v2 incorporates the self-audit in docs/phase1a_audit.md. Changes from v1:
  F1  the inorganic count is scored once, not twice (the paper's 678/674 self-contradiction
      is a property of the paper and is reported as a note, not as a second test).
  F2  ratios derived from a failing count are marked DERIVED and excluded from the score.
  F3  invented numeric tolerances removed; a claim that only holds under rounding is
      reported as CONSISTENT, never as MATCH.
  F4  every cross-ontology count is reported under all three defensible readings
      (synonym lines / distinct (category, id) pairs / distinct ids), because a single
      synonym line may carry several bracketed identifiers.
  F5  sibling checks use consistent units.
  F6  root-term inclusion is stated explicitly; input file provenance is pinned above.
  F7  claims are tiered STRUCTURAL / BOOKKEEPING / TRIVIAL and scored separately.
  F8  the single-parent property is asserted before depth traversal relies on it;
      alt_id is retained and reported.

Conventions, all forced by the paper's own text (Methods -> Component 1):
  * "Kingdom, SuperClass, Class, and SubClass ... denote the first, second, third and
    fourth levels", so depth(Kingdom) = 1 and depth(root) = 0.
  * A *category* is any [Term] other than the root CHEMONTID:9999999, because the paper
    counts categories "in addition to the root category (Chemical entities)".
"""

import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = "CHEMONTID:9999999"
HERE = Path(__file__).resolve().parent.parent
OBO = HERE / "data" / "ChemOnt_2_1.obo"
RESULTS = HERE / "results"
EXPECTED_SHA = "8616a6ecb96c8aeb204739a4de045cd290eab9a0e164d0ceaa5d888d0751fe22"

STRUCTURAL, BOOKKEEPING, TRIVIAL = "STRUCTURAL", "BOOKKEEPING", "TRIVIAL"
MATCH, MISMATCH, CONSISTENT, DERIVED = "MATCH", "MISMATCH", "CONSISTENT", "DERIVED"


def parse_obo(path):
    header, terms, cur, in_header = {}, {}, None, True
    for raw in path.open(encoding="utf-8"):
        line = raw.rstrip("\n")
        if line.startswith("["):
            in_header = False
            if cur is not None and "id" in cur:
                terms[cur["id"]] = cur
            cur = {"syn": [], "xref": [], "alt_id": []} if line == "[Term]" else None
            continue
        if not line.strip():
            continue
        key, _, val = line.partition(": ")
        if in_header:
            header.setdefault(key, []).append(val)
            continue
        if cur is None:
            continue
        if key in ("id", "name", "def"):
            cur[key] = val
        elif key == "is_a":
            cur.setdefault("parents", []).append(val.split("!")[0].strip())
        elif key == "synonym":
            cur["syn"].append(val)
        elif key in ("xref", "alt_id"):
            cur[key].append(val)
    if cur is not None and "id" in cur:
        terms[cur["id"]] = cur
    return header, terms


SYN_RE = re.compile(r'^"(.*?)"\s+(\w+)(?:\s+(\w+))?\s*(?:\[(.*)\])?\s*$')


def parse_synonym(s):
    """-> (text, scope, source|None, [external ids]).  A line may carry several ids."""
    m = SYN_RE.match(s)
    if not m:
        return s, None, None, []
    ids = [x.strip() for x in (m.group(4) or "").split(",") if x.strip()]
    return m.group(1), m.group(2), m.group(3), ids


def main():
    if not OBO.exists():
        sys.exit(f"missing {OBO}")
    RESULTS.mkdir(exist_ok=True)

    # ---- F6: provenance is verified, not assumed -----------------------------
    digest = hashlib.sha256(OBO.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA:
        sys.exit(f"input file digest changed:\n  expected {EXPECTED_SHA}\n  got      {digest}")

    header, terms = parse_obo(OBO)
    categories = [i for i in terms if i != ROOT]

    # ---- F8: assert the tree property before depth relies on it --------------
    multi = {i: t["parents"] for i, t in terms.items() if len(t.get("parents", [])) > 1}
    dangling = sorted({p for t in terms.values() for p in t.get("parents", [])
                       if p not in terms})
    parentless = [i for i in terms if not terms[i].get("parents")]
    is_tree = (not multi) and (not dangling) and parentless == [ROOT]
    assert not multi, f"depth traversal assumes a single parent; found {len(multi)} exceptions"

    children = defaultdict(list)
    for i, t in terms.items():
        for p in t.get("parents", []):
            children[p].append(i)

    depth, seen = {}, set()

    def get_depth(i):
        if i in depth:
            return depth[i]
        assert i not in seen, f"cycle through {i}"
        seen.add(i)
        ps = terms[i].get("parents")
        depth[i] = 0 if not ps else get_depth(ps[0]) + 1
        seen.discard(i)
        return depth[i]

    for i in terms:
        get_depth(i)

    by_level = Counter(depth[i] for i in categories)

    def subtree(r):
        n, st = 0, [r]
        while st:
            x = st.pop()
            n += 1
            st.extend(children[x])
        return n

    org = next(i for i in children[ROOT] if terms[i]["name"] == "Organic compounds")
    inorg = next(i for i in children[ROOT] if terms[i]["name"] == "Inorganic compounds")
    n_org, n_inorg = subtree(org), subtree(inorg)

    def kingdom_of(i):
        while terms[i].get("parents") and terms[i]["parents"][0] != ROOT:
            i = terms[i]["parents"][0]
        return terms[i]["name"]

    sc_split = Counter(kingdom_of(i) for i in categories if depth[i] == 2)

    # ---- F3: report depth statistics under every reasonable reading ----------
    leaves = [i for i in categories if not children[i]]
    depth_stats = {
        "mean_over_categories": sum(depth[i] for i in categories) / len(categories),
        "mean_over_leaves": sum(depth[i] for i in leaves) / len(leaves),
        "mean_incl_root": sum(depth.values()) / len(depth),
        "max": max(by_level),
    }

    # ---- F4/F5: cross-ontology counts under all three readings ---------------
    mapping = {}
    n_syn_all = sum(len(t["syn"]) for t in terms.values())
    n_syn_cat = sum(len(terms[i]["syn"]) for i in categories)
    scopes, sources = Counter(), Counter()
    for i in categories:
        for s in terms[i]["syn"]:
            text, scope, source, ids = parse_synonym(s)
            scopes[scope] += 1
            sources[source or "UNTYPED"] += 1
            d = mapping.setdefault(source or "UNTYPED",
                                   {"lines": 0, "pairs": set(), "ids": set(),
                                    "texts": set(), "cats": set()})
            d["lines"] += 1
            d["cats"].add(i)
            d["texts"].add(text.strip().lower())
            for x in ids:
                d["pairs"].add((i, x))
                d["ids"].add(x)
    readings = {k: {"lines": v["lines"], "pairs": len(v["pairs"]), "ids": len(v["ids"]),
                    "texts": len(v["texts"]), "categories": len(v["cats"])}
                for k, v in mapping.items()}

    pp = next(i for i, t in terms.items()
              if t.get("name") == "Phenylpropanoids and polyketides")
    no_chebi = [i for i in categories if not any("ChEBI_TERM" in s for s in terms[i]["syn"])]
    alt_ids = {i: t["alt_id"] for i, t in terms.items() if t["alt_id"]}

    C = []

    def claim(cid, tier, paper, obs, verdict, note="", alt=None):
        C.append({"id": cid, "tier": tier, "paper": paper, "observed": obs,
                  "verdict": verdict, "note": note, "other_readings": alt})

    # -- STRUCTURAL: single-definition counts, no interpretive freedom ---------
    claim("C10", STRUCTURAL, "strict tree", "strict tree" if is_tree else "NOT a tree",
          MATCH if is_tree else MISMATCH,
          f"{len(multi)} multi-parent nodes, {len(dangling)} dangling parents, "
          f"parentless = {parentless} (expected exactly the root)")
    claim("C3", STRUCTURAL, 2, by_level[1], MATCH if by_level[1] == 2 else MISMATCH,
          "Kingdom level: " + ", ".join(sorted(terms[i]["name"] for i in children[ROOT])))
    claim("C4", STRUCTURAL, "31 = 26 organic + 5 inorganic", by_level[2],
          MATCH if (by_level[2] == 31 and sc_split["Organic compounds"] == 26
                    and sc_split["Inorganic compounds"] == 5) else MISMATCH,
          f"observed split {dict(sc_split)}")
    claim("C6", STRUCTURAL, 1729, by_level[4], MATCH if by_level[4] == 1729 else MISMATCH,
          "SubClass level (depth 4)")
    claim("C5", STRUCTURAL, 764, by_level[3], MATCH if by_level[3] == 764 else MISMATCH,
          "Class level (depth 3). Unresolved +1 — see docs/phase1a_audit.md; only the "
          "paper's Additional file 1 Table S1 can settle drift vs manuscript slip.")
    claim("C7", STRUCTURAL, 2296, sum(v for k, v in by_level.items() if k >= 5),
          MATCH if sum(v for k, v in by_level.items() if k >= 5) == 2296 else MISMATCH,
          "levels 5-11. Same unresolved +1 as C5.")
    claim("C8", STRUCTURAL, 11, depth_stats["max"],
          MATCH if depth_stats["max"] == 11 else MISMATCH, "maximum depth")
    claim("C2-split", STRUCTURAL, "4146 organic / 678 inorganic", f"{n_org} / {n_inorg}",
          MATCH if (n_org, n_inorg) == (4146, 678) else MISMATCH,
          "F1: the paper's Conclusion separately says 674 inorganic, contradicting its own "
          "Methods figure of 678. That is an internal inconsistency in the paper, not a "
          "second measurement; the OBO gives 678. Scored once.")
    claim("C2-total", STRUCTURAL, 4825, len(categories),
          MATCH if len(categories) == 4825 else MISMATCH,
          f"{len(terms)} [Term] stanzas = 1 root + {len(categories)} categories. "
          f"4146 + 678 = {n_org + n_inorg} independently confirms {len(categories)}, so the "
          "paper's 4825 counts the root as a category despite saying 'in addition to' it.")
    claim("C12", STRUCTURAL, "34 children / 273 descendants",
          f"{len(children[pp])} / {subtree(pp) - 1}",
          MATCH if (len(children[pp]), subtree(pp) - 1) == (34, 273) else MISMATCH,
          "Phenylpropanoids and polyketides, the paper's worked example")
    claim("D7", STRUCTURAL, len(categories),
          sum(1 for i in categories if terms[i].get("def")),
          MATCH if sum(1 for i in categories if terms[i].get("def")) == len(categories)
          else MISMATCH, "every category carries a text definition")

    # -- BOOKKEEPING: counts whose definition the paper leaves ambiguous -------
    claim("C9", BOOKKEEPING, "~5 (average depth per node)",
          round(depth_stats["mean_over_categories"], 3), CONSISTENT,
          "F3: no tolerance is invented. 4.636 rounds to 5 but is not 5; reported as "
          "consistent with the paper's stated rounding, not as a reproduction.",
          {k: round(v, 3) for k, v in depth_stats.items()})
    claim("D1", BOOKKEEPING, 9012, n_syn_cat, MISMATCH,
          f"F6: {n_syn_all} synonym lines file-wide; {n_syn_cat} on non-root categories "
          f"(the root carries 1 of its own). Neither equals 9012. The paper says "
          f"'English synonyms', which may exclude the {sources['UNTYPED']} untyped lines "
          "or some tagged subset; no subset we can define lands on 9012.")
    claim("D3", BOOKKEEPING, "6014 ChEBI mappings", readings["ChEBI_TERM"]["lines"], MISMATCH,
          "F4: three readings reported; all exceed 6014.",
          readings["ChEBI_TERM"])
    claim("D3-avg", BOOKKEEPING, 1.24,
          round(readings["ChEBI_TERM"]["lines"] / readings["ChEBI_TERM"]["categories"], 4),
          DERIVED,
          "F2: a ratio built on D3, which failed. Not independent evidence and not scored. "
          "Note the paper's own pair is near-inconsistent: 6014 / 4825 = 1.2465, which the "
          "paper rounds down to 1.24.")
    claim("D4", BOOKKEEPING, "every category has >=1 ChEBI term",
          f"{readings['ChEBI_TERM']['categories']} of {len(categories)}", MISMATCH,
          "exception: " + ", ".join(f"{i} ({terms[i]['name']})" for i in no_chebi))
    claim("D5", BOOKKEEPING, "789 categories / 307 LIPID MAPS terms",
          f"{readings['LIPIDMAPS_TERM']['categories']} categories / "
          f"{readings['LIPIDMAPS_TERM']['texts']} distinct terms", MISMATCH,
          "F5: 'terms' read as distinct case-folded synonym text, the reading nearest the "
          "paper; by distinct bracketed ID it is "
          f"{readings['LIPIDMAPS_TERM']['ids']}.", readings["LIPIDMAPS_TERM"])
    claim("D6", BOOKKEEPING, "844 categories / 945 MeSH mappings",
          f"{readings['MeSH_TERM']['categories']} categories / "
          f"{readings['MeSH_TERM']['lines']} mappings", MISMATCH,
          "F4: 'mappings' read as synonym lines because the paper's 945 exceeds its 844 "
          "categories; by distinct (category, MeSH-ID) pair it is "
          f"{readings['MeSH_TERM']['pairs']}, a 20% difference. The reading matters here.",
          readings["MeSH_TERM"])

    # -- TRIVIAL: satisfied by almost any file of the right shape --------------
    claim("C1", TRIVIAL, ">4800", len(categories),
          MATCH if len(categories) > 4800 else MISMATCH, "F7: near-vacuous")
    claim("D2", TRIVIAL, ">9000", n_syn_cat, MATCH if n_syn_cat > 9000 else MISMATCH,
          "F7: near-vacuous")
    ok_scopes = set(k for k in scopes if k) <= {"EXACT", "NARROW", "BROAD", "RELATED"}
    claim("D8", TRIVIAL, "scopes subset of EXACT/NARROW/BROAD/RELATED",
          sorted(k for k in scopes if k), MATCH if ok_scopes else MISMATCH,
          "F7: vocabulary check")

    # ---- report --------------------------------------------------------------
    print(f"\nChemOnt {header.get('data-version', ['?'])[0]}  ({header.get('date', ['?'])[0]})")
    print(f"  file    {OBO.name}")
    print(f"  sha256  {digest}  [verified]")
    print(f"  {len(terms)} [Term] stanzas = 1 root + {len(categories)} categories"
          f"; {len(alt_ids)} retired alt_id -> {sum(len(v) for v in alt_ids.values())} "
          f"extra identifier(s)\n")

    for tier in (STRUCTURAL, BOOKKEEPING, TRIVIAL):
        rows = [r for r in C if r["tier"] == tier]
        scored = [r for r in rows if r["verdict"] in (MATCH, MISMATCH)]
        n_ok = sum(r["verdict"] == MATCH for r in scored)
        print(f"--- {tier}  ({n_ok}/{len(scored)} reproduced exactly) " + "-" * 24)
        for r in rows:
            print(f"  {r['id']:<10} paper: {str(r['paper']):<38} "
                  f"observed: {str(r['observed']):<30} {r['verdict']}")
        print()

    print("Notes on every non-MATCH:")
    for r in C:
        if r["verdict"] != MATCH:
            print(f"  [{r['id']}] {r['verdict']}: paper {r['paper']} vs observed {r['observed']}")
            if r["note"]:
                for chunk in [r["note"][i:i + 92] for i in range(0, len(r["note"]), 92)]:
                    print(f"      {chunk}")
            if r["other_readings"]:
                print(f"      other readings: {r['other_readings']}")

    scored = [r for r in C if r["verdict"] in (MATCH, MISMATCH)]
    struct = [r for r in scored if r["tier"] == STRUCTURAL]
    print(f"\nHeadline: {sum(r['verdict'] == MATCH for r in struct)}/{len(struct)} "
          f"STRUCTURAL claims reproduced exactly; "
          f"{sum(r['verdict'] == MATCH for r in scored)}/{len(scored)} overall "
          f"(1 derived and 1 rounding-only claim excluded from both).")

    out = {
        "provenance": {"file": OBO.name, "sha256": digest,
                       "data_version": header.get("data-version", [None])[0],
                       "date": header.get("date", [None])[0],
                       "generated_by": header.get("auto-generated-by", [None])[0]},
        "n_term_stanzas": len(terms), "n_categories": len(categories),
        "alt_ids": alt_ids,
        "tree": {"is_tree": is_tree, "multi_parent": len(multi),
                 "dangling": dangling, "parentless": parentless},
        "level_histogram": {str(k): v for k, v in sorted(by_level.items())},
        "level_histogram_by_kingdom": {
            str(L): dict(Counter(kingdom_of(i) for i in categories if depth[i] == L))
            for L in sorted(by_level)},
        "depth_stats": depth_stats,
        "n_organic": n_org, "n_inorganic": n_inorg,
        "synonyms": {"lines_file_wide": n_syn_all, "lines_on_categories": n_syn_cat,
                     "by_scope": dict(scopes), "by_source": dict(sources),
                     "readings_per_source": readings},
        "categories_without_chebi": {i: terms[i]["name"] for i in no_chebi},
        "claims": C,
    }
    (RESULTS / "taxonomy_claims.json").write_text(json.dumps(out, indent=2))

    with (RESULTS / "taxonomy_nodes.tsv").open("w") as fh:
        fh.write("chemontid\tname\tparent\tdepth\tkingdom\tn_children\tn_descendants\n")
        for i in sorted(terms):
            p = terms[i]["parents"][0] if terms[i].get("parents") else ""
            fh.write(f"{i}\t{terms[i].get('name','')}\t{p}\t{depth[i]}\t"
                     f"{kingdom_of(i) if i != ROOT else ''}\t{len(children[i])}\t"
                     f"{subtree(i) - 1}\n")
    print(f"\nwrote {RESULTS/'taxonomy_claims.json'}, {RESULTS/'taxonomy_nodes.tsv'}")


if __name__ == "__main__":
    main()
