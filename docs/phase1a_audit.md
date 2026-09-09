# Phase 1a — self-audit of the replication logic

Findings from auditing `scripts/01_verify_taxonomy.py` against the OBO. Each flaw is
recorded with its effect on the reported verdicts. Fixed in `01_verify_taxonomy.py` v2.

## F1 — the same quantity was counted twice, once as a pass and once as a fail
The inorganic category count appears as **C2-inorg** (paper 678 → MATCH) and again as
**C11** (paper 674 → MISMATCH). These are one measurement, not two. The paper contradicts
*itself* (Methods 678, Conclusion 674); that is a property of the paper, not two
independent tests of ours. Counting both inflated the numerator *and* the denominator.
**Effect:** the honest denominator is 25, not 26.

## F2 — a derived ratio was scored as an independent pass
**D3-avg** (1.24 ChEBI synonyms per category) was scored MATCH at 1.2511, while its own
numerator **D3-count** (6014) had already FAILED at 6034. A ratio cannot reproduce when
the count it is built from does not. It is a *dependent* check and must not be counted as
independent corroboration.
**Effect:** D3-avg reclassified as DERIVED, excluded from the score.

## F3 — two pass criteria were tolerances I invented, not the paper's
- **C9** "average depth of five levels per node": passed via `abs(mean − 5) < 0.5`. The
  observed mean is 4.636. The 0.5 window was chosen by me and has no basis in the paper.
- **D3-avg**: passed via `< 0.02`, likewise invented.
The correct treatment is to report the number and say plainly whether it is consistent
with the paper's stated rounding, not to manufacture a threshold that makes it pass.
**Effect:** C9 reclassified from MATCH to CONSISTENT (4.636 rounds to 5; it is not 5).

## F4 — "mappings" was measured as synonym *lines*, but a line can carry several IDs
A `synonym:` line ends in a bracketed xref list which may hold **more than one** external
identifier. Measured three ways, the MeSH figures are:

| reading | value |
|---|---|
| synonym lines tagged MeSH_TERM | 966 |
| distinct (category, MeSH-ID) pairs | **1165** |
| distinct MeSH IDs | 514 |

I reported 966 as if it were the only reading. The paper says "844 ClassyFire categories
have been mapped to at least one corresponding MeSH term, accounting for a total of 945
mappings" — 945 > 844, so *lines* is the closest reading, but this had to be argued, not
assumed. ChEBI is barely affected (6034 lines vs 6031 pairs); MeSH is affected by 20 %.
**Effect:** all three readings are now reported for every mapped source.

## F5 — sibling checks used inconsistent units
**D5-terms** (LIPID MAPS "307 terms") was measured as distinct case-folded synonym *text*,
while D3/D6 used *lines*. Different units for the same kind of question. LIPID MAPS gives
832 lines / 832 pairs / **318 distinct IDs** / 309 distinct texts. The paper's 307 is
nearest the *text* reading, which is a real (and reportable) result — but only once all
readings are on the table.

## F6 — an undocumented inclusion choice, and no provenance for the input file
- The root term `CHEMONTID:9999999` carries **one** synonym line of its own
  (`"a chemical entity" RELATED CHEMONT_TERM`). The file holds 9025 synonym lines; 9024
  are on non-root categories. I silently used 9024 without saying so. Neither resolves the
  paper's 9012 (Δ = 12 or 13).
- More seriously: the OBO was copied from an earlier working directory with **no recorded
  source URL and no checksum**. The whole of Phase 1a rests on it. Now pinned:
  `sha256 = 8616a6ecb96c8aeb204739a4de045cd290eab9a0e164d0ceaa5d888d0751fe22`,
  `data-version: 2.1`, `date: 26:08:2016 16:28`, `auto-generated-by: OBO-Edit 2.3.1`.
  Re-download from the ClassyFire site and confirm this digest before the report is
  submitted.

## F7 — trivial claims padded the score
`C1` (">4800"), `D2` (">9000") and `D8` (scope vocabulary ⊆ EXACT/NARROW/BROAD/RELATED)
are satisfied by almost any file with the right shape. Scoring them alongside "the
taxonomy is a strict tree of depth 11" overstates the replication.
**Effect:** claims are now tiered STRUCTURAL / BOOKKEEPING / TRIVIAL and scored separately.

## F8 — latent correctness bugs (no effect on current results, but unsound)
- `get_depth` follows `parents[0]` only. Sound only because zero multi-parent nodes exist;
  that fact was checked but never *asserted* before the traversal used it.
- `alt_id` is parsed and discarded. `CHEMONTID:0000284` ("Aniline and substituted
  anilines") carries `alt_id: CHEMONTID:0003923`, so the ID space is 4825 primary + 1
  retired = 4826 identifiers. This does **not** explain the Class off-by-one: a merge makes
  the released file *smaller*, whereas the paper's 764 is *smaller* than the observed 765.

## What the audit did NOT overturn
The exact structural matches stand: strict tree, 4146/678 organic-inorganic partition
(which sums to 4824 and thereby pins the category count independently), 2 kingdoms,
31 SuperClasses = 26 + 5, 1729 SubClasses, max depth 11, and the Phenylpropanoids
worked example at 34 children / 273 descendants. Those are single-definition counts with
no ambiguity, and they are the claims the paper's argument actually rests on.

## The one gap that cannot be closed from the manuscript text
The Class (+1) and below-SubClass (+1) discrepancies are unresolvable against the prose
alone. The paper points to **Additional file 1: Table S1** ("A more complete description of
this taxonomic hierarchy"), which is the authors' own per-level listing. That supplement is
part of the paper and is the correct instrument to settle whether the +1 is version drift
or a manuscript slip. It should be obtained before the final report claims either.
