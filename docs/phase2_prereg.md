# Phase 2 — pre-registration

Written before any validation code (`scripts/12_validate_flagged.py`) existed or was run.
Implements `docs/phase2_plan.md` §3–§7. **Frozen by hash, not by commit**: the user has not
authorised git commits, so the sha256 of this file and an ISO timestamp are computed
immediately before the primary test and written into `results/phase2/validation.json`
(`prereg.sha256`, `prereg.frozen_at`). Any later edit changes the hash and is detectable.
One section (§9) is appended after calibration and before validation; the hash of the
file *without* §9 is recorded too (`results/phase2/prereg_hash_initial.txt`).

## 1. Objects

* Tree: `data/chemont_dictionary.tsv` via `chemont.from_dictionary_tsv` (release ids,
  root 9999999 at depth 0, Kingdom = 1 … SubClass = 4). Never modified (C4).
* S1: `results/dp_alt_pairs.tsv`. The 5 rows with an empty direct parent (35 labels) are
  dropped. `n(d,a)` = rows with direct parent `d` listing `a` as alternative parent.
* Undirected pair `{u,v}`: `w_raw = n(u,v) + n(v,u)`.
* `PMI(d,a) = log2[ (n(d,a)/Σ_x n(d,x)) / (Σ_x n(x,a)/Σ_{x,y} n(x,y)) ]`.
  `PMI_sym{u,v} = min` over the directions that exist (one-directional pairs use their
  single direction — required, else every parent→child R5 pair, whose reverse is never
  observed, would be excluded by construction).
* `w_pmi = max(0, PMI_sym)`, `w_cond(d→a) = n(d,a)/Σ_x n(d,x)`.
* LCA depth from the tree; `sibling` = same parent.
* S2 (alt–alt co-occurrence, `results/alt_cooccurrence_pairs.tsv`, if available):
  `PMI(a,b) = log2[ C(a,b)·N / (f(a) f(b)) ]`, `f(a) = Σ_x n(x,a)`, `N = 73,105,281`.
  Filtered by the same A∩B rule. **Sensitivity channel only**; never enters H1a.
* S3: R5 ancestor pairs, R6 (intermediate node, direct parent) violations, and
  non-null mismatched slots of PATH-violating shapes. Separate tier, never merged into
  S1 weights.

## 2. Separation (§3)

* Criterion A(θ): `LCA_depth{u,v} ≥ θ`.
* Criterion B(τ, n_min): `PMI_sym{u,v} ≥ τ` **and** `w_raw ≥ n_min`.
* **Primary overlay `G_c` = S1 edges in A(θ=4) ∩ B(τ*, n_min*)**, weight `w_pmi`.
  A-only and B-only are sensitivity variants.
* τ*, n_min* are chosen **only** from the calibration set (§3), never from the 21/20
  expert classes, by the rule below.

## 3. Calibration rule (fixed before calibration is computed)

Calibration set (plan §3.4): (i) the R5-violating (dp, alt) pairs, recomputed from
`dp_alt_pairs.tsv` + tree (expected 93); (ii) every pair (reported node, expected node)
at a non-null mismatched slot of a PATH-violating label shape in `results/path_counts.tsv`
(expected 32 shapes). Collapsed to undirected pairs; only pairs present in S1 can be
evaluated.

For each of the 9 cells τ ∈ {1,2,3} × n_min ∈ {10,100,1000}, criterion **B alone**:
`recall = |cal ∩ B| / |cal ∩ S1|`, `retain = |B| / |S1|`, `lift = recall / retain`.
**Choose the cell with maximum lift among cells with recall ≥ 0.5** (precision proxy
under a positive-unlabelled design, with a recall floor). Ties → higher recall, then
smaller τ. If no cell reaches 0.5, choose maximum recall.

## 4. Primary confirmatory test — H1a (run once)

* Flagged set: every distinct FP/FN name in Additional file 4 sheet 2 (cols D/E),
  resolved name → id by exact name, then case-insensitive name, then OBO synonym, then
  documented manual adjudication. Pre-checked join (names only, no scores): 26 distinct
  names, 25 exact + 1 case-insensitive (`Transition Metal Hydrides` → 536), 26/26
  resolve; **21** act as a direct parent. **Primary n = 21.** The plan's exact-join
  n = 20 is reported as a secondary row.
* Population: the 3,604 classes that act as a direct parent in S1.
* Score: `H_f(d) = −Σ p log2 p` over `p(a|d)` renormalised to alternatives `a` with
  `{d,a} ∈ G_c`; `H_f = 0` when fewer than two retained alternatives.
* Statistic: `T = mean H_f` over the flagged classes.
* Null: size-stratified permutation. Size = compounds with `d` as direct parent
  (`path_counts.tsv`). Deciles by rank over the 3,604 classes. Each control set draws,
  for each flagged class, one non-flagged class from the same decile (without
  replacement within a set). R = 10,000, seed 798.
* `p = (1 + #{T_null ≥ T_obs}) / (1 + R)`, one-sided. **α = 0.05.**
  Supported iff p < 0.05.

Secondary, labelled as such: uncontrolled p (controls from all non-flagged classes);
the same two tests on the unfiltered entropy (2×2 table, §6.3); n = 20 exact-join set;
alternative score `1 − max_a p(a|d)` on `G_c`; bootstrap 95 % CI of
`T_obs − E[T_null]` (10,000 resamples of the flagged classes); percentiles of the
positive controls Sesquiterpenoids / Diterpenoids; leakage delta with the test
compounds removed; ChEBI axis (Additional file 5, 20 *compounds*, Spearman between
disagreement `(missed + suggestions)/(ChEBI + suggestions)` and `H_f` of the compound's
release direct parent; exploratory).

## 5. H1b (structure)

On `G_c` minus isolated nodes: Onnela weighted mean clustering (`w_pmi`) and Louvain
modularity (networkx, weighted, resolution 1, seed 0). Null: degree-preserving
rewiring (`double_edge_swap`, 10·|E| swaps) with weights permuted over rewired edges,
200 replicates (fewer only if runtime forces it; stated). Supported iff z > 3 for both;
not supported iff both observed values lie within the null's central 95 %; else
inconclusive. Binary transitivity, ER and Chung–Lu (label-frequency) nulls reported
alongside, not decisive.

## 6. H1c (alignment)

NMI (arithmetic normalisation) and ARI between Louvain communities of `G_c` and the
SuperClass of each node (Kingdom-level nodes labelled by their kingdom). Supported iff
NMI < 0.7; not supported iff NMI ≥ 0.9; inconclusive in between.

## 7. Sensitivity sweep (run only after the primary test)

θ ∈ {2,3,4,5} × τ ∈ {1,2,3} × n_min ∈ {10,100,1000} = 36 cells of A(θ) ∩ B(τ,n_min).
Per cell: edges, Jaccard(A,B), |A∩B|/|A|, |A∩B|/|B|, Poisson-bootstrap stability
(200 replicates; share of edges retained in ≥ 95 %), face validity of the three named
cases, H1a p (2,000 permutations, same design), NMI vs SuperClass, Louvain modularity
and its z against a small rewiring null. Reported as a whole surface, never a best cell.

## 8. Monte Carlo error propagation (objective d)

Sample compounds ∝ compounds-per-direct-parent; replace the direct parent by a
competitor `a ~ p(a|d)` restricted to the overlay; record the shallowest of
kingdom/superclass/class/subclass whose (padded) label slot changes, else `none`.
Nulls: wrong class uniform over the other direct-parent classes; and ∝ class size.
200 replicates (Poisson-resampled counts) × 20,000 compounds; 95 % percentile CIs.
Run on `G_c` (θ=4 forces slots 1–4 equal, so this is reported as a construction check),
on B-only (tree-free, the informative one), and on unfiltered S1.

## 9. Addendum — calibration outcome (appended after `10_build_overlay.py`, before any validation code ran)

Source: `results/phase2/calibration.json`, `overlay_meta.json`. No expert-flagged class
score had been computed when this section was written.

* Calibration set: 93 R5 pairs (69 direct parent–child) + 13 slot pairs from the 32
  PATH-violating shapes, of which 1 exists in S1 → **94 evaluable pairs**.
* Rule outcome: **τ* = 3, n_min* = 10** (recall 0.755, retain 7.0 %, lift 10.8; the
  n_min = 100/1000 cells fall below the 0.5 recall floor).
* Primary `G_c` = A(4) ∩ B(3, 10): **2,754 edges**, 1,069 non-isolated nodes,
  density 0.0048. |A| = 6,386, |B| = 73,184, Jaccard 0.036.
* **Plan §8 gate.** "NO-GO if A and B agree on < 10 % of either edge set." Observed
  |A∩B|/|A| = 43.1 %, |A∩B|/|B| = 3.8 %. Because |A|/|B| = 8.7 %, the B-side share cannot
  exceed 8.7 % even if A ⊂ B, so read literally the gate tests relative set size, not
  agreement; the overlap is 6.2× the independence expectation (|A||B|/|S1| = 447).
  Decision, taken blind to any validation result: **the primary stays A∩B** as frozen in
  §2; the gate's fallback A(θ=4) alone is added as a labelled secondary H1a row (same
  design, same seed). The literal gate reading is reported as triggered.
* **Size definition.** The pilot's r = 0.788 is reproduced only with size = total
  alternative-parent labels under `d` (r = 0.7877); with compounds-as-direct-parent it is
  0.669. The primary keeps compounds (§4). Decile matching on label totals is added as a
  labelled secondary, and is the definition used for the pilot-reproduction rows.
* Descriptive (not validation): 86.7 % of the 3,604 direct-parent classes have
  `H_f = 0` on `G_c`.
