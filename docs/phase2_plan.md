# Phase 2 — a weighted confusability overlay on ChemOnt

CLL798D Network Science · Group 2 · Pathway 1, extension component
Baseline: Djoumbou Feunang *et al.* (2016), *J Cheminform* **8**:61 (ClassyFire/ChemOnt).
Prerequisite: Phase 1a–1d, complete. Every figure below is either measured in this repo
(`results/scan.json`, `results/dp_alt_pairs.tsv`, `results/graph_data.json`) or quoted from
the pilot recorded in the Phase-2 brief. Numbers not yet measured are marked **TBD**.

---

## 0. Where Phase 1 leaves us, and the one-line gap

Phase 1 established that ChemOnt is exactly the object the paper describes — a strict tree,
4824 categories + root, max depth 11 — and that the paper's stated label invariants fail on
a small but structured minority of the 73,105,281 released rows:

| invariant | violating rows | source |
|---|---:|---|
| PATH (label is a root-path of its direct parent) | 198,959 | `scan.json.invariants` |
| — fully populated 5-slot labels only | 79,171 | `path_violation_rows_fully_populated` |
| R5 (alt parents unrelated to direct parent), label instances | 596,367 (upper bd.) | `scan.json` |
| — distinct violating (dp, alt) pairs | 93, of which **69 direct parent–child** | `scan.json` |
| R6 intermediate nodes, inclusive reading | 85,036 | `scan.json` |
| R6, strict Table S1 reading | 9,589 | `scan.json` |

Phase 1's verdict is *descriptive*: the invariants fail. It says nothing about **where in the
taxonomy** the failures concentrate, nor whether their location is predictable, nor whether
that location carries information a downstream consumer could act on. That is the gap
Phase 2 fills, and it is the gap that turns a topologically trivial object (a tree: no
cycles, no triangles, degree-assortativity dominated by the root, every centrality a
function of depth and fan-out) into a weighted graph on which the course's machinery is
non-vacuous.

---

## 1. Research question and falsifiable hypothesis

**Question.** ChemOnt forces every compound into exactly one direct parent, but ClassyFire
emits, alongside that choice, a mean of 12.6 *alternative parents* per compound
(922,418,764 labels over 73.1M rows). Do the alternative-parent co-assignment statistics
encode a recoverable, class-level signal of *taxonomic ambiguity* — a signal that is
(a) not merely a restatement of class size, (b) concentrated on classes an independent
expert panel actually got wrong, and (c) structured as a graph rather than as noise?

**Hypothesis H1 (falsifiable).** Define a weighted, undirected overlay graph
`G_c = (V, E_c, w)` on the ChemOnt class set from *competing* (not descriptive)
alternative-parent co-assignments (§3–4). Then:

* **H1a (validation).** The 20 distinct ChemOnt classes implicated in the paper's own
  expert-reviewed errors (Additional file 4: 17 false positives + 13 false negatives) score
  higher on a per-class ambiguity statistic derived from `G_c` than size-matched controls,
  at p < 0.05 under a size-stratified permutation null.
* **H1b (structure).** `G_c` is not explainable by a degree-preserving null: its weighted
  clustering coefficient and modularity under Louvain exceed the configuration-model
  expectation by more than 3 null standard deviations.
* **H1c (alignment).** The community partition of `G_c` is *not* a relabelling of the
  ChemOnt tree: NMI between Louvain communities on `G_c` and the SuperClass partition is
  bounded away from 1 (target: < 0.7), i.e. the overlay carries information the tree does
  not already encode.

**Falsification.** H1a fails if the size-matched permutation p-value exceeds 0.05 under the
pre-registered primary criterion of §3. H1b fails if the observed clustering/modularity lie
within the null's central 95 %. H1c fails if NMI ≈ 1, which would mean the overlay is a
weighted redrawing of the tree and the extension has produced nothing new. **All three are
reportable outcomes**; §6 says how a negative result is written up.

**What H1 does *not* claim.** See §2.

---

## 1A. The data this plan operates on (reference)

Placed early because none of what follows is intelligible without it. This section is
*only* about the data; the methodology resumes at §2. For the source-file formats, the
repeat-padding convention and the OBO stanza layout, see `docs/how_this_was_built.md`
Part I — this section covers only what is specific to Phase 2 and does not repeat it.

### 1A.1 `results/dp_alt_pairs.tsv` — the file the overlay is built from

Produced by the single 73.1M-row pass. Three columns, tab-separated, one header line:

```
direct_parent   alternative_parent   rows
2309            3940                 1878047
```

One row means: **across all 73,105,281 compounds, `rows` of them were assigned
`direct_parent` as their direct parent while also carrying `alternative_parent` in their
alternative-parent list.** It is a weighted, directed co-occurrence count between two
ChemOnt classes.

| Property | Value |
|---|---|
| Rows | **1,047,341** |
| Total label instances (`Σ rows`) | **922,418,764** |
| Rows with an empty `direct_parent` | **5**, carrying 35 label instances |
| Distinct undirected pairs after self-loop removal | **1,045,955** |

Resolved through `data/chemont_dictionary.tsv`, the heaviest pairs are:

| Label instances | direct parent | alternative parent |
|---:|---|---|
| 1,878,047 | Alpha amino acid amides | Organic oxides |
| 1,877,618 | Alpha amino acid amides | Hydrocarbon derivatives |
| 1,765,771 | Alpha amino acid amides | Carbonyl compounds |
| 1,351,480 | Alpha amino acid amides | Secondary carboxylic acid amides |
| 1,261,556 | Alpha amino acids and derivatives | Hydrocarbon derivatives |

**Read that table before designing anything.** The top pairs are `Organic oxides`,
`Hydrocarbon derivatives`, `Carbonyl compounds` — generic descriptors that co-occur with
almost everything. Raw co-occurrence weight is dominated by ubiquity, not by competition.
That single observation is why §3 exists.

### 1A.2 The other inputs

| File | Schema / contents | Role in Phase 2 |
|---|---|---|
| `results/path_counts.tsv` | `kingdom, superclass, class, subclass, direct_parent, rows` — 3,631 rows | Evidence source S2: the 79,171 rows whose label path is not a root-path of its own direct parent |
| `results/inode_sets.tsv` | `subclass, direct_parent, intermediate_nodes, rows` — 5,178 rows | Cross-check only; R6 violations are lineage errors, not competition |
| `results/scan.json` | All Phase 1b/1c figures + `sufficient_statistics` | Denominators; the 93 R5-violating pairs (evidence source S3) |
| `data/chemont_dictionary.tsv` | `numeric_id, chemont_id, name, parent_numeric_id, parent_name` — 4,824 categories + root | Resolves every numeric ID; defines the tree used for LCA and ancestry |
| `data/ChemOnt_2_1.obo` | The authors' own release, `sha256 8616a6ec…51fe22` | Second opinion on the tree. **Disagrees with the dictionary on 4 parent edges**, so per-level counts differ (OBO 765/1,729/2,297 vs dictionary 766/1,729/2,296). Always state which artefact a number came from |
| `data/paper_supplementary/…MOESM4_ESM.xlsx` | Sheet 1: 800 CIDs + SMILES. Sheet 2: 22,043 assignment rows with `FALSE POSITIVES` / `FALSE NEGATIVES` columns | **The only ground truth that exists.** See §1A.5 |

No rescan of the 2 GB source is required for anything in this plan.

### 1A.3 Why this problem is hard — the four label types

From the paper's Additional file 1, Table S1:

| Label | Paper's definition | Competing candidate? |
|---|---|---|
| **Direct parent** | "the category corresponding to the largest skeleton or most dominant feature" | This *is* the assignment |
| **Alternative parents** | "other categories that describe the classified compound and **do not display a parent-child relationship** to each other or to the direct parent" | **No — descriptors by design** |
| **Intermediate nodes** | "descendants of the subclass and ascendants of the direct parent" | No — lineage filler |
| **Substituents** | functional groups present, minus anything the classification already implies | No |

The extension rests on a claim the paper does not make: that a *minority* of alternative
parents are not descriptors but genuine rivals for the direct-parent slot. Nothing in the
data labels which is which. Establishing that separation is the scientific core (§3).

Note also that the paper **contradicts itself** on the R5 rule: the Results section says
alternative parents have no *ancestor–descendant* relation to the direct parent, while
Table S1 says the weaker *parent-child*. Phase 1c tested both readings; 69 of the 93
violating pairs are direct parent-child and therefore breach both.

### 1A.4 The graph object being built

| Property | Value | Consequence |
|---|---:|---|
| Classes acting as a **direct parent** | 3,604 | Entropy is only defined for these — not all 4,824 |
| Classes acting as an **alternative parent** | 3,172 | |
| **Union — the overlay's node set** | **3,843** | Not 4,824. 981 classes never appear in this file at all |
| Undirected edges, unfiltered | 1,045,955 | |
| **Density on the union** | **0.1417** | |

**A density of 0.14 is dense, not sparse.** Small-world statistics, betweenness centrality
and community detection are close to meaningless on a graph this dense — every node is
within a hop or two of everything. Filtering is not cosmetic; it is what makes any network
measure interpretable at all. §5 triages the measures on exactly this basis.

### 1A.5 Ancestry cannot be the separation criterion — measured

The intuitive filter is "keep pairs that are taxonomically close". Measuring the depth of
the lowest common ancestor of every pair shows how little room that leaves:

| LCA depth | Pairs | % of pairs | Label instances | % of mass |
|---|---:|---:|---:|---:|
| 0 (root) | 9,431 | 0.9 % | 348,129 | 0.0 % |
| **1 (Kingdom)** | **846,118** | **80.8 %** | **711,476,573** | **77.1 %** |
| 2 (SuperClass) | 153,460 | 14.7 % | 119,737,972 | 13.0 % |
| 3 (Class) | 31,896 | 3.0 % | 70,113,451 | 7.6 % |
| 4 (SubClass) | 4,789 | 0.5 % | 14,995,302 | 1.6 % |
| 5–9 | 1,642 | 0.2 % | 5,747,302 | 0.6 % |
| *siblings (same parent)* | *14,313* | *1.4 %* | *40,157,023* | *4.4 %* |

**Four-fifths of all pairs meet only at the Kingdom.** They share nothing but "both are
organic". And outright ancestry — the R5 violation — covers just **93 pairs out of
1,047,336**, so it cannot serve as a general criterion, though those 93 are a
high-confidence evidence stream of their own.

### 1A.6 The validation set, and a trap in it

Additional file 4, Sheet 2, holds the seven-expert review of the 800-compound test set:
**17 false positives** and **13 false negatives**, marked in two dedicated columns.

The resolution chain, exactly:

| Step | Count |
|---|---:|
| FP + FN markings | 17 + 13 = 30 |
| Distinct category **names** among them | 26 |
| Names resolving to a ChemOnt numeric ID | **25** |
| Of those, classes that ever act as a direct parent (so have an entropy score) | **20** |

The pilot's `n = 20` is that last row. Report the whole chain, not just the 20 — an
examiner asking "20 out of what?" should get an answer.

> **Naming trap.** Additional file 5 also involves the number 20 — but those are 20
> **compounds** in the ChEBI comparison (claim E9), a completely different set from the
> 20 **classes** here. The two are unrelated and must never be conflated in the report.

**This is the only ground truth available, and it is small.** Twenty points is a pilot,
not a result. §6 states what a negative result looks like and how the circularity risk —
these 20 classes come from the same pipeline that produced the labels being scored — is
handled.

---

## 2. Four caveats that constrain the design (stated before the method, not after)

**C1 — This measures ambiguity, not error.** There is no ground truth for 73.1M compounds.
The overlay weight is a measure of *how often ClassyFire hedged*, not of *how often
ClassyFire was wrong*. Any sentence of the form "edge weight = probability of
misclassification" is unsupportable and must not appear in the report. The defensible claim
is: *high-weight edges mark class pairs the system itself declines to separate cleanly*.
The operational value proposition — a downstream LLM or annotation pipeline reading the
weight as a "verify this one" flag — is a **use-case argument**, evaluated on whether the
flag concentrates on the only 30 known errors, not on a calibrated error rate.

**C2 — n = 20.** The only labelled validation set in existence is the 20 distinct ChemOnt
classes recoverable from Additional file 4's 17 FP + 13 FN. p = 0.0101 on 20 points is a
pilot, and will remain one. Three things can be done and none of them fixes it:
1. Report a confidence interval on the effect, not just p (bootstrap over the 20).
2. Add a *second, independent* validation axis (§6.2, ChEBI disagreement from Additional
   file 5 — 20 compounds, again small, but a *different* 20 and a different error mode).
3. Pre-register the primary statistic and the null before running, so the p-value is not
   the survivor of a search. **This is the main defence and it is free.**
   Everything else is secondary and labelled as such.
The honest framing for the report: *this is a hypothesis-generating result with one
pre-registered confirmatory test at n = 20, not an established effect.*

**C3 — Alternative parents are descriptive by design.** The paper states R5 explicitly:
alternative parents "describe the compound but do not have an ancestor–descendant
relationship ... with the Direct Parent". They are *additional true descriptors*, not
rejected hypotheses. Treating all 922M as competing hypotheses is simply wrong. The
measured evidence for how wrong: over all 1,047,336 mappable (dp, alt) pairs, the
lowest common ancestor sits at

| LCA depth | distinct pairs | share | label instances | share |
|---:|---:|---:|---:|---:|
| 0 (root only) | 9,431 | 0.9 % | 348,129 | 0.04 % |
| **1 (Kingdom)** | **846,118** | **80.8 %** | **711,476,573** | **77.1 %** |
| 2 (SuperClass) | 153,460 | 14.7 % | 119,737,972 | 13.0 % |
| 3 (Class) | 31,896 | 3.0 % | 70,113,451 | 7.6 % |
| 4 (SubClass) | 4,789 | 0.5 % | 14,995,302 | 1.6 % |
| 5–9 | 1,642 | 0.2 % | 5,747,302 | 0.6 % |

*(measured directly from `results/dp_alt_pairs.tsv` against `data/chemont_dictionary.tsv`;
self-pairs excluded, 0 unmappable ids.)* Four fifths of the mass is pairs that meet only at
Kingdom — a functional-group descriptor attached to a structurally unrelated skeleton.
Those are the paper's design intent and must be filtered out. Separating them from the
genuinely competing minority is §3, and it is the scientific core of this extension.

**C4 — It is an overlay.** `data/ChemOnt_2_1.obo`, `data/chemont_dictionary.tsv` and
`scripts/chemont.py` are frozen. No edge is added to, removed from, or reweighted in the
tree. The overlay is a *second* edge set over the same nodes, written to new files. The
Phase-1 test suite (`scripts/test_checks.py`, 31 tests) must still pass unchanged after
Phase 2; that is the mechanical guarantee of non-modification, and is a go/no-go gate.

**C5 — Circularity.** The 20 expert-flagged classes come from the same authors' pipeline
that emitted the alternative-parent labels being scored. Two distinct circularity risks:
* *Shared-substrate*: the panel reviewed ClassyFire output, so both signals derive from one
  system. Unfixable; must be stated. It bounds the claim to "the system's hedging predicts
  the system's reviewed errors", not "predicts chemical truth".
* *Same-instance leakage*: the 800 test compounds are (almost certainly) inside the 73.1M
  release, so the flagged classes' entropies are computed partly from the very rows the
  panel judged. **This is measurable and must be measured**: 800 rows out of 73.1M is
  ≈ 1.1 × 10⁻⁵ of the corpus, so leakage is numerically negligible, but the plan should
  *demonstrate* that rather than assert it — recompute the ambiguity score for the 20
  flagged classes with those compounds' InChIKeys excluded and show the score is unchanged
  to 3 decimals. Cost: one extra grouped pass over the sample, minutes.
* Partial mitigation for shared-substrate: Additional file 5 (ChEBI comparison, 20
  compounds) provides errors adjudicated against an *external* ontology. Weak (n = 20,
  different error definition) but genuinely independent of ClassyFire's own panel.

---

## 3. The descriptive-vs-competing separation

This is the central methodological section. The task: given a (direct parent `d`,
alternative parent `a`) pair with row count `n(d,a)`, decide whether `a` is a *descriptor*
of the compound (paper's intent, filter out) or a *competing candidate* for the direct-
parent slot (retain, weight).

Note first what **cannot** be used: ancestry. R5 already forbids an ancestor–descendant
relation, and Phase 1 found only **93** pairs that violate it. So "is `a` an ancestor of
`d`?" partitions 93 pairs out of 1,047,336 and is useless as a general criterion — though
those 93 (69 of them direct parent–child, dominated by
`Aryl thioethers → Alkylarylthioethers`, 551,743 rows) are a separate, high-confidence
evidence stream in their own right (§4, source S3).

### 3.1 Candidate criterion A — taxonomic proximity (structural)

`a` competes with `d` iff they are close in the tree: `depth(LCA(d,a)) ≥ θ`.
Rationale: a genuine alternative *subclass* assignment must share most of the skeleton
lineage; a descriptor need not. Pre-registered primary setting **θ = 4** (LCA at or below
SubClass), which retains 6,431 pairs / 20.7M label instances (0.6 % / 2.2 %). Strictest
variant, **sibling-only** (`parent(d) = parent(a)`): 14,313 pairs, 40,157,023 label
instances (1.4 % of pairs, 4.4 % of mass) — note siblings are *not* a subset of θ = 4,
since siblings at depth 3 have LCA depth 2, which is why both are reported.

* Pros: parameter-free given θ, computable in seconds from `dp_alt_pairs.tsv` +
  `chemont.py`, no distributional assumptions, interpretable to a chemist.
* Cons: θ is a knob; assumes competition is local, which fails for a genuinely misplaced
  compound (the `Depsipeptides` / `Peptidomimetics` case sits across a Class boundary).

### 3.2 Candidate criterion B — asymmetry / conditional surprise (statistical)

A descriptor co-occurs with *everything*; a competitor co-occurs specifically with `d`.
Score each ordered pair by pointwise mutual information against the marginal alternative-
parent frequency:

```
PMI(d, a) = log2 [ P(a | dp = d) / P(a) ]
P(a | dp = d) = n(d,a) / Σ_x n(d,x)          # rows of dp_alt_pairs.tsv grouped by d
P(a)          = Σ_x n(x,a) / Σ_{x,y} n(x,y)  # global marginal of a as an alternative parent
```

Retain the pair iff `PMI(d,a) ≥ τ` **and** `n(d,a) ≥ n_min` (a support floor; τ, n_min
pre-registered from the calibration in §3.4, not tuned on the validation set). A high-PMI
pair means "when the system lands on `d`, it reaches for `a` far more than chance" —
the operational reading of competition. Symmetrise for the overlay by
`PMI_sym = min(PMI(d,a), PMI(a,d))` where both directions exist, which suppresses the
one-sided descriptor case.

* Pros: purely data-driven; no appeal to the tree, so criterion A and B are *independent*
  evidence and their agreement is meaningful; naturally down-weights the ubiquitous
  descriptors ("Organooxygen compounds"-type nodes with huge marginals).
* Cons: two knobs; PMI is unstable in the low-count tail (hence `n_min`); the marginal `P(a)`
  is itself contaminated by whatever competition exists.

### 3.3 Choosing between them — and testing that the choice is not arbitrary

Do **not** choose by validation-set performance; that would spend the only labels we have on
model selection and invalidate §6. Choose by three prior criteria, in this order:

1. **Agreement.** Compute the Jaccard overlap of the retained edge sets A(θ=4) ∩ B(τ,n_min)
   over a grid of (θ, τ, n_min). If the two criteria — one structural, one statistical, with
   no shared assumption — agree on a large core, that core is the defensible edge set.
   **Pre-registered primary overlay = A ∩ B** (the intersection), with A-only and B-only
   reported as sensitivity variants. Rationale: intersection is conservative and its errors
   are of the "missed edge" kind, which weakens rather than manufactures H1a.
2. **Stability under resampling.** Bootstrap the corpus at the row level (Poisson-bootstrap
   the `n(d,a)` counts, 200 replicates — cheap, since it operates on the 1.05M-row
   sufficient statistic, not on 73.1M rows). Report the fraction of edges retained in ≥ 95 %
   of replicates. A criterion whose edge set is not stable is not usable.
3. **Face validity on known cases.** The overlay must contain, with high weight, the cases
   Phase 1 independently identified as genuine confusion:
   `Aryl thioethers`–`Alkylarylthioethers` (551,743 rows), the
   `Depsipeptides`/`Peptidomimetics` axis (58,345 rows), `Oxosteroids`–`3-oxosteroids`
   (6,218). If a criterion drops these, it is wrong regardless of its statistics. This is a
   *necessary* check, not a fitting procedure — the three cases were fixed by Phase 1
   before Phase 2 began.

**Arbitrariness test.** Sweep θ ∈ {2,3,4,5} × τ ∈ {1,2,3} × n_min ∈ {10,100,1000} (36
settings). For each, recompute the full H1a/H1b/H1c pipeline and tabulate the results.
Report the *whole surface*, not the best cell. The claim is only credible if the sign and
rough magnitude of the effect survive most of the surface; if it survives only near one
corner, say so and downgrade the conclusion to "the effect is contingent on the separation
rule, which we could not fix on independent grounds". That sentence is an acceptable
deliverable.

### 3.4 Calibration set (kept separate from validation)

Reserve the **93 R5-violating pairs** and the **32 distinct path-violation shapes** as a
*calibration* set for choosing τ and n_min: these are cases where Phase 1 already proved,
from the paper's own rules, that something is structurally wrong — independently of the
expert panel. Tuning on them and validating on the 20 expert classes keeps the two disjoint.
Document this split before any parameter is chosen.

---

## 4. Construction of the overlay

**Node set `V`.** All 4,824 ChemOnt categories (root excluded), keyed by the release's
numeric ids so the overlay joins to `results/graph_data.json` and `taxonomy_nodes.tsv`
without remapping. Measured facts to state in the report: only **3,604** classes ever occur
as a direct parent in the release, **3,172** occur as an alternative parent, their union is
**3,843**; 4,120 classes appear anywhere at all (`scan.json`). Isolated nodes are retained
in `V` — dropping them would inflate every density and centrality statistic.

**Evidence sources → edges.** Three sources, deliberately kept as separate weight channels
so their contributions are auditable rather than blended into one opaque score:

| src | evidence | file | edge | raw weight | scale |
|---|---|---|---|---|---|
| **S1** | direct-parent ↔ alternative-parent co-assignment, filtered by §3 | `results/dp_alt_pairs.tsv` | (d, a) | `n(d,a)` rows | 1,047,336 mappable pairs; 1,045,955 after self-loop removal + undirected collapse (matches `graph_data.json.meta.confus_edges_total`) |
| **S2** | alternative-parent **co-occurrence** (a, b both listed on the same compound) | `results/alt_cooccurrence_pairs.tsv` — **not yet generated** | (a, b) | rows | **TBD**; scan comment estimates ≈ 78 pairs/row × 73.1M, so expect 10⁷–10⁸ raw increments over ≤ 4824² = 23.3M possible pairs |
| **S3** | rule-violation edges: R5 pairs, R6 (intermediate node, direct parent) pairs, and PATH shape mismatches (reported vs expected slot) | `scan.json`, `path_counts.tsv`, `inode_sets.tsv` | as listed | rows | 93 R5 pairs; 16 top R6 pairs enumerated, 85,036 rows total; 32 path shapes |

S1 is the primary channel and the one the pilot used. S2 is the channel that makes the
object a genuine *class–class network* rather than a star around direct parents, and is the
one piece of new computation over the 73.1M rows. S3 is small, high-precision, and is used
for face validity (§3.3) and as a separate high-confidence edge tier — never merged into the
S1 weight, because its edges are *proved* violations while S1's are statistical.

**Weights and normalisation.** Raw counts span ≈ 6 orders of magnitude and are dominated by
class size, so three normalised forms are computed and all three are carried through the
analysis (with the primary pre-registered):

| form | definition | property | use |
|---|---|---|---|
| `w_raw` | `n(d,a)` | interpretable, size-dominated | reporting, edge tiering |
| `w_cond` | `n(d,a) / Σ_x n(d,x)` | row-stochastic; removes `d`'s size | asymmetric/directed views, PageRank |
| **`w_pmi`** (primary) | `max(0, PMI_sym(d,a))` from §3.2 | size-controlled both ends, additive in bits | centrality, clustering, Louvain |

`w_pmi` is the primary because §5's size confound (r = 0.788 between log₁₀ class size and
entropy) is the single largest threat to the result, and only a both-ends-normalised weight
addresses it at the *edge* level rather than only in the null model.

**Per-class ambiguity score.** Keep the pilot's definition for continuity —
`H(d) = −Σ_a p(a|d) log₂ p(a|d)` over the *filtered* edge set — and additionally report the
unfiltered version, since the pilot's numbers (flagged mean 5.81 bits vs 4.67 overall) are
unfiltered. Confirmed reproducible: the unfiltered per-class entropy over all 3,604
direct-parent classes has mean **4.671**, median 5.023, max 7.525, computed here from
`dp_alt_pairs.tsv` — identical to the pilot's 4.67 and to the `"e"` field already emitted by
`scripts/05_build_graph_data.py`. **The filtered score is a new quantity and the pilot's
p-values do not transfer to it**; they must be recomputed.

**Prototype already in the repo.** `scripts/05_build_graph_data.py` already builds an
undirected, unfiltered version of S1 (1,045,955 edges, top 14,000 shipped to the frontend)
and the unfiltered entropy. Phase 2 must *supersede* it in a new script rather than edit it,
so the Phase-1 frontend artefacts stay byte-reproducible.

---

## 5. Network analysis: what is informative here, and what would be ritual

The tree has 4,822 edges, zero triangles, and centralities that are pure functions of depth
and fan-out — every course measure is degenerate on it. The overlay has ~10⁶ weighted edges
over 3,843 non-isolated nodes; unfiltered undirected density is **0.1417**, which is *dense*,
not sparse — a fact that by itself disqualifies several standard analyses unless the §3
filter is applied first. Honest triage:

| measure | on the overlay | verdict |
|---|---|---|
| Degree / strength centrality | ranks classes by how many distinct classes they are confusable with | **informative**, but must be reported against class size or it re-derives size |
| Eigenvector centrality | finds the mutually-confusable core | **informative**; expect it to localise on the steroid / terpenoid / peptide blocks — a testable prediction |
| PageRank on `w_cond` (directed d→a) | "where does the system's uncertainty flow" | **informative**; the directed form is the natural one and the tree cannot supply it |
| Betweenness | 3,843 nodes × dense weighted graph = expensive, and the interpretation ("bridge between confusion regions") is thin | **ritual** unless the filtered graph is sparse; run only on the A∩B graph, else omit |
| Clustering / transitivity (weighted, Onnela) | does confusability come in cliques (whole neighbourhoods mutually ambiguous) or chains? | **informative and central** — this is H1b, and it is exactly the property a tree cannot have |
| Assortative mixing by **depth** and by **class size** | is confusion within-level or across-level? | **informative** — a strongly depth-assortative overlay says ambiguity is a *sibling* phenomenon; a disassortative one implicates the `Depsipeptides`-type cross-level failures |
| Assortative mixing by kingdom / SuperClass | 2 kingdoms, 31 SuperClasses; a categorical mixing matrix with a modularity-like `r` | **informative**, cheap, and directly interpretable |
| Erdős–Rényi `G(n,m)` null | wrong null: does not preserve the enormous degree heterogeneity (out-degree min 1, median 151, mean 290.6, max 1,749) | **ritual** — include one line showing it fails, then discard it, and say why |
| Configuration model / degree-preserving rewiring | preserves the heterogeneity that ER destroys | **the correct null**; see §6.1 |
| Small-world (σ, ω vs rewired lattice) | at density 0.14 the unfiltered graph has L ≈ 2 by construction | **ritual on the unfiltered graph**; meaningful only on the filtered one, and only if that is sparse (**TBD** after §3) |
| Power-law fit to the weighted degree distribution | `powerlaw`/Clauset–Shalizi–Newman with a likelihood-ratio test against log-normal and exponential | **run it, and expect it to fail** — report the LR test honestly; "the degree distribution is heavy-tailed but log-normal is not rejected in favour of a power law" is a legitimate and more likely result than a scale-free claim. Do not report an α without `x_min`, its uncertainty, and the LR comparison |
| Modularity + Louvain, and spectral bisection | H1c: do confusion communities recut the taxonomy? | **the payoff analysis** — compare the partition against SuperClass via NMI / ARI, and inspect the blocks that cut across SuperClasses; spectral bisection of the Fiedler vector gives an independent, deterministic second opinion on the top-level split |

All of it in NetworkX; 3,843 nodes is small enough that only betweenness and the Louvain
sweep need care.

---

## 6. Nulls, significance, and validation

### 6.1 Structural nulls (for H1b/H1c)

* **Degree-preserving rewiring** (`nx.double_edge_swap`, ≥ 10·|E| swaps, 200 replicates) on
  the binarised filtered graph. Reports z-scores for transitivity, mean clustering, and
  Louvain modularity. Preserves the degree sequence, which is the confound that ER misses.
* **Weighted null**: for weighted clustering and weighted modularity, rewire the *topology*
  as above and reshuffle weights over the rewired edges (preserving the weight
  distribution). State plainly that this is a heuristic weighted configuration model, not an
  exact uniform sample.
* **Size-preserving edge null**: sample each edge with probability ∝ `size(u)·size(v)`
  (a Chung–Lu / gravity null on class size). This tests the specific alternative "the
  overlay is just a function of how many compounds each class has", which is the alternative
  the r = 0.788 correlation makes most dangerous.
* **Report the null's own spread**, not just the point z-score; give the empirical p as
  `(1 + #{null ≥ obs}) / (1 + R)`.

### 6.2 Validation (H1a) — and its weaknesses, stated

**Primary test (pre-registered, one test).** Size-stratified permutation. Recover the 20
distinct ChemOnt classes from Additional file 4 (columns for false positives and false
negatives on worksheet 2, already parsed by `scripts/03_verify_evaluation.py`, which reads
them **by category name** — the name→numeric-id join against
`data/chemont_dictionary.tsv` is a real failure point and must be asserted 20/20 or the
misses listed explicitly). Compute the filtered ambiguity score for each. Draw 10,000
control sets, each control matched to its flagged class by class-size decile. p-value =
fraction of control sets whose mean score ≥ observed. Pilot (unfiltered): **p = 0.0101**,
against **p = 0.0001** uncontrolled — the gap between the two is the size confound made
visible, and both must be reported side by side.

**Reported alongside, always:** flagged mean 5.81 bits vs 4.67 overall = 73.7th percentile,
and a bootstrap CI on the effect size (mean difference in bits), because with n = 20 the CI
is the honest summary and the p-value is not.

**Positive controls.** `Sesquiterpenoids` and `Diterpenoids` — both expert-flagged, both at
the 99.9th percentile in the pilot — should remain extreme under the filtered score. If the
§3 filter destroys them, the filter is suspect.

**Secondary, independent axis.** Additional file 5's 20-compound ChEBI comparison
(reproduced in Phase 1d: ~31 ClassyFire categories vs ~33 manual ChEBI classes, ~14 missing
from ChEBI). Define per-compound ChEBI/ClassyFire disagreement and test its correlation with
the mean overlay ambiguity of that compound's assigned classes. Independent of the authors'
own panel, therefore partially escapes C5. **Small (n = 20, and a *different* 20 from the
flagged-class count — the numerical coincidence must not be allowed to blur the two in the
writeup), and the disagreement definition is our construction, so this is exploratory.**

**Leakage check** (from C5): recompute the 20 flagged classes' scores excluding the 800 test
compounds' InChIKeys. Expect no change beyond the third decimal; report the delta.

**Named weaknesses, to appear in the report verbatim:**
1. n = 20; one pre-registered test; no held-out set; no replication cohort exists.
2. Shared substrate (C5) bounds the claim to within-system consistency.
3. The 20 classes are themselves size-biased (expert review finds errors where compounds
   are); the size-matched permutation addresses this but the matching is by decile, which is
   coarse — report the size distributions of flagged vs matched controls so a reader can
   judge the residual imbalance.
4. `H` is one of several possible ambiguity scores; report at least one alternative
   (max-alternative share `1 − max_a p(a|d)`, i.e. Gini/collision-based) and show whether
   the conclusion depends on the choice.

### 6.3 What a negative result looks like, and how it is reported

| outcome | reported as |
|---|---|
| H1a fails (p > 0.05 filtered, though the unfiltered pilot passed) | "the association with expert-flagged classes is not robust to the descriptive/competing separation" — a substantive finding: it means the pilot's signal was carried by *descriptive* alternative parents, i.e. by class size and label verbosity, not by competition. Report the filtered and unfiltered p-values, the size-matched and uncontrolled p-values, in one 2×2 table. |
| H1b fails | the overlay is a configuration-model graph: heterogeneous degrees and nothing more. Report the z-scores. Still a valid course deliverable — the null-model machinery is exactly the syllabus content. |
| H1c fails (NMI ≈ 1) | the overlay is a weighted redrawing of the tree; the extension yields no new partition. Report the NMI and the confusion matrix against SuperClass. |
| the §3 sweep is unstable | report the full 36-cell surface and downgrade to "contingent on the separation rule" (§3.3). |

No outcome is permitted to be dropped. The pre-registration file (§7, `docs/phase2_prereg.md`)
is written and committed *before* the validation script is run, and is cited in the report.

---

## 7. Deliverables mapped to the repo

Runtime anchor: the full 73.1M-row pass takes **395.2 s ≈ 6.6 min** at 184,987 rows/s
(`scan.json`). Everything that works off the sufficient statistics is seconds to minutes.

| # | new file | reads | writes | est. runtime |
|---|---|---|---|---|
| 0 | `docs/phase2_prereg.md` | — | — | — (written first, committed before step 4) |
| 1 | *(re-run, no new code)* `scripts/02_scan_dataset.py --pairs` | 2 GB dataset via `$ATLAS_DATA` | `results/alt_cooccurrence_pairs.tsv`, refreshed `scan.json` (fills `distinct_alt_cooccurrence_pairs`, currently `null`) | **~8–12 min** — the 6.6 min baseline plus ~78 pair increments/row; **memory is the risk**, up to 23.3M Counter keys (§9 R1) |
| 2 | `scripts/10_build_overlay.py` | `dp_alt_pairs.tsv`, `alt_cooccurrence_pairs.tsv`, `scan.json`, `path_counts.tsv`, `inode_sets.tsv`, `chemont.py` | `results/overlay_edges.tsv` (u, v, w_raw, w_cond, w_pmi, src, lca_depth, sibling), `results/overlay_nodes.tsv` (id, name, depth, superclass, size, H_filtered, H_unfiltered, degree, strength) | ~60–120 s |
| 3 | `scripts/11_separate_edges.py` | `results/overlay_edges.tsv` | `results/separation_sweep.tsv` (36 cells × retained-edge count, Jaccard(A,B), bootstrap stability, face-validity flags) | ~5–10 min (200 Poisson-bootstrap replicates on 1.05M rows) |
| 4 | `scripts/12_validate_flagged.py` | `data/paper_supplementary/…MOESM4_ESM.xlsx`, `…MOESM5_ESM.xlsx`, `overlay_nodes.tsv`, `chemont_dictionary.tsv` | `results/validation.json` (name→id join audit, primary p, uncontrolled p, bootstrap CI, leakage delta, ChEBI secondary) | ~30–60 s (10,000 permutations over 3,604 classes) |
| 5 | `scripts/13_network_measures.py` | `results/overlay_*.tsv` | `results/network_measures.json` (§5 table, all measures + null z-scores), `results/communities.tsv` | ~10–20 min (betweenness and the 200-replicate rewiring dominate; cap betweenness with `k`-sampling and say so) |
| 6 | `scripts/14_figures.py` | all of the above | `results/figs/*.png` (see §10) | ~2 min |
| 7 | extend `scripts/test_checks.py` | — | — | +L1 tests: PMI on a toy taxonomy; LCA-depth on the toy tree; permutation-test correctness on a synthetic set with a planted effect; **assertion that no Phase-1 output changes** |
| 8 | `docs/phase2_results.md` | — | — | the ledger, in the style of `docs/paper_claims.md`: every Phase-2 number with the file it came from |

New dependencies: `networkx`, `scipy`, `powerlaw`, `openpyxl` (already listed;
**note `openpyxl` is not currently installed in the working environment** — step 4 will fail
until it is). Add to `requirements.txt`; do not vendor. Add a `phase2` target to the
`Makefile` mirroring the existing `data` target's ordering.

**Frontend (optional, only if time remains after step 6).** `results/network_layout.json`
and the existing WebGL renderer already draw the top 14,000 unfiltered S1 edges. A minimal,
high-value addition is a toggle between "all alternative-parent edges" and "competing edges
only (A∩B)", which visualises §3 directly. Explicitly **out of scope** for the scientific
result and must not consume time budgeted for §6.

---

## 8. Schedule and go/no-go gates

| week | work | checkpoint | go/no-go |
|---|---|---|---|
| 1 | §7 step 0 (pre-registration) + step 1 (`--pairs` run) | `alt_cooccurrence_pairs.tsv` exists; `distinct_alt_cooccurrence_pairs` non-null; row totals reconcile to 73,105,281 | **NO-GO if the run OOMs or exceeds 30 min** → fall back to S1+S3 only, drop S2, and say so in the report. The extension survives without S2. |
| 2 | steps 2–3: overlay + separation sweep | `separation_sweep.tsv` complete; A∩B retains the three face-validity cases; ≥ 90 % of A∩B edges stable across bootstrap replicates | **NO-GO if A and B agree on < 10 % of either edge set** → the two criteria measure different things; report that as the finding, pick A(θ=4) as primary on grounds of interpretability, and reframe §3 as an open question |
| 3 | step 4: validation, run **once**, after pre-registration is committed | `validation.json`; name→id join is 20/20 | proceed regardless of sign — §6.3 covers both branches. **Gate on process, not on result**: if the join is < 20/20, resolve the names before running the test, not after seeing it |
| 4 | step 5: network measures + nulls | `network_measures.json`, `communities.tsv` | **NO-GO on betweenness only** if the filtered graph is still dense (density > 0.05) → omit betweenness and small-world, per §5 |
| 5 | steps 6–8: figures, results ledger, report draft; prelim presentation | draft report + slides | — |
| 6 | sensitivity surface, writeup hardening, final presentation | full 36-cell table in an appendix | — |

Hard rule: **step 4 is run exactly once on the primary configuration.** The 36-cell sweep is
run afterwards and reported as sensitivity, clearly labelled as such. Violating this rule
would convert p = 0.01 at n = 20 into nothing at all.

---

## 9. Risks

| id | risk | likelihood | impact | mitigation |
|---|---|---|---|---|
| R1 | `--pairs` OOM: ~78 pairs/row over 73.1M rows, up to 23.3M distinct keys in a Python `Counter` | medium | blocks S2 | Two-pass: first pass counts only pairs among the top-K most frequent alternative parents (K ≈ 500, covers most mass); or shard by `hash(a,b) % 8` and merge. Gate at week 1; S2 is droppable. |
| R2 | The §3 separation is arbitrary and the result depends on it | medium | undermines the central claim | §3.3 agreement + stability + face validity; the 36-cell surface; the pre-registered A∩B primary; §6.3 permits "contingent" as an outcome |
| R3 | Size confound survives normalisation (r = 0.788 is large) | medium | H1a becomes uninterpretable | `w_pmi` normalises at both ends; size-matched permutation; Chung–Lu size null; report the flagged-vs-control size distributions |
| R4 | Name→id join on Additional file 4 fails for some of the 20 classes (synonyms, renames between the 2016 supplement and the release dictionary) | **high** | shrinks n below 20 | Resolve via `chemont_dictionary.tsv` names, then OBO synonyms (`data/ChemOnt_2_1.obo`), then manual adjudication documented per class. Report the join table in full. Resolve *before* the test is run |
| R5 | `openpyxl` missing in the environment; `zstd` CLI needed for step 1 | high / certain | steps 1 and 4 fail | `pip install -r requirements.txt`; confirm `zstd` on PATH before week 1 |
| R6 | Circularity objection at viva (C5) | certain | credibility | pre-empt it in §2 of the report, quantify the leakage, present the ChEBI axis as the partial answer, and concede the shared-substrate limit rather than defending it |
| R7 | Scope creep into the frontend | medium | eats §6 time | frontend explicitly gated behind step 6 completion |
| R8 | Accidental modification of Phase-1 artefacts | low | breaks reproduction | new scripts only (`10_`–`14_`); `test_checks.py` must pass unchanged; verify Phase-1 outputs byte-identical before submission |
| R9 | A negative result is read as project failure | medium | grade | §6.3 pre-commits the reporting of every branch; the null-model and community machinery is the syllabus deliverable and is produced either way |

---

## 10. Report and presentations

**Final report** (structure; the extension section is the load-bearing half):
1. Reproduction summary — one page, pointing at `docs/paper_claims.md` and the Phase-1c
   invariant table above. Do not re-litigate Phase 1.
2. Motivation: the tree/uncertainty mismatch; 922,418,764 alternative-parent labels;
   the 79,171 fully populated path violations; the `Aryl thioethers` case.
3. **Caveats first** (§2, all five). Placing them before the method, not in a limitations
   paragraph at the end, is the single strongest signal of rigour available here.
4. Overlay construction (§4), with the three-source table and the weight table.
5. **The separation** (§3) — longest section, with the LCA-depth table, both criteria, the
   agreement/stability analysis and the 36-cell surface.
6. Network analysis (§5) with the informative/ritual triage stated openly, including the ER
   null shown to be the wrong null.
7. Validation (§6) with the pre-registration cited and the outcome reported per §6.3.
8. `docs/phase2_results.md` as an appendix ledger, one row per number, each with its source file.

**Figures** (six, no more): (F1) tree vs overlay side by side, same layout, showing that the
tree has no triangles and the overlay does; (F2) the LCA-depth mass distribution (§2's table
as a bar chart) with the θ cut marked; (F3) entropy vs log₁₀ class size scatter, r = 0.788,
with the 20 flagged classes highlighted — this figure *shows the confound* rather than
hiding it; (F4) the size-matched permutation null with the observed mean marked; (F5)
observed vs configuration-model clustering and modularity; (F6) Louvain communities coloured
onto the radial tree layout, with cross-SuperClass blocks called out.
Chart conventions per the repo's existing visual style.

**Prelim presentation** (≈ week 5): the question, the LCA-depth table (why 80.8 % of the
evidence must be thrown away), the pilot's two p-values side by side (0.0001 uncontrolled vs
0.0101 size-matched) as the honest headline, and the pre-registration commitment. Ask the
audience for the §3 criterion they find least arbitrary — that is a genuinely useful question
and it demonstrates the methodological point.

**Final presentation**: F1, F3, F4, F6; one slide on what a negative result would have looked
like and whether we got one; one slide on the n = 20 ceiling and what dataset would lift it
(an independently adjudicated, multi-thousand-compound ChemOnt gold set — which does not
exist, and saying so is the correct answer).

---

## 11. Discrepancies found between the Phase-2 brief and the repository

Recorded here because Phase 1's standard is that discrepancies are quantified, not smoothed.

1. **`README.md` Phase-1c point 3 prints 79,164** fully populated path violations, while
   `results/scan.json` gives `path_violation_rows_fully_populated = 79171` — the value the
   README's own next sentence says is emitted by the scan. 79,164 = 198,959 − 119,795 is
   precisely the mixed-population derivation the README warns against. The README figure
   should be corrected to 79,171 (a Phase-1 documentation fix, not a Phase-2 change).
2. **`dp_alt_pairs.tsv` sums to 922,418,729 alternative-parent labels, 35 short of
   `scan.json`'s 922,418,764.** The gap is exactly 5 pair rows carrying an empty
   `direct_parent` field (7 rows in the corpus have a null direct parent). Any Phase-2
   statistic recomputed from the TSV will be 35 labels short; this is negligible but must be
   handled explicitly rather than discovered later.
3. **The file has 1,047,341 data rows, of which 1,047,336 are mappable**; after self-loop
   removal and undirected collapse it yields **1,045,955** edges, which matches
   `graph_data.json.meta.confus_edges_total` exactly.
4. **The pilot is partly already in the repo.** `scripts/05_build_graph_data.py` computes
   the same per-class entropy (verified: mean 4.671 over 3,604 classes, matching the brief's
   4.67) and an unfiltered undirected overlay, shipping entropy as node field `"e"` and the
   top 14,000 edges to the frontend. The permutation tests are not in the repo. Phase 2
   should supersede this in new scripts rather than edit it (§4, §7).
5. **README says the scan took 397 s; `scan.json` says 395.2 s.** Immaterial, noted for the
   ledger.
6. **The "20" collision.** The 20 distinct classes from Additional file 4's 17 FP + 13 FN is
   a different set from `E9-n = 20`, the ChEBI-comparison *compounds* in Additional file 5.
   Both appear in this plan; they must never be conflated in the writeup.
7. **`openpyxl` is declared in `requirements.txt` but not installed** in the current
   environment, and `pdftotext` is unavailable, so `rules.pdf` could not be read directly
   for this plan; the course requirements used here are those restated in the brief and the
   README. Confirm the deliverable list against `rules.pdf` before the prelim.
