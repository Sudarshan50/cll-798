# Phase 2 — results ledger

Implements `docs/phase2_plan.md` §3–§7 under the pre-registration `docs/phase2_prereg.md`.
Every number below is copied from a file in `results/phase2/`, named in the **source**
column. Nothing in Phase 1 was modified (verified by hash, §T).

**Pre-registration was frozen by hash, not by commit.** The user has not authorised git
commits. Sections 1–8 of the prereg were written first: sha256 `9197d853…cf4a7c` at
2026-11-12T06:18:02Z (`prereg_hash_initial.txt`). §9 (calibration outcome, gate decision)
was appended after `10_build_overlay.py` and before any validation code ran. The whole
file's sha256 `acbe48ec…a827c` was computed in the same process immediately before the
primary test (run 2026-11-12T06:26:34Z) and is stored in `primary_test.json` and
`validation.json.prereg`. The file is unchanged since then (`unchanged_since_primary: true`).
The primary test was run **once**; re-runs of `12_` reload it and assert that the score
vector's hash is unchanged.

## Verdicts

| hypothesis | pre-registered criterion | observed | verdict |
|---|---|---|---|
| **H1a** validation | size-matched permutation p < 0.05 | **p = 0.453** (T = 0.268 vs null 0.266 ± 0.142, n = 21) | **not supported** |
| **H1b** structure | z > 3 for weighted clustering **and** modularity vs degree-preserving null | z = 181.5 and 134.2 (tree-constrained null: 13.3 and 19.5) | **supported** |
| **H1c** alignment | NMI(Louvain, SuperClass) < 0.7 | NMI = 0.638; 0 of 171 communities cross a SuperClass | rule says *supported*; reported as **inconclusive** (§N3) |

In one line: the overlay is strongly structured (H1b). But the pre-registered competing
overlay does **not** pick out the expert-flagged classes (H1a). Its community structure is
nested in the tree by construction (H1c). The pilot's unfiltered signal reproduces
(p = 0.0102 vs the brief's 0.0101). It does not survive the θ = 4 filter, which zeroes out
86.7 % of classes, including both positive controls. The tree-free B-only score separates
the flagged classes strongly (p = 0.0001), but that comparison was not pre-registered and
is reported as exploratory.

---

## O. Overlay construction (`10_build_overlay.py`)

| # | quantity | value | source |
|---|---|---:|---|
| O1 | S1 directed (dp, alt) rows used / labels | 1,047,336 / 922,418,729 | `overlay_meta.json` s1 |
| O2 | rows dropped (empty direct parent) | 5 rows, 35 labels | plan §11.2 |
| O3 | S1 undirected pairs; pairs observed in both directions | 1,045,955; 1,381 | s1 |
| O4 | direct-parent classes / alt classes / union | 3,604 / 3,172 / 3,843 | s1 |
| O5 | unfiltered density on the union | 0.1417 | s1 |
| O6 | LCA-depth table (directed rows): root / K / SC / C / SubC / 5–11 | 9,431 / 846,118 / 153,460 / 31,896 / 4,789 / 1,642 pairs — matches plan §1A.5 exactly | `lca_depth_directed` |
| O7 | siblings (directed / undirected) | 14,313 / 14,193 | same |
| O8 | unfiltered entropy over 3,604 dp classes: mean / median / max | 4.671 / 5.023 / 7.525 (= pilot, = `05_`'s `"e"`) | `entropy` |
| O9 | r(log10 size, H_unfiltered), size = compounds as dp | 0.669 | `entropy` |
| O10 | same, size = alt-label total Σₓ n(d,x) | 0.7877 (= pilot's 0.788) | computed from `overlay_nodes.tsv` |
| O11 | S2 alt–alt pairs / co-occurrence increments | 1,098,113 / 6,204,984,080 | `s2` |
| O12 | S2 input check | `results/alt_cooccurrence_pairs.tsv` (sha256 `5b17cd4c…`) is identical, pair for pair, to C − n(a,b) − n(b,a) from the final `class_projection.npz` | verified in session |
| O13 | S2 competing (A∩B) pairs; not in S1 primary | 2,953; 1,053 | `s2` |
| O14 | S3 tier pairs: R5 / R6 / PATH | 121: 93 / 16 / 13 | `s3` |

## S. Separation (`10_` calibration, `11_separate_edges.py`)

| # | quantity | value | source |
|---|---|---:|---|
| S1 | calibration set: R5 pairs (direct parent–child) | 93 (69) — matches `scan.json` | `calibration.json` |
| S2 | PATH-violating shapes; non-null slot pairs; present in S1 | 32; 13; 1 | same |
| S3 | evaluable calibration pairs | 94 | same |
| S4 | **chosen τ*, n_min*** (max lift, recall ≥ 0.5) | **3, 10** (recall 0.755, retain 7.0 %, lift 10.8) | same |
| S5 | \|A(4)\|, \|B(3,10)\|, **\|A∩B\|** | 6,386; 73,184; **2,754** (0.26 % of pairs, 16,298,289 labels = 1.77 % of mass) | `separation_sweep.tsv` primary row |
| S6 | Jaccard(A,B); share of A; share of B; enrichment over independence | 0.036; 0.431; 0.038; 6.2× | same |
| S7 | plan §8 gate "agree on < 10 % of either set" | triggered on the B side (3.8 %). This is bounded by \|A\|/\|B\| = 8.7 % even if A ⊂ B. A∩B kept as primary per prereg §2; A-only is reported as a secondary (V9) | prereg §9 |
| S8 | bootstrap stability (200 Poisson replicates): edges kept in ≥ 95 % | **90.2 %** (gate ≥ 90 %: pass); mean retention 97.4 % | same |
| S9 | cells with stable95 ≥ 0.9 | 29 / 36 | `separation_summary.json` |
| S10 | face validity, Aryl thioethers–Alkylarylthioethers | present; 22.9th pct by w_pmi, 99.9th by w_raw | sweep |
| S11 | face validity, Oxosteroids–3-oxosteroids | present; 82.7th pct by w_pmi, 87.0th by w_raw | sweep |
| S12 | face validity, Depsipeptides/Peptidomimetics axis | **absent from S1 entirely** (0 rows for any of its pairs), so no S1 criterion can retain it; present in the S3 tier | `separation_summary.json` |
| S13 | θ = 5 cells | lose both in-S1 face cases (fail face validity) | sweep |

The full 36-cell surface, with the H1a and H1b/H1c columns merged in, is
`results/phase2/sweep_full.tsv`.

## V. Validation — H1a (`12_validate_flagged.py`)

**Join chain** (`validation.json.join_chain`, full table in `join_audit`):
30 FP + FN markings (17 + 13) → 26 distinct names → 25 exact + 1 case-insensitive
(`Transition Metal Hydrides` → 536 `Transition metal hydrides`) = **26/26 resolved** →
**21 direct-parent classes** (primary). The exact-join subset is the plan's n = 20. The
5 resolved classes that never act as a direct parent are Amino acids, peptides, and analogues;
Carbohydrates and carbohydrate conjugates; Lipids and lipid-like molecules; Organic
zwitterions; Prenol lipids. (These are flagged *classes*. They are unrelated to the 20
ChEBI-comparison *compounds* of V14.)

| # | test | T_obs | null mean ± sd | p | source |
|---|---|---:|---:|---:|---|
| **V1** | **PRIMARY: H_f on G_c, n = 21, size-decile matched, R = 10,000** | **0.268** | 0.266 ± 0.142 | **0.453** | `primary` |
| V2 | H_f, uncontrolled | 0.268 | 0.184 ± 0.119 | 0.225 | `secondary` |
| V3 | H_unfiltered, matched (compound deciles) | 5.638 | 5.354 ± 0.201 | 0.069 | same |
| V4 | H_unfiltered, uncontrolled | 5.638 | 4.664 ± 0.321 | 0.0004 | same |
| V5 | **pilot definition**: n = 20, H_unfiltered, label-total deciles | 5.813 | 5.464 ± 0.165 | **0.0102** (brief 0.0101) | same |
| V6 | pilot definition, uncontrolled | 5.813 | 4.661 ± 0.328 | 0.0002 (brief 0.0001) | same |
| V7 | n = 20, H_f matched / uncontrolled | 0.281 | 0.273 / 0.183 | 0.438 / 0.198 | same |
| V8 | alternative score 1 − maxₐ p(a\|d) on G_c, matched / uncontrolled | 0.070 | 0.073 / 0.053 | 0.497 / 0.281 | same |
| V9 | gate fallback A(4)-only, matched / uncontrolled | 0.415 | 0.400 / 0.293 | 0.437 / 0.201 | same |
| V10 | **exploratory**: tree-free B(3,10)-only, matched | 3.817 | 2.890 ± 0.227 | **0.0001** | same |
| V11 | H_f, label-total deciles | 0.268 | 0.253 | 0.417 | same |

| # | quantity | value | source |
|---|---|---|---|
| V12 | effect (bits) with 95 % bootstrap CI: filtered / unfiltered | +0.002 [−0.226, 0.331] / +0.281 [−0.125, 0.642] | `descriptive` |
| V13 | flagged mean vs overall: H_unfiltered (n = 21; n = 20) | 5.64 vs 4.67 (70.0th pct); 5.81 (76.3rd pct; the brief says 73.7th, percentile definition not stated there) | same |
| V14 | flagged classes with any H_f > 0 | 4 / 21 (overall 13.3 %): Dipeptides 3.01, Hexoses 1.07, Ketoximes 0.93, Aromatic monoterpenoids 0.62 | `join_audit` |
| V15 | **positive controls** Sesquiterpenoids / Diterpenoids | H_unfiltered at the 100.0th / 99.9th pct; **H_f = 0 for both**. By plan §6.2 the filter is therefore suspect | `positive_controls` |
| V16 | residual size imbalance (log10 compounds, quartiles) | flagged [2.98, 3.84, 4.62] vs matched pool [2.98, 3.86, 4.63] vs all [1.82, 2.78, 3.77]; 10 of 21 flagged sit in deciles 8–9 | `size_balance` |
| V17 | **leakage** (C5): 765/800 test compounds found by CID and subtracted (8,616 labels) | max \|ΔH_f\| over flagged = 0.00026, ΔT = 1.3e-5, unchanged to 3 dp | `leakage` |
| V18 | ChEBI axis (Additional file 5, 20 compounds, 18 found in release) | Spearman(disagreement, H_f of dp) = 0.02 (p = 0.93); with H_unfiltered 0.17 (p = 0.50) | `chebi_secondary` |
| V19 | 36-cell sensitivity surface, H1a | 0 / 36 cells with p < 0.05 (min 0.096 at θ=2, τ=2, n_min=10); T_obs > null mean in 25 / 36; θ = 5 cells have no flagged class with H_f > 0 (p = 1) | `sensitivity_sweep_h1a` |

**§6.3 2×2 (the table the plan requires for a failed H1a):**

| | size-matched | uncontrolled |
|---|---:|---:|
| filtered (G_c, n = 21) | 0.453 | 0.225 |
| unfiltered (n = 21, compound deciles) | 0.069 | 0.0004 |
| unfiltered, pilot definition (n = 20, label deciles) | 0.0102 | 0.0002 |

Reading, per §6.3: *the association with expert-flagged classes is not robust to the
descriptive/competing separation*. Even unfiltered, it weakens from 0.0102 to 0.069 when
the matching uses compounds instead of label totals and the 21st resolved class is
included. The pilot's signal is largely carried by descriptive alternative parents, that
is, by class size and label verbosity. The one exception is V10. The PMI criterion alone,
with no tree constraint, gives p = 0.0001 with a large effect. So the structural θ = 4 cut
is what destroys the signal, not the down-weighting of descriptors. This is a
post-hoc observation from a secondary row. It is a hypothesis for a future pre-registered
test, not a finding.

## N. Network measures — H1b, H1c (`13_network_measures.py`)

Primary graph G_c: 1,069 non-isolated nodes (of 4,824 retained in V), 2,754 edges,
density 0.0048, 171 components, largest component 62 nodes (benzoic acids).
Source: `network_measures.json`.

### N1. Observed vs nulls (z = (obs − mean)/sd)

| measure | observed | degree-preserving (R=200) | tree-constrained (R=100) | Erdős–Rényi (R=50) | Chung–Lu label freq (R=50) |
|---|---:|---:|---:|---:|---:|
| weighted clustering (Onnela) | 0.137 | 0.004 ± 0.001, **z 181.5** | 0.111 ± 0.002, z 13.3 | 0.001, z 351 | 0.086, z 10.4 |
| Louvain modularity | 0.954 | 0.501 ± 0.003, **z 134.2** | 0.936 ± 0.001, z 19.5 | 0.548, z 159 | 0.304, z 123 |
| transitivity | 0.473 | 0.019, z 309 | 0.447, z 7.6 | 0.005, z 575 | 0.196, z 60 |
| degree assortativity | 0.474 | −0.015, z 25.5 | 0.470, z 0.5 | — | −0.276, z 89 |
| max degree / degree CV | 33 / 1.10 | preserved | preserved | 13.9 / 0.44 | 218.8 / 2.36 |

ER is the wrong null, and the table shows why. It cannot reproduce the degree spread
(max degree 13.9 against 33, CV 0.44 against 1.10). It is shown once and then discarded.
Chung–Lu on label frequency fails in the opposite direction: size-driven hubs reach degree
219. Size alone therefore does not produce this graph. The tree-constrained null keeps
degrees and confines every swap to one SubClass block, so the θ constraint holds. Even
then clustering and modularity exceed it (z 13.3, 19.5). H1b therefore holds beyond
what the filter itself imposes.

### N2. Centrality, assortativity, other measures

| # | measure | value |
|---|---|---|
| N2a | strength top-3; Spearman vs log10 label freq | Aminosaccharides, O-glycosyl compounds, Oligopeptides; ρ = 0.08 |
| N2b | eigenvector (largest component) top-3 | Halobenzoic acids, Aminobenzoic acids, M-phthalic acid and derivatives; ρ = 0.10 |
| N2c | PageRank on directed w_cond top-3 | Monosaccharides, Amino acids, Aldehydes; ρ = 0.42 |
| N2d | betweenness (density 0.0048 < 0.05, so computed exactly) top-3 | Oligopeptides, Benzoic acids, 3-halobenzoic acids and derivatives; ρ = 0.39 |
| N2e | assortativity: degree / depth / log label freq | 0.474 / 0.382 / 0.290 |
| N2f | within-SuperClass / within-SubClass share; SuperClass assortativity | 1.0 / 1.0 / 1.0, all **by construction** (θ = 4) |
| N2g | same-depth / sibling / parent–child edge shares; mean \|Δdepth\| | 0.462 / 0.331 / 0.014; 0.72. Ambiguity is a sibling-and-cousin phenomenon, not cross-level |
| N2h | small world (largest component only): L, C vs rewired; σ | 2.17, 0.526 vs 2.00, 0.314; σ = 1.54 |
| N2i | power law, strength: α, x_min, n_tail; LR vs lognormal (p); vs exponential (p) | 2.43 ± 0.06, 28.6, 490/1,069; R = −5.49 (p = 4e-8); R = −6.19 (p = 6e-10). **Lognormal favoured; no scale-free claim** |
| N2j | power law, degree | α 1.85 ± 0.03, x_min 2; R vs lognormal −8.18 (p = 3e-16). Lognormal favoured |
| N2k | spectral bisection (Fiedler, largest component) | 32 / 30 split, both inside Benzoic acids and derivatives; ARI vs Louvain 0 (Louvain keeps that component whole) |

### N3. H1c in detail

NMI(Louvain, SuperClass) = 0.638, ARI = 0.196, 171 communities. The pre-registered rule
(NMI < 0.7) returns *supported*. But **0 of 171 communities cross a SuperClass** and
NMI against the SubClass partition is **0.983**. Criterion A(θ) admits only pairs whose
LCA is at depth ≥ θ. For every θ ≥ 2 (all 36 cells), no edge can leave a SuperClass. The
partition is therefore a refinement of the tree by construction, and NMI < 0.7 measures
granularity, not re-cutting. Across the surface, NMI runs from 0.636 to 0.984. H1c is
reported as **inconclusive**: the pre-registered test cannot fail in the direction it
was meant to detect. This is a flaw in the plan's design of H1c, recorded rather than
patched.

The only informative version is the tree-free B(3,10)-only overlay: 3,478 nodes,
73,184 edges, 28.5 % of edges within a SuperClass. It gives NMI = 0.383, ARI = 0.216,
and **7 of 7 communities cross SuperClasses**. Example: a 976-node block mixing
Organoheterocyclics, Organic acids and Organometallics. It also beats its own
degree-preserving null (z_wclust 330, z_Q 191, R = 10). On that graph the overlay
re-cuts the taxonomy, but this is secondary.

Variant S1 ∪ S2 (A∩B on both channels): 1,170 nodes, 3,807 edges, NMI 0.632, z_wclust
131, z_Q 175 (R = 20). Adding S2 changes nothing qualitative.

Sweep (H1b columns, R = 5 rewires per cell): every one of 36 cells has z_Q ≥ 17.6 and
z_wclust ≥ 12.7. H1b is not contingent on the separation rule.

## M. Monte Carlo error propagation (`13b_error_propagation.py`)

200 replicates, each with Poisson-resampled counts and sizes, × 20,000 compounds; mean and
95 % CI. Levels: shallowest slot that changes when the direct parent is swapped
(`error_propagation.json`).

| overlay / null | kingdom | superclass | class | subclass | none |
|---|---:|---:|---:|---:|---:|
| **B-only competitor** (tree-free) | 0.0005 | **0.548** [0.541, 0.554] | 0.202 [0.197, 0.207] | 0.126 [0.122, 0.131] | 0.123 [0.118, 0.128] |
| uniform wrong class | 0.118 | 0.741 | 0.115 | 0.022 | 0.004 |
| size-proportional wrong class | 0.002 | 0.744 | 0.164 | 0.073 | 0.016 |
| unfiltered S1 (every alt parent) | 0.001 | 0.777 | 0.129 | 0.070 | 0.022 |
| primary A∩B | 0 | 0 | 0 | 0 | 1.000 (by construction) |

* Coverage: 99.75 % of compounds have a direct parent with ≥ 1 B-competitor; 31.1 % have
  one on A∩B.
* Reading: a swap to a B-competitor still changes the SuperClass in 54.8 % of draws.
  The uniform null gives 85.9 % for SuperClass-or-above, the size-proportional null 74.6 %.
  The competitor stays inside the SubClass (subclass or none) 24.9 % of the time, against
  2.6 % uniform and 8.9 % size-proportional. Competing alternatives are tree-local well
  beyond chance, yet most of them would still move the label at SuperClass level.
* The exact (non-MC) expectations agree with the MC means to < 0.001
  (`expected_exact`). The CIs reflect Monte Carlo and Poisson noise only. They are tight
  because the corpus is large, and they say nothing about whether a competitor is right
  (C1).

## T. Tests and non-modification

* `scripts/test_phase2.py`: 26/26 pass. Covers PMI against hand values, LCA depth and
  padded slots on a toy tree, PMI_sym and undirected collapse, filtered entropy, NMI/ARI,
  a permutation test that recovers a planted +1.5 effect (and stays quiet without one,
  while the uncontrolled version is fooled by size), H0 calibration, and bootstrap CI
  coverage.
* Phase-1 hashes before vs after (`phase1_hashes_before.txt`, checked by L3 of
  `test_phase2.py`): `scan.json`, `path_counts.tsv`, `taxonomy_claims.json` and
  `graph_data.json` are byte-identical.
* `python3 scripts/test_checks.py`: **31 passed, 0 failed** (L1, L2 and L3 all ran)
  after the full Phase-2 pipeline. The file itself is unmodified.

## Deviations from the plan, and why

1. **Primary n = 21, not 20.** The plan's own R4 resolution chain resolves
   `Transition Metal Hydrides` case-insensitively. This was fixed in the prereg before the
   test; n = 20 is reported (V5–V7).
2. **PMI_sym for one-directional pairs** uses the single direction. Requiring both would
   exclude every parent→child R5 pair, including the face-validity cases.
3. **§8 agreement gate** is literally triggered on the B side. A∩B was kept, as frozen in
   the prereg; A-only was added as a secondary row (S7, V9). Decided before validation.
4. **Size definition.** The pilot's r = 0.788 uses alternative-label totals, not compounds.
   The primary keeps compounds, as pre-registered. Label-total matching appears as
   secondary rows.
5. **H1c** is reported as *inconclusive* although the frozen rule says *supported*
   (§N3). Both are stated.
6. **Error propagation** headline uses the B-only overlay. On A∩B the outcome is 100 %
   "none" by construction. That result is still reported.
7. **Extra nulls and variants**, not in the plan: tree-constrained rewiring, B-only and
   S1 ∪ S2 graphs, and the 36-cell H1b/H1c columns with only R = 5 rewires per cell (the
   primary null keeps R = 200).
8. **S2 was not dropped.** It came from the parallel bipartite scan (see O12), not from
   `02_scan_dataset.py --pairs`. It is a sensitivity channel only and never enters H1a.
9. **Plan file names** `results/overlay_*.tsv` etc. are written under `results/phase2/`
   (ownership rule). No `14_figures.py`: the frontend payload `ui.json` carries the chart
   data. The `Makefile` and `requirements.txt` were not touched.
10. **Leakage** used the release's CID column in one 46 s `zstd | rg` pass. The 800 test
    compounds are listed by CID, not InChIKey; 35 were not found.

## Caveats (plan §2), as they bear on these results

* **C1 — ambiguity, not error.** Every weight here measures how often ClassyFire hedged.
  No edge weight is a misclassification probability. The Monte Carlo profile (M) shows
  where a label *would* move under a swap, not how often it is wrong.
* **C2 — n = 21 (20).** One pre-registered test, which failed. The CI on the filtered
  effect, [−0.23, +0.33] bits, spans zero comfortably. V10 is a hypothesis-generating
  observation at the same n.
* **C3 — alternative parents are descriptors by design.** 77.1 % of label mass meets only
  at Kingdom (O6). The calibration recovered 75.5 % of known-structural pairs while
  retaining 7 % of pairs (S4). Removing descriptors is what V3→V1 measures.
* **C4 — overlay only.** No tree edge was touched. Phase-1 outputs are byte-identical and
  `test_checks.py` is unchanged.
* **C5 — circularity.** Same-instance leakage is measured and negligible (V17,
  |ΔH_f| ≤ 0.0003). Shared-substrate circularity is unfixable. The independent ChEBI axis
  shows no association (V18, n = 18 compounds, exploratory).

## How to reproduce

```bash
cd 798_d_op
.venv/bin/python scripts/12_validate_flagged.py --extract   # once: 2 GB pass, ~1 min (needs zstd, rg, $ATLAS_DATA or ../chemont_project/data)
.venv/bin/python scripts/10_build_overlay.py        #   ~7 s  overlay_edges/nodes.tsv, calibration.json
.venv/bin/python scripts/11_separate_edges.py       #  ~25 s  separation_sweep.tsv (200 bootstrap reps)
.venv/bin/python scripts/12_validate_flagged.py     #  ~30 s  validation.json (primary reloaded from primary_test.json)
.venv/bin/python scripts/13_network_measures.py     # ~4.5 min network_measures.json, communities.tsv
.venv/bin/python scripts/13b_error_propagation.py   #  ~37 s  error_propagation.json
.venv/bin/python scripts/13c_build_ui.py            #   ~1 s  ui.json, sweep_full.tsv
.venv/bin/python scripts/test_phase2.py             #   ~2 s  26 tests
python3 scripts/test_checks.py                      # Phase-1 suite, must still pass
```

Total ≈ 6 min on an 8-core, 8 GB Mac. Deleting `results/phase2/primary_test.json` would
re-run the primary test. Doing so would violate the run-once rule (plan §8), and the new
run would not be the pre-registered one.
