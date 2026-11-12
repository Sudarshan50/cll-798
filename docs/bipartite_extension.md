# Extension — the bipartite nature of ClassyFire and its one-mode projections

CLL798D · Group 2 · requested by the course instructor: *explore the bipartite nature and
the one-mode projection*. Textbook frame: Newman, *Networks: An Introduction*, ch. 6
("Bipartite networks"), and the chapters on affiliation networks and their projections.

**The short version.**

1. The release is an affiliation network: 73,105,281 compounds × 3,851 ChemOnt categories,
   with 995,524,038 incidences.
2. Its compound-mode projection is almost complete (density 0.9917), and the clustering is
   manufactured by construction. Both are caused by a handful of hub labels.
3. The class-mode projection has density 0.188, and its weights mostly measure how common
   each label is. A maximum-entropy null that fixes both degree sequences (BiCM), plus an
   effect-size floor (lift ≥ 2), keeps 12.5% of its edges. This is the primary backbone.
4. Those edges carry real structure. Clustering is 2.8× the degree-preserving null.
   Communities agree only moderately with the taxonomy's own SuperClasses (NMI 0.41).
   Under a bipartite null that does not reward projection cliques (Guimerà M_B), the
   co-assignment partition scores 0.074, against 0.042 for the official SuperClass partition
   and 0.018 for the Class partition: the data's own grouping beats the tree's.
5. Functional-group profiles track the tree only weakly (Mantel r = −0.20). They are not
   nested beyond what their degrees imply. Hazard-relevant groups (nitro, azo, diazo,
   N-nitroso, perchlorate salts) are spread across hundreds of classes: knowing a
   compound's skeleton class does not tell you whether it carries a reactive group.

Every number below is read from `results/bipartite/bipartite_results.json` (abbreviated
**R**; the section key is given as `R:A_B…`, `R:C…`, `R:D…`) unless another file is named.
Inputs are the sufficient statistics of `scripts/09_scan_bipartite.py`
(`results/bipartite/*.npz`, `bipartite_scan.json`), `results/path_counts.tsv`,
`results/dp_alt_pairs.tsv`, `data/chemont_dictionary.tsv` and the 200k-row sample.

---

## 0. The bipartite objects

| symbol | rows × columns | definition | incidences |
|---|---|---|---:|
| **B** | compounds × categories | M(r) = {direct parent} ∪ alternative parents | 995,524,038 |
| **S** | compounds × substituents | the release's `substituents` list | 679,126,726 |
| **D** | direct-parent classes × substituents | D = Pᵀ S, P = one-hot direct parent | 679,126,698 |

Projections: C = BᵀB (class mode; diagonal = k_c), W = Bᵀ diag(1/(|M_r|−1)) B (Newman 2001
collaboration weights), F = SᵀS (substituent mode), and B Bᵀ (compound mode), which is never
materialised.

**Sum identities**, asserted in `14_bipartite_projections.py`:

- Σ_c k_c = Σ_L L·deg_M[L] = 995,524,038.
- Σ_s k_s = Σ_L L·deg_S[L] = 679,126,726.
- Both compound histograms sum to 73,105,281 rows.

Σ D is 28 short of Σ k_s. Those 28 incidences are exactly the substituents on the 7 rows with
no direct parent (`R:D…matrix.gap_vs_kS_sum`).

## A. Two-mode structure (`14_bipartite_projections.py`, `R:A_B….two_mode`)

| mode | n | mean degree | median | max |
|---|---:|---:|---:|---:|
| compounds in B (\|M\|) | 73,105,281 | 13.618 | 13 | 87 |
| categories in B (k_c) | 3,851 | 258,511 | 2,025 | 72,763,536 |
| compounds in S | 73,105,281 (245,989 have none) | 9.290 | 9 | 37 |
| substituents in S | 396 | 1,714,966 | 5,510 | 72,779,463 |

The compound side is narrow and unimodal. The category side spans eight decades, from 1 up
to 99.5% of all compounds. The top hubs of B, as a share of compounds:

| hub | share |
|---|---:|
| `Hydrocarbon derivatives` | 99.5% |
| `Organic oxides` | 73.3% |
| `Azacyclic compounds` | 64.3% |
| `Heteroaromatic compounds` | 49.7% |
| `Carbonyl compounds` | 40.3% |

**Power-law fit** of k_c, following Clauset, Shalizi & Newman (2009), with the `powerlaw`
package, discrete:

- Fit: x_min = 221,276, α = 1.73 ± 0.04, KS D = 0.050.
- The tail is only 403 of 3,851 classes (10.5%).

Likelihood-ratio tests (normalised R; R > 0 favours the power law):

| alternative | R | p | verdict |
|---|---:|---:|---|
| lognormal | −1.93 | 0.053 | lognormal marginally preferred |
| exponential | +5.49 | 4e-8 | power law preferred |
| truncated power law (nested LR) | −6.98 | 2e-4 | truncated preferred |

**Verdict: not a power law.** The distribution is heavy-tailed, but the tail is better
described as lognormal or truncated. The 396 substituents behave the same way: α = 1.41,
lognormal p = 0.15 (inconclusive), truncated preferred at p = 0.002. This is a negative
result and is reported as such. Hub labels are a design feature of ChemOnt (the generic
"chemical entities" layer), not evidence of preferential attachment.

## B. The compound-mode projection (`R:A_B….compound_projection_by_level`, `…_multilabel`)

### B.1 One level of the hierarchy: a disjoint union of cliques, in closed form

At a single level every compound carries exactly one label, so the projection is ⋃_c K_{n_c}.
With N = 73,105,281 nodes and class sizes n_c taken from `path_counts.tsv`:

- edges E = Σ_c C(n_c, 2), and density = E / C(N, 2);
- a compound in class c has degree n_c − 1, so that degree occurs with probability n_c/N;
- clustering is 1 for n_c ≥ 3 and 0 otherwise, so the average clustering is
  Σ_{n_c≥3} n_c / N, and transitivity is 1;
- components = number of labels + number of null rows;
- degree assortativity is exactly 1: both ends of every edge sit in the same clique.

Rows with a null slot at a level are isolated nodes. `scripts/test_bipartite.py` checks
every one of these formulas against networkx on an explicit toy graph.

| level | cliques | edges | density | mean degree | giant share | null rows (isolated) |
|---|---:|---:|---:|---:|---:|---:|
| kingdom | 2 | 2.666e15 | 0.9976 | 72.9M | 0.999 | 0 |
| superclass | 31 | 6.892e14 | 0.2579 | 18.9M | 0.364 | 53 |
| class | 641 | 2.522e14 | 0.0944 | 6.90M | 0.260 | 119,788 |
| subclass | 1,921 | 5.764e13 | 0.0216 | 1.58M | 0.108 | 119,795 |
| direct parent | 3,617 | 1.412e13 | 0.0053 | 386k | 0.026 | 7 |

At the subclass level, 11,368,062 rows have a direct parent at depth 3. Their subclass slot
repeats the class id (repeat-padding, `how_this_was_built.md` §B.2), and they are kept as
their own clique keyed by that id. The giant component at the superclass level is
`Organoheterocyclic compounds` (26.6M compounds). Average clustering is ≥ 0.998 at every
level, and that value is put in by construction, not measured.

### B.2 The multi-label projection (compounds linked if they share any label in M)

**Multigraph edges, exact:** Σ_c C(k_c, 2) = **8,683,247,409,594,716**. That is 3.25× the
number of distinct compound pairs (2.672e15), so a typical connected pair shares several
labels.

**Simple-graph density.** This has no closed form because pairs share labels in correlated
ways, so it was estimated on the 200k sample.

- **Density: 0.99168**, from 999,997 uniformly random pairs; 99% Wilson CI
  [0.99144, 0.99191].
- **Transitivity: 0.9996**, from 999,979 random triples: P(jk | ij, ik), which is unbiased
  for 3·triangles / connected triples.
- **Degree:** for 1,000 random focal compounds, the median compound is adjacent to 99.6% of
  the sample, and the 1st percentile to 99.5%.

The multi-label projection is therefore slightly *less* dense than the Kingdom projection
(0.9976), because "Organic compounds" itself is not a member of M.

**Hub dependence** (`R:…multilabel.hub_removal`; chart `ui.charts.hub_removal`). Labels were
removed in order of k_c, and the density and transitivity re-estimated on the same pairs and
triples:

| top-k labels removed | 0 | 1 | 2 | 3 | 5 | 10 | 20 | 50 | 100 | 200 |
|---|---|---|---|---|---|---|---|---|---|---|
| density | 0.992 | 0.828 | 0.728 | 0.653 | 0.497 | 0.325 | 0.219 | 0.077 | 0.025 | 0.006 |
| transitivity | 1.000 | 0.913 | 0.857 | 0.792 | 0.671 | 0.546 | 0.461 | 0.402 | 0.408 | 0.45 (71 triples) |

- 16.5% of all compound–compound edges exist *only* through `Hydrocarbon derivatives`, and
  67% only through the top 10 labels (`edges_only_via_top`).
- After removing 100 hubs, density is 0.025 but transitivity is still 0.41, 16× the density.
  In a random graph of that density, clustering would equal the density. That gap is the
  manufactured clustering Newman describes: every label is a clique.
- Removing 500 hubs leaves only 2 sampled connected triples, so transitivity is reported as
  null there.
- The substituent projection of compounds is equally saturated (density 0.9916, 99% CI
  0.9914–0.9919).

### B.3 Structural equivalence: what the projection forgets (`R:….structural_equivalence`)

`degrees.npz` holds the multiplicities of identical label sets M.

- **24,264,828 distinct label sets** describe 73,105,281 compounds, a compression of 3.0×.
- 14,745,280 sets are singletons, so **58,360,001 compounds (79.8%) share their entire label
  set with at least one other compound**.
- These compounds are structurally equivalent, with identical neighbourhoods in BBᵀ. No
  projection-based measure can tell them apart.
- The largest set is shared by 120,678 compounds, and 959 sets each cover ≥ 1,000 compounds
  (2.69M compounds in total). Chart: `ui.charts.signature_sizes`.

The projection also loses *which* label created an edge. Combined with density 0.99, the
compound projection contains essentially no information beyond the bipartite degree
sequence. That is the textbook argument, now quantified. **Analyses should be done on B,
not on BBᵀ.**

## C. The class-mode projection and validated backbones (`15_bipartite_backbones.py`, `R:C…`)

**Raw projection** (off-diagonal C, `R:C….raw_projection`):

- 3,851 nodes and 1,393,402 edges: density 0.188, mean degree 724.
- Average clustering 0.835, transitivity 0.606.
- Weights run from a median of 17 up to 53,530,353, for `Organic oxides`–`Hydrocarbon
  derivatives`.
- Onnela weighted clustering is 6.0e-6. It normalises by the maximum weight, so on this
  projection it is uninformative.

Strength and weight histograms are in the same JSON. The heaviest edges are all hub–hub
pairs (`heaviest_raw_edges`).

### Methods

**1. Hypergeometric SVN** (Tumminello et al. 2011, PLoS ONE 6:e17994).

- Test: p_ab = P(X ≥ C_ab), with X ~ Hypergeom(N_rows = 73,105,281, k_a, k_b).
- Tests: m = K(K−1)/2 = 7,413,175, so the Bonferroni cut is p < 0.01/m = 1.35e-9. A
  BH-FDR variant at 0.01 is also reported.
- Over-expression ratio: C_ab / E[C_ab], with E[C_ab] = k_a k_b / N.

`scipy.stats.hypergeom.logsf` takes about 1 ms per pair at this N, so the tail is computed
directly instead (`bipartite_common.hypergeom_logsf`). It takes the exact log-pmf from
`gammaln`, sums 64 terms of the pmf ratio recursion, and bounds the remainder geometrically,
which is valid because the pmf is log-concave. The result is a conservative upper bound on
p. It matches scipy to 1e-6 relative on every significant toy pair and is never below it
(test).

**2. BiCM** (Saracco et al. 2015, Sci Rep 5:10595; projection test Saracco et al. 2017,
New J Phys 19:053022). This is the null behind the primary backbone (method 2b).

- Model: p_rc = x_r y_c / (1 + x_r y_c), with Σ_c p_rc = d_r and Σ_r p_rc = k_c.
- Rows with equal degree share x_r, so the fit needs only the deg_M histogram: 81 degree
  classes × 3,851 categories instead of 73M rows.
- Solver: for fixed x, each column equation is monotone in log y_c and is solved by
  safeguarded Newton. The 80 free row parameters (one is a gauge) are then solved by
  `scipy.optimize.root` (hybr).
- Constraints hold to 5.8e-11 (rows) and 3.6e-14 (columns) relative, asserted < 1e-6.
  `max p_rc` = 0.99974.
- Test statistic: V_ab = Σ_r B_ra B_rb is Poisson-binomial with q_r = p_ra p_rb.
- **Poisson approximation**, λ = Σ_r q_r. Its upper tail dominates the Poisson-binomial
  tail, so it is conservative (tested against the exact convolution). BH-FDR at 0.01 over m.
- Sensitivity check: a refined normal approximation (Volkova 1996) with exact mean,
  variance and skewness. It validates 234,218 edges, Jaccard 0.95 with the Poisson set.

**3. Disparity filter** (Serrano, Boguñá & Vespignani 2009, PNAS 106:6483), on the
Newman-weighted W.

- α_ij = (1 − w_ij/s_i)^(k_i − 1).
- An edge is kept if α < 0.01 from either endpoint.

**2b. BiCM + lift ≥ 2 — the PRIMARY backbone** (`extension_research.md` §2.6). An edge is
kept only if it is BiCM-validated (BH-FDR 0.01) *and* has lift = C_ab / E_BiCM[C_ab] ≥ 2.
At N = 7.3e7, a 0.4% excess on a hub pair is about 30 sd, so significance alone keeps edges
of no substantive size. `ui.edges`, `ui.nodes[].comm` and `ui.communities` all come from
this backbone. Plain BiCM is kept as a row in the table.

**For every backbone**:

- Sparse clustering (average and transitivity).
- Louvain on the unweighted backbone (networkx `louvain_communities`, seed 42, resolution 1)
  and its modularity Q.
- For the primary backbone only, a second Louvain run **weighted by ln(lift)**. The weight
  is always > ln 2, so it is positive. The Newman weight was not used because it re-imports
  the size effect the lift removes.
- **Bipartite modularity M_B** (Guimerà, Sales-Pardo & Amaral 2007, PRE 76:036102, Eq. 5),
  evaluated on the full off-diagonal C, never on a backbone:

  M_B(P) = Σ_s [ Σ_{a≠b∈s} C_ab / Σ_d n_d d(d−1)  −  Σ_{a≠b∈s} k_a k_b / (Σ_d n_d d)² ]

  - The sums run over ordered pairs; n_d is the number of compounds with \|M\| = d.
  - The null is bipartite: the expected number of ordered co-assignments is k_a k_b scaled
    by the compound-side pair count. Unlike Newman–Girvan Q on C, it therefore does not
    count the clique each compound induces as community structure.
  - Nodes absent from a backbone are singleton modules and contribute 0.
  - `bc.bimod`. A toy test checks it on a block-diagonal bipartite graph: the true partition
    beats all 200 shuffles, one module gives Σk²/(Σk)², and all-singletons gives 0.
- The modularity of the SuperClass partition itself on the same graph.
- NMI (arithmetic normalisation) and ARI against the SuperClass (depth 2) and Class
  (depth 3) partitions, with each node mapped to its ancestor at that depth; shallower
  nodes map to themselves.
- A degree-preserving rewiring null: 10 replicates, one accepted double-edge swap per edge
  each. The swap is networkx's algorithm, reimplemented on arrays for speed (tested).

### Results — the backbone table (`R:C….backbones`; `ui.backbones`)

| backbone | edges | density | ⟨C⟩ | ⟨C⟩ rewired (10 reps) | Louvain Q | comms | NMI SC | ARI SC | NMI Class | ARI Class | **M_B** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| raw C | 1,393,402 | 0.188 | 0.835 | — | — | — | — | — | — | — | — |
| SVN Bonferroni | 201,779 | 0.0292 | 0.458 | 0.209 ± 0.001 | 0.300 | 7 | 0.402 | 0.268 | 0.387 | 0.080 | 0.0809 |
| SVN BH-FDR | 260,078 | 0.0357 | 0.438 | 0.214 ± 0.001 | 0.287 | 6 | 0.393 | 0.294 | 0.358 | 0.063 | 0.0754 |
| BiCM (FDR) | 222,927 | 0.0308 | 0.405 | 0.174 ± 0.001 | 0.302 | 8 | 0.421 | 0.340 | 0.401 | 0.093 | 0.0875 |
| **BiCM + lift ≥ 2 (primary)** | **174,525** | **0.0242** | **0.336** | **0.120 ± 0.001** | **0.335** | **8** | **0.412** | **0.322** | **0.394** | **0.080** | **0.0739** |
| primary, Louvain weighted ln(lift) | 174,525 | 0.0242 | 0.336 | 0.120 | 0.434 (weighted) | 8 | 0.404 | 0.252 | 0.402 | 0.083 | 0.0674 |
| disparity (W) | 94,736 | 0.0187 | 0.807 | 0.750 ± 0.001 | 0.223 | 7 | 0.313 | 0.198 | 0.338 | 0.063 | 0.0658 |

**What the table shows:**

- Every backbone is a single connected component. The primary has 3,799 nodes, and its
  clustering z is 383.
- On the primary backbone, the SuperClass partition reaches Newman–Girvan modularity only
  0.124, against Louvain's 0.335. The weighted run gives 0.191 vs 0.434.
- The weighted and unweighted partitions of the primary agree at NMI 0.65: same scale,
  partly different boundaries.
- Primary communities by dominant SuperClass (`ui.communities`):

  | community | size | purity |
  |---|---:|---:|
  | Organoheterocyclic | 900 | 0.63 |
  | Lipids | 660 | 0.59 |
  | Organoheterocyclic (mixed) | 577 | 0.19 |
  | Mixed metal/non-metal | 424 | 0.81 |
  | Phenylpropanoids | 380 | 0.54 |
  | Nucleosides | 353 | 0.31 |
  | Benzenoids | 300 | 0.62 |
  | Organic acids | 205 | 0.42 |

**What the lift floor removes** (`R:C….lift_floor`):

- It removes 48,402 of the 222,927 BiCM edges (21.7%), but 70.0% of their total weight.
  The removed edges join big classes (median endpoint k = 255,782, against 26,149 for kept
  edges).
- **Hub–hub pairs essentially disappear.** Among the top 20 labels by k_c there are 186 raw
  pairs. BiCM validates 73 of them, the primary keeps 3:

  | surviving hub–hub pair | lift |
  |---|---:|
  | `Aryl chlorides`–`Organochlorides` | 5.3 |
  | `Alkyl aryl ethers`–`Phenoxy compounds` | 3.3 |
  | `Trialkylamines`–`Aralkylamines` | 2.8 |

  All three are synonym-like pairs, not plumbing.
- Of the 15 heaviest raw edges, BiCM validates 12 and the primary validates 0.
  `Organic oxides`–`Hydrocarbon derivatives` (lift 1.004) is gone.

### Bipartite modularity: does co-assignment beat the tree? (`R:C….taxonomy_bimod`; `ui.taxonomy_bimod`)

| partition of the 3,851 categories | modules | M_B | shuffled-label null (100 reps) | z |
|---|---:|---:|---|---:|
| Kingdom | 2 | 0.0179 | — | — |
| **SuperClass** (official tree) | 31 | **0.0424** | 0.0018 ± 0.0022 | 18.6 |
| Class (official tree) | 641 | 0.0182 | — | — |
| **primary Louvain** (BiCM + lift ≥ 2) | 8 | **0.0739** | 0.0033 ± 0.0016 | 43.3 |
| primary, weighted Louvain | 8 | 0.0674 | — | — |
| plain-BiCM Louvain | 8 | 0.0875 | — | — |

- **Yes.** Every backbone's Louvain partition beats the official SuperClass partition: 0.066–0.088
  against 0.042, and the primary by 1.7×. The Class partition scores lower than SuperClass.
  The 641 classes are too fine, and they split co-assigned families.
- Both partitions are far above shuffled labels, so the tree *is* aligned with
  co-assignment, just less than the data's own grouping.
- Absolute M_B values are small, because the hub labels co-occur with everything and their
  co-assignments necessarily cross any module boundary.
- **Honest caveat.** Plain BiCM's partition scores *higher* than the primary (0.0875 vs
  0.0739). M_B is computed on the full weighted C, so a partition that keeps the heavy
  hub–hub edges inside a module collects their weight. The lift floor deliberately drops
  those edges. M_B rewards capturing mass, not chemistry. Louvain was not run on the M_B
  modularity matrix itself (`extension_research.md` §3, option 2), so these are lower bounds
  on the best attainable M_B.

**Edge-set Jaccard** (`R:C….jaccard`):

|  | SVN-B | SVN-FDR | BiCM | BiCM+lift | disparity | dp–alt overlay |
|---|---:|---:|---:|---:|---:|---:|
| SVN-B | 1 | 0.776 | 0.806 | 0.625 | 0.253 | 0.166 |
| SVN-FDR | | 1 | 0.857 | 0.671 | 0.211 | 0.206 |
| BiCM | | | 1 | 0.783 | 0.203 | 0.177 |
| BiCM+lift | | | | 1 | 0.153 | 0.136 |
| disparity | | | | | 1 | 0.084 |

**The Phase-2 dp–alt overlay** (`results/dp_alt_pairs.tsv`) is a *different* projection:

- It is directed (direct parent → alternative parent).
- It omits alternative-parent–alternative-parent co-assignment.
- Here it is compared undirected, with self-pairs dropped: 1.05M pairs, all of them present
  in C.
- 85% of BiCM edges are dp–alt pairs, but the BiCM keeps only 18% of the overlay
  (`dp_alt_overlap`). For the primary backbone the figures are 84% and 14%.

### What is genuine and what is artefact in C

- **Genuine.** Validated backbones are more clustered than their own degree sequence
  implies: 2.8× for the primary (z = 383), 2.3× for plain BiCM, 2.2× for SVN (z ≈ 170–230).
  They resolve into 6–8 communities that line up with the major SuperClasses (NMI
  0.39–0.42, ARI 0.27–0.34). Under the bipartite null those communities beat the official
  SuperClass partition (M_B 0.074 vs 0.042).
- **Genuine.** The strongest BiCM edges are *synonym-like* pairs that label the same
  substructure twice:

  | pair | co-occurring compounds | obs/exp |
  |---|---:|---:|
  | `Aryl chlorides`–`Organochlorides` | 10.5M | 5.3 |
  | `Aryl fluorides`–`Fluorobenzenes` | | 7.1 |
  | `Aryl bromides`–`Organobromides` | | 11.8 |
  | `Anisoles`–`Methoxybenzenes` | | 7.2 |

  This is real redundancy in ChemOnt's alternative-parent layer, the kind of thing the
  Phase-2 confusability analysis targets.
- **Artefact or design.** Agreement with the tree is moderate, not high. Co-assignment
  communities cut across SuperClasses: the "Organometallic" community is only 23% pure.
  Alternative parents are descriptors by design (Table S1), so co-assignment reflects shared
  substructures, not shared lineage. The low ARI against 641 Classes (0.09) is partly a
  resolution mismatch: 8 communities cannot match 641 classes.
- **Significance is not effect size.** At N = 7.3e7, BiCM validates 2,669 edges with
  obs/exp < 1.1, 48,402 with obs/exp < 2, and 111,482 (half) with obs/exp < 5
  (`R:C….bicm.validated_by_effect_size`). Some are hub–hub, for example `Organic oxides`–
  `Hydrocarbon derivatives` (53.53M observed vs 53.33M expected) at roughly 30 sd.
  - Only 34% of raw pairs have obs/exp > 1 under the hypergeometric null. The Bonferroni
    SVN keeps 14.5% of raw pairs, with a validated median ratio of 5.2 against 0.54 for all
    raw pairs.
  - This is why the primary backbone adds the lift ≥ 2 floor. `ui.edges` holds its top
    6,000 edges by −log₁₀ p.
- **Disparity is a different object.** It keeps the locally dominant edges of each node:
  hub-and-spoke structure, whose clustering (0.81) is mostly reproduced by its degree
  sequence (null 0.75). It shares only ~20–25% of edges with the global nulls. As a local
  filter it does not control for compound-degree heterogeneity, and its taxonomy agreement
  is the lowest of the four.

## D. Class × functional-group bipartite (`16_functional_groups.py`, `R:D…`)

**The matrix.** D has 3,598 direct-parent classes × 396 substituents and 221,536 nonzeros
(fill 0.155).

**Binarisation.** Revealed comparative advantage (Balassa 1965; Hidalgo & Hausmann 2009):

RCA_ds = (D_ds / Σ_s D_ds) / (Σ_d D_ds / Σ D) ≥ 1

This gives 97,819 ones (fill 0.069). RCA is the primary choice because it removes
class-size effects; raw presence (D > 0) is reported alongside.

**Nestedness.** NODF (Almeida-Neto et al. 2008), in its order-free form: for each pair with
k_i > k_j, the share of j's ones that i also has; pairs with equal degree score 0. The null
is curveball, fixed-fixed (Strona et al. 2014): 30 replicates of 5 × rows trades each.

| matrix | NODF | rows / cols | null mean ± sd | z |
|---|---:|---:|---:|---:|
| RCA ≥ 1 | 27.08 | 27.18 / 18.76 | 27.57 ± 0.09 | −5.4 |
| presence | 65.40 | 65.63 / 46.05 | 77.01 ± 0.09 | −136 |

**Negative result: not nested beyond the degree sequence.** The presence matrix is
significantly *less* nested than degree-preserving randomisations, which points to a
compartmentalised structure: groups are specific to chemistries. Its raw NODF of 65 is
produced almost entirely by the degree heterogeneity of class sizes.

**Fitness–Complexity** (Tacchella et al. 2012, Sci Rep 2:723), on the RCA matrix:

- Updates: F_c ← Σ_p M_cp Q_p and Q_p ← 1/Σ_c M_cp F_c⁻¹, each normalised to mean 1.
- 200 iterations. The ranking is stable: Spearman 0.999999 against 400 iterations.
- Fitness correlates with diversification (Spearman 0.81) and with compound count (0.69).
  Complexity anti-correlates with ubiquity (−0.89).
- The top 1% of classes hold 95% of total fitness.

**Artefact.** The top of the ranking is driven by groups held by a single class:

- 37 groups have RCA ubiquity 1, and all 20 top-complexity groups are among them.
- Some are self-incidences: `Other non-metal oxides` is a substituent of itself. These are
  3.16M compounds over 330 classes, 0.47% of the mass.
- The top classes (`Other non-metal oxides`, `Halogen oxides`, `Benzyloxycarbonyls`,
  `Metalloid oxides`, `Fluorobenzenes` …) are therefore holders of exclusive labels, not
  "chemically capable" classes.

Sensitivity variants:

| variant | Spearman vs primary | new top class |
|---|---:|---|
| self-incidences removed | 0.945 | `Benzene and substituted derivatives` |
| also groups with ubiquity < 5 removed | 0.944 | `Oligopeptides`, `Very long-chain fatty acids` … |

The bulk ranking is robust, but the head is not. Read F-C here as a diversity ranking.

**Taxonomy vs functional-group profiles** (Mantel test):

- Sample: 1,000 random classes with ≥ 100 compounds (from 2,546 eligible), giving 499,500
  pairs.
- Similarity: cosine of the RCA rows. Distance: tree path length, depth_a + depth_b −
  2·depth(LCA). 999 permutations, one-sided for a negative association.

| similarity | Pearson r | Spearman | p | null |
|---|---:|---:|---:|---|
| RCA cosine | **−0.203** | −0.194 | **0.001** | 0.000 ± 0.026 |
| raw share profiles (sensitivity) | −0.308 | | 0.001 | |

Further, for the RCA cosine:

- Spearman correlation with LCA depth: 0.36.
- Mean cosine is 0.285 within a SuperClass and 0.210 between SuperClasses.
- By LCA depth 0 → 6, mean cosine rises: 0.05, 0.23, 0.28, 0.32, 0.37, 0.39, 0.39.

**Genuine but weak.** The skeleton taxonomy does order functional-group chemistry, but
explains about 4% of the variance in profile similarity (9% for raw profiles). The release
also strips substituents already implied by the classification, which pushes the signal
down by construction.

**The full cosine projection** (`R:D….cosine_projection`):

- 79% of the 6.47M class pairs have cosine > 0.
- At cosine ≥ 0.5 there are 152,589 edges: density 0.024, average clustering 0.45.

**Process-safety reading** (`R:D….hazard_groups`, `ui.fg.hazard_groups`).

Functional groups, not skeletons, carry thermal instability and reactivity hazards. For the
energetic and reactive groups in the substituent vocabulary:

| group | compounds | classes with RCA ≥ 1 | complexity rank /396 | top SuperClass (share) |
|---|---:|---:|---:|---|
| Organic nitro compounds | 2,764,697 | 348 | 314 | Benzenoids (0.53) |
| Nitroaromatic compounds | 2,773,362 | 259 | 292 | Benzenoids (0.51) |
| Azo compounds | 246,498 | 784 | 338 | Benzenoids (0.34) |
| Diazo compounds | 5,568 | 319 | 265 | Benzenoids (0.43) |
| Organic diazonium salts | 8,167 | 261 | 228 | Benzenoids (0.65) |
| Organic N-nitroso compounds | 10,812 | 401 | 266 | Organoheterocyclics (0.29) |
| Organic perchlorate salts | 9,399 | 297 | 256 | Organoheterocyclics (0.49) |
| Organometallic compounds | 64,676 | 708 | 340 | Benzenoids (0.48) |
| Fluorodinitromethyl compounds | 368 | 13 | 104 | Organic 1,3-dipolar (0.83) |

Organic nitro compounds are present in 1,960 direct-parent classes, and their single largest
class, `Nitrobenzenes`, holds only 16%.

**Engineering conclusion.** Hazard-bearing groups are over-represented in hundreds of
skeleton classes each, and most sit in the low-complexity, high-ubiquity part of the
ranking. Across all 396 groups, the median normalised entropy over 28 SuperClasses is 0.35,
against 0.48 for the compound mass itself, and 55% of groups have a majority in one
SuperClass. So groups *are* somewhat concentrated by the taxonomy, which is consistent with
the weak Mantel signal. But the hazard groups above have top-SuperClass shares of only
0.29–0.65.

A class-level screen ("avoid class X") would therefore miss most compounds carrying a given
reactive group. Safety screening must query substituents directly. The two bipartite
networks (B and S) are complementary, and neither projection replaces the other.

---

## Limitations

1. **Alternative parents are descriptors by design** (Table S1). C measures co-description,
   not competition or similarity of compounds. The hub labels (`Hydrocarbon derivatives`,
   `Organic oxides`) are taxonomy plumbing, and they dominate every unnormalised measure.
2. **The hypergeometric null ignores compound-degree heterogeneity.** A compound with 40
   labels co-assigns 780 pairs, one with 5 labels only 10. The SVN treats all compounds
   alike, which is why BiCM, fixing both degree sequences, is the primary null. In practice
   the two agree closely here (Jaccard 0.81–0.86), and BiCM is stricter than the SVN at
   equal FDR: 222,927 vs 260,078 edges.
3. **Poisson approximation.** The Poisson upper tail is heavier than the Poisson-binomial
   tail, so it is conservative. For hub pairs with q_r near 1 it overstates the variance
   considerably. The refined-normal sensitivity check validates 5% more edges (Jaccard
   0.95).
4. **Huge N.** Statistical significance is cheap at 7.3e7 rows. Validated ≠ strong (§C),
   hence the lift ≥ 2 floor on the primary backbone. The threshold 2 is a convention
   (`extension_research.md` §2.6), not a derived optimum. Plain BiCM remains in every table
   for comparison.
5. **The 200k sample** (`paths.find('sample_200k.tsv')`) is used only for the compound-mode
   density and transitivity. Its representativeness was checked (`R:A_B….representativeness`):
   - It is **not the head of the file**: its CIDs run from 372 to 177.8M and are unsorted.
   - Its direct-parent distribution has total variation distance 0.0260 from
     `path_counts.tsv`. Multinomial draws of the same size from the full distribution give
     0.0262 ± 0.0007 (p = 0.60).
   - Mean \|M\| is 13.603 in the sample vs 13.618 in full. Per-label shares correlate at
     r = 0.99999.

   It is consistent with a uniform random sample.
6. **Louvain** is stochastic (seed 42). It is run on the unweighted backbones, and also
   weighted by ln(lift) on the primary. A handful of large communities is typical at
   resolution 1 on graphs this dense; the NMI against 641 Classes is limited by that
   resolution mismatch.
   - M_B is evaluated for these partitions, not optimised directly.
   - The M_B shuffled null preserves module sizes only. It is a floor, not a test of the
     tree against the best partition.
7. **Rewiring** used one accepted swap per edge per replicate (10 replicates). With more
   mixing the null clustering would fall slightly further, so the reported gap is a lower
   bound.
8. **Dictionary taxonomy.** SuperClass and Class mapping uses `data/chemont_dictionary.tsv`,
   not the OBO. The two differ on 4 parent edges (`how_this_was_built.md` §B.4), which is
   immaterial at these aggregate levels.
9. **Label-set equality** uses a 64-bit hash (`09_scan_bipartite.py`). The collision
   probability over 2.4e7 sets is about 1.6e-5.
10. **Fitness–Complexity** is dominated by exclusive labels (§D), and its head should not be
    read as chemistry.

## How to reproduce

```bash
# inputs: results/bipartite/*.npz + bipartite_scan.json (from scripts/09_scan_bipartite.py,
#         one 4-minute pass over the 2 GB file), results/path_counts.tsv, results/dp_alt_pairs.tsv,
#         data/chemont_dictionary.tsv, sample_200k.tsv (located by scripts/paths.py / $ATLAS_DATA)
PY=.venv/bin/python
$PY scripts/test_bipartite.py              # 16 toy tests, < 1 s
$PY scripts/14_bipartite_projections.py    # A + B, ~10 s
$PY scripts/15_bipartite_backbones.py      # C, ~60 s (BiCM fit 7 s; rewiring nulls ~40 s; M_B nulls 2 s)
$PY scripts/16_functional_groups.py        # D, ~50 s (curveball nulls 20 s, Mantel 25 s)
```

Each script merges its section into `results/bipartite/bipartite_results.json` and rebuilds
`results/bipartite/ui.json`, so they can run in any order. The other outputs are:

- `backbone_edges_{svn_bonferroni,svn_fdr,bicm,bicm_lift2,disparity}.tsv`, with columns
  u, v, weight (C_ab, or W_ab for disparity) and score (−log₁₀ p, or −log₁₀ α for
  disparity);
- `communities_<method>.tsv`, plus `communities_bicm_lift2_weighted.tsv`.

All randomness is seeded (`numpy.random.default_rng(42)`, Louvain seed 42).

## References

- Almeida-Neto M. et al. (2008) A consistent metric for nestedness analysis. *Oikos* 117:1227.
- Balassa B. (1965) Trade liberalisation and revealed comparative advantage. *Manchester School* 33:99.
- Clauset A., Shalizi C.R., Newman M.E.J. (2009) Power-law distributions in empirical data. *SIAM Rev* 51:661.
- Guimerà R., Sales-Pardo M., Amaral L.A.N. (2007) Module identification in bipartite and directed networks. *Phys Rev E* 76:036102.
- Hidalgo C.A., Hausmann R. (2009) The building blocks of economic complexity. *PNAS* 106:10570.
- Newman M.E.J. (2001) Scientific collaboration networks II. *Phys Rev E* 64:016132.
- Newman M.E.J. *Networks: An Introduction*, Oxford UP, ch. 6 (bipartite networks, one-mode projections).
- Saracco F. et al. (2015) Randomizing bipartite networks: the case of the World Trade Web. *Sci Rep* 5:10595.
- Saracco F. et al. (2017) Inferring monopartite projections of bipartite networks. *New J Phys* 19:053022.
- Serrano M.Á., Boguñá M., Vespignani A. (2009) Extracting the multiscale backbone of complex weighted networks. *PNAS* 106:6483.
- Strona G. et al. (2014) A fast and unbiased procedure to randomize ecological binary matrices with fixed row and column totals. *Nat Commun* 5:4114.
- Tacchella A. et al. (2012) A new metrics for countries' fitness and products' complexity. *Sci Rep* 2:723.
- Tumminello M. et al. (2011) Statistically validated networks in bipartite complex systems. *PLoS ONE* 6:e17994.
- Volkova A.Y. (1996) A refinement of the central limit theorem for sums of independent random indicators. *Theory Probab Appl* 40:791.
