# Extension research — the bipartite view of ClassyFire/ChemOnt and its one-mode projections

CLL798D · Group 2 · literature review for the professor's suggestion: "explore the bipartite
nature and the one-mode projection". Compiled 12 Nov 2026. Every reference in §8 was checked by
web search against the publisher, PubMed/PMC, arXiv or Crossref. Items not checked are marked
**[unverified]**. Section numbers for Newman's textbook refer to the **2nd edition** (OUP 2018),
whose table of contents was checked directly.

---

## 0. What we already hold, and what it is sufficient for

From `scripts/09_scan_bipartite.py` → `results/bipartite/` (one 231 s pass over 73,105,281 rows;
`bipartite_scan.json` currently says `"complete": false`, because a rerun is fixing the degree
histogram that was cut off at 64. Every figure below comes from that rerun-in-progress file):

| object | symbol | size | note |
|---|---|---|---|
| label incidence | `B` (compound × category), `M(r) = {DP} ∪ AP(r)` | 995,524,038 ones | not stored per compound |
| category degree | `k_c = Σ_r B_rc` | 3,851 non-zero of 4,832 ids | exact |
| compound-degree histogram | `n_d = #{r : d_r = d}` | d = 1…92, mean 13.62 | exact; **~92 distinct values** |
| class projection | `C = BᵀB` | 1,393,402 non-zero off-diagonal pairs (density ≈ 0.19 over 3,851 nodes) | exact |
| Newman-weighted projection | `W = Bᵀ diag(1/(d_r−1)) B` | same support | exact |
| DP × substituent | `D = PᵀS` | 221,536 non-zero cells, 396 substituents | exact |
| substituent projection | `F = SᵀS` | 28,649 pairs | exact |
| label-set signatures | hash of sorted `M(r)` | 24,264,828 distinct; 14,745,280 singletons; largest 120,678 | **counts only, not the table** |
| per-compound rows | — | 200,000-row sample only | — |

**The two facts that set feasibility for everything below.**

1. *A maximum-entropy null on this bipartite graph needs only the two degree sequences.*
   In the BiCM (§2.3), every compound with the same degree gets the same Lagrange multiplier,
   so the 73.1M-node side collapses to **~92 degree classes**. The whole null model therefore
   has ~92 + 3,851 unknowns and is solvable on a laptop in seconds.
2. *A rescan is cheap.* Re-running the pass to accumulate an extra bounded statistic (for
   example `n_{d,c}`, the number of compounds of degree d carrying category c, a 92 × 4,832 matrix;
   or `Σ_r B_rc B_rc'/d_r`) costs about 4 minutes. Anything that is "not stored" but is bounded by
   the category count costs a rescan, not a redesign.

---

## 1. One-mode projection: definition and known pitfalls

**Definition** (Newman 2018, §6.6.1, p.115ff). For an incidence matrix `B` (groups × items in
Newman's orientation; here compounds × categories), the weighted projection onto categories is
`P = BᵀB` with the diagonal zeroed. `P_cc'` is the number of shared compounds; `P_cc = k_c`.
Newman states the loss directly: the projection "discards a lot of the information present in
the original bipartite network", and even the weighted version cannot recover which pairs were
co-members of *the same* group. The underlying object is a hypergraph (§6.5, p.114); affiliation
networks are introduced in §4.5 (p.60); and §7.2.1 notes that every group becomes a clique.

**Pitfalls, with what each does to our data.**

| pitfall | mechanism | source | our case |
|---|---|---|---|
| clique artefact | each compound with `d_r` labels injects a `K_{d_r}` | Newman §6.6.1, §7.2.1 | mean `d_r` = 13.6, max ≈ 92 → every compound adds up to 4,186 class pairs |
| inflated clustering | triangles come free from shared groups, not from triadic closure | Newman, Strogatz & Watts 2001, Eq. 90; Opsahl 2013 | `C` will be highly clustered *by construction*; any H1b-style clustering claim (`phase2_plan.md` §1) must beat a **bipartite** null, not a one-mode configuration model |
| density inflation by hubs | a few high-degree "groups" connect almost everything | Guillaume & Latapy 2004; Latapy et al. 2008 | compounds with d ≥ 40 (~0.2M) touch most frequent classes; generic classes (e.g. "Organic oxygen compounds"-type nodes) act as hubs in `C` |
| weight ≠ tie strength | raw count `C_cc'` scales like `k_c k_c'` | Newman 2001 PRE; Zhou et al. 2007 | `phase2_plan.md` already measured r = 0.788 between overlay weight and size |
| information loss | different bipartite graphs → identical projections | Guillaume & Latapy 2004; Latapy et al. 2008 | 24.3M distinct label sets collapse into 1.39M pair counts |

**Closed forms for the compound-mode projection (already planned; exact from class sizes).**
At a single hierarchy level the compounds are partitioned, so the projection is a disjoint union
of cliques. With class sizes `n_c` and `N = Σ n_c`:

```
edges          E   = Σ_c n_c (n_c − 1) / 2
degree         k_r = n_{c(r)} − 1                    (every compound in class c)
degree dist    P(k) = Σ_{c : n_c = k+1} n_c / N      (size-biased class-size distribution)
mean degree    <k> = Σ_c n_c² / N − 1
clustering     C_local = 1 for every r with n_{c(r)} ≥ 3;   components = #classes
```

For the full `M(r)` label set the cliques overlap, and only bounds are available without
per-compound data: `max_{c∈M(r)} n_c − 1 ≤ k_r ≤ Σ_{c∈M(r)} n_c − |M(r)|`. This is the textbook
demonstration that the clustering and the heavy tail of a projection can be pure artefacts.

**Expected clustering of the class projection under a bipartite configuration model**
(Newman, Strogatz & Watts 2001, Eqs. 90–91, verified in the arXiv text; textbook §12.11.2).
Map "actors" to categories and "movies" to compounds:

```
1/C − 1 = (μ2 − μ1)(ν2 − ν1)² / [ μ1 ν1 (2ν1 − 3ν2 + ν3) ]
μ_n = n-th moment of category degree k_c        (from k_c)
ν_n = n-th moment of compound degree d_r        (from n_d)
```

This is computable now. **Caveat:** the derivation is a sparse-limit result: it counts every
triangle and ignores repeated co-membership. Our projection has density ≈ 0.19, so treat it as a
qualitative baseline, and use the BiCM expectations of §2 for quantitative z-scores.

**Weighting schemes.** Each one changes what an edge means:

```
raw count          C_cc'    = Σ_r B_rc B_rc'
Newman (2001)      W_cc'    = Σ_r B_rc B_rc' / (d_r − 1)               (d_r ≥ 2)
resource alloc.    WRA_cc'  = (1/k_c') Σ_r B_rc B_rc' / d_r            (Zhou et al. 2007; asymmetric)
```

Newman's identity gives a free correctness check on our `W`: `Σ_{c'≠c} W_cc' = #{r ∈ c : d_r ≥ 2}`,
the number of compounds in c that carry at least one other label. `WRA` needs `Σ_r B_rc B_rc'/d_r`.
This is not derivable from `W`, so it costs a rescan (~4 min) or a 200k-sample estimate.
`WRA` is the natural *directed* weight for the plan's "where does uncertainty flow" PageRank
(`phase2_plan.md` §5).

**Two-mode clustering** (Opsahl 2013). Clustering is redefined on the bipartite graph itself as
`C* = (#closed 4-paths) / (#4-paths)`, where a 4-path is class–compound–class–compound–class and
it is "closed" if it lies on a 6-cycle. It needs per-compound adjacency, so it is feasible only on
the 200k sample. Latapy et al. (2008) give the matching two-mode density, redundancy and
degree-correlation statistics, and benchmark them against random bipartite graphs with the same
degree sequences.

---

## 2. Statistically validated projections (null models)

Goal: keep edge (c, c') only if c and c' share **more compounds than the degrees alone predict**.
All methods below differ only in the null model.

### 2.1 Hypergeometric SVN (Tumminello et al. 2011)

The null is that c's compounds are a uniform random subset of all `N_T` compounds:

```
p(c,c') = P(X ≥ C_cc'),   X ~ Hypergeom(N_T, k_c, k_c')
        = Σ_{X ≥ C_cc'} binom(k_c, X) binom(N_T − k_c, k_c' − X) / binom(N_T, k_c')
scipy:  hypergeom.sf(C_cc' − 1, N_T, k_c, k_c')
```

Tumminello et al. control the compound-side heterogeneity by **stratifying**: they test separately
within each subsystem of B-nodes with the same degree, then apply Bonferroni (`0.01/#tests`) or FDR
over all tests. *For us:* the unstratified version is feasible now (needs `k_c`, `C`, `N_T`). The
stratified version needs `n_{d,c}` and the per-degree pair counts `C^(d)_cc'`; the first is a cheap
rescan, the second is a rescan of up to 92 sparse matrices of the same support as `C`. The
unstratified version **over-validates** pairs that are concentrated in high-degree compounds.
That is the failure BiCM fixes.

### 2.2 FDR (Benjamini & Hochberg 1995)

Sort the M p-values `p_(1) ≤ … ≤ p_(M)`, find the largest `î` with `p_(î) ≤ î·t/M`, and reject all
`p ≤ p_(î)`. Saracco et al. (2017) use `t = 0.01`. Here `M = C(3851, 2) ≈ 7.41 × 10⁶` (pairs with
`C_cc' = 0` can be skipped for over-representation tests, but they still count in M).

### 2.3 BiCM (Saracco et al. 2015) — the recommended null

The maximum-entropy ensemble of bipartite graphs whose **expected** degrees on both layers equal
the observed ones:

```
p_rc = x_r y_c / (1 + x_r y_c)
constraints:  Σ_c p_rc = d_r  for every compound r;    Σ_r p_rc = k_c  for every category c
```

**Grouped form for our data** (exact, because the solution depends on r only through `d_r`):

```
d   = Σ_c  x_d y_c / (1 + x_d y_c)            for each d with n_d > 0   (~92 equations)
k_c = Σ_d  n_d x_d y_c / (1 + x_d y_c)        for each c with k_c > 0   (3,851 equations)
fixed point:  x_d ← d / Σ_c [ y_c / (1 + x_d y_c) ] ;   y_c ← k_c / Σ_d [ n_d x_d / (1 + x_d y_c) ]
check:        max relative residual on both constraint sets < 1e−8
```

Vallarano et al. (2021) compare Newton, quasi-Newton and fixed-point solvers for these models and
ship the NEMtropy package, which handles binary bipartite models. Whether NEMtropy accepts a
*grouped* degree sequence directly is **[unverified]**; the 20-line numpy loop above avoids the
question. One caveat: an exact degree `d_r = |C|` (a compound carrying every category) would be
degenerate, but max d ≈ 92 ≪ 3,851, so it cannot arise here.

### 2.4 V-motif projection test (Saracco et al. 2017)

The similarity of c, c' is the number of V-motifs (shared compounds), `V_cc' = Σ_r B_rc B_rc' = C_cc'`.
Under BiCM the per-compound events are independent Bernoulli trials with probability
`p_rc p_rc'`, so `V_cc'` is **Poisson–binomial**. Grouped by compound degree, it is a sum of
binomials:

```
V_cc' ~ Σ_d Binomial(n_d, q_d),     q_d = p(d,c) · p(d,c')
mean  μ  = Σ_d n_d q_d
var   σ² = Σ_d n_d q_d (1 − q_d)
p-value  = P(V ≥ C_cc')            then BH-FDR over all pairs (t = 0.01)
```

Saracco et al. (Appendix B) discuss approximations: Poisson when μ is small; Gaussian with
continuity correction when μ is large; a skewness-corrected Gaussian in between. They warn that
the approximations can fail on real data. Practical rule: exact tail by convolution over the
~92 binomials for pairs with μ < ~50, and a skewness-corrected normal above that.
`E_BiCM[C_cc'] = μ` is also the **correct expected weight** to subtract in projection modularity
and to report next to every overlay edge.

### 2.5 Other backbones worth one comparison each

* **SDSM** (Neal 2014). Same idea as BiCM (degree sequences fixed *in expectation*). It
  originally used Monte-Carlo-estimated cell probabilities; the R `backbone` package (Neal 2022)
  now computes them with BiCM. Its **signed** variant also flags pairs that co-occur
  *significantly less* than expected. The **SDSM-EC** extension (Neal & Neal 2024) adds
  *prohibited edges*. For ChemOnt, any pair the rules forbid (R5: no ancestor/descendant pair
  inside `M(r)`) is a structural zero and must be excluded from the null and the test.
* **Disparity filter** (Serrano, Boguñá & Vespignani 2009). A local null on the *projection*,
  not on the bipartite graph: `α_cc' = (1 − w_cc'/s_c)^(k_c − 1)`, where `k_c` is the projected
  degree. Keep the edge if `α < α*` at either endpoint. It is cheap and runs on `C` or `W`
  directly. It ignores the bipartite origin, so use it as a contrast, not as the primary backbone.
* **Gualdi et al. (2016)** apply a BiCM-type validated projection to portfolio overlaps, a
  template for "validated co-holding" analogous to "validated co-assignment". **Cimini et al.
  (2019)** is the review that ties all of these maximum-entropy null models together; cite it for
  the framework.

### 2.6 Two caveats specific to a 73-million-row bipartite graph

1. **Significance is cheap at this N.** Two classes with `k = 10⁶` each have a hypergeometric
   expected overlap of `10¹²/7.31×10⁷ ≈ 13,700` with σ ≈ 117. A +3 % excess (≈ 14,100) already gives
   z ≈ 3.4. FDR will validate many substantively trivial edges, so **always pair the p-value with
   an effect size**, for example `lift = C_cc'/μ` or `log2 lift`, and report backbones at
   lift ≥ 2 as well as at FDR alone.
2. **Validation does not separate descriptive from competing labels.** Most validated
   edges will be *definitional* (a compound that is an aryl thioether is also a thioether-type
   compound). BiCM removes size effects, not chemistry. `phase2_plan.md` §3 (descriptive vs
   competing) remains necessary; the validated projection is the better *input* to it.

---

## 3. Bipartite community detection

* **Barber (2007)** — bipartite modularity with a bipartite null
  `Q = (1/m) Σ_ij (B_ij − k_i d_j / m) δ(g_i, h_j)`, plus the **BRIM** algorithm. BRIM alternately
  assigns each side given the other, and Q never decreases. It needs the bipartite matrix,
  so it is infeasible on compounds × categories at 73M rows. It is **feasible on DP × substituent
  `D`** (3,851 × 396) and on the 200k sample. Python: pyCONDOR.
* **Guimerà, Sales-Pardo & Amaral (2007)** — bipartite modularity *for a partition of one layer*,
  defined **entirely from the projection plus degrees**:

  ```
  M_B(P) = Σ_s [ Σ_{i≠j∈s} c_ij / Σ_a m_a(m_a − 1)  −  Σ_{i≠j∈s} t_i t_j / (Σ_a m_a)² ]
  ours:   c_ij = C_cc'   t_i = k_c   m_a = d_r   Σ_a m_a(m_a−1) = Σ_d n_d d(d−1)   Σ_a m_a = Σ_d n_d d
  ```

  **Every term is already stored.** This is the right modularity for H1c. Plain Newman–Girvan Q on
  `C` uses a one-mode null `s_c s_c'/2L`, which double-counts the clique structure. M_B can be
  optimised by feeding Louvain/Leiden the modularity matrix
  `B_cc' = C_cc'/Σd(d−1) − k_c k_c'/(Σd)²`.
* **LPAwb+ / DIRTLPAwb+** (Beckett 2016) — weighted bipartite modularity via label propagation.
  It is fast, so it suits the weighted `D` matrix. Code: Julia, MATLAB/Octave, R.
* **biLouvain** (Pesántez-Cabrera & Kalyanaraman 2016) — Louvain for the Murata+ bipartite
  modularity. It is one to three orders of magnitude faster than the alternatives the authors
  tested; C++ code.

**Projection vs bipartite communities: what to compare.** Three partitions of the 3,851 categories:

1. Louvain on `C` with Newman–Girvan Q (the naive projection);
2. Louvain on the M_B modularity matrix (bipartite null, still projection-only data);
3. Louvain on the BiCM-validated backbone (§2.4).

Compare each against the ChemOnt SuperClass partition (NMI/ARI) and against each other. The
direct precedent is **Glass & Girvan (2015)**. They projected a gene–GO-term bipartite network
onto terms, ran community detection, and showed that the term communities differ strongly from
the GO DAG's own branches while carrying biological signal. That is exactly the H1c question
(`phase2_plan.md` §1) on a different ontology.

---

## 4. Nestedness and economic complexity on the class × functional-group matrix

**Input.** The weighted matrix `D` (direct parent × substituent). Binarise it with the revealed
comparative advantage used by Hidalgo & Hausmann (2009):

```
RCA_cs = (D_cs / Σ_s D_cs) / (Σ_c D_cs / Σ_cs D_cs);     M_cs = 1  iff  RCA_cs ≥ 1
```

This means "class c over-uses substituent s relative to the whole release". Raw `D_cs > 0` would
mostly measure class size.

**Nestedness — NODF** (Almeida-Neto et al. 2008). Order the rows and columns by marginal total.
For every row pair (i above j), `DF_ij = 100` if `MT_j < MT_i`, else 0, and `PO_ij` is the percentage
of j's 1s that are also in i. Then `N_paired = PO` if `DF = 100`, else 0. Do the same for column
pairs, and compute

```
NODF = Σ N_paired / [ n(n−1)/2 + m(m−1)/2 ]          (0–100; row and column parts separately available)
```

*Interpretation:* high NODF means specialist classes use **subsets** of the substituents used by
generalist classes, and rare substituents occur only in classes that also carry the common ones.
*Trap:* ChemOnt builds some nestedness in by design, because a child class inherits its parent's
defining groups. NODF must be tested against a null that preserves the margins (BiCM on `M`, or a
fixed-fixed swap), and the hierarchy-only nestedness should be reported separately. Mariani et
al. (2019) review both the measures and the null-model pitfalls.

**Method of Reflections** (Hidalgo & Hausmann 2009):

```
k_c,0 = Σ_s M_cs  (diversity)          k_s,0 = Σ_c M_cs  (ubiquity)
k_c,N = (1/k_c,0) Σ_s M_cs k_s,N−1     k_s,N = (1/k_s,0) Σ_c M_cs k_c,N−1
```

**Fitness–Complexity** (Tacchella et al. 2012). This is a non-linear map; normalise by the mean at
every step, starting from `F⁰ = Q⁰ = 1`:

```
F̃_c(n) = Σ_s M_cs Q_s(n−1)              Q̃_s(n) = 1 / Σ_c M_cs / F_c(n−1)
F_c(n) = F̃_c(n) / <F̃(n)>               Q_s(n) = Q̃_s(n) / <Q̃(n)>
```

(The structure of the two steps is verified in the paper's text. The algebra above is the
standard published form; check it against Eq. 1 of the paper before coding.)

*Interpretation for chemistry.* A class's "fitness" is the **functional-group versatility** of
its members. A substituent's "complexity" is high when it appears only in versatile classes, and
it is driven down by a single low-diversity class that also uses it. A complexity ranking
separates "generic" groups (hydroxyl- or carbonyl-type substituents that even simple classes
carry) from "sophisticated" ones (groups that occur only in polyfunctional natural-product-like
classes). This hypothesis is testable against class depth in the tree and against the
Phenylpropanoid/terpenoid/peptide blocks. MoR's fixed point is an eigenvector of the class–class
transition matrix `M̃ = D_c⁻¹ M D_s⁻¹ Mᵀ`, which links it to eigenvector centrality (Newman §7.1).

---

## 5. Hypergraph view: what pairwise projection loses

Each compound is a hyperedge over its label set `M(r)`. Battiston et al. (2020) review hypergraph
and simplicial representations and show that pairwise projection loses group-level
information. In particular, a 3-way co-occurrence cannot be reconstructed from its three pairwise
counts. Newman (§6.5) makes the same point: the bipartite graph *is* the hypergraph; the projection
is not.

*For us:*

* the hyperedge-size distribution is exactly `n_d` (already stored);
* hyperedge multiplicities are the label-signature counts: 24.3M distinct label sets,
  14.7M seen once, and the largest seen 120,678 times. Only these summary counts are stored, so the
  exact hypergraph would need the signature table written out (bounded by 24.3M rows, ≈ GB-scale;
  a larger rescan, not a 4-min one);
* a feasible higher-order test on the 200k sample: for category triples (a, b, c), compare the
  observed triple count with the count predicted from pairwise counts alone (for example, the
  maximum-entropy estimate given the three pair counts). Over-represented triples are genuine
  three-way "chemotypes" that the projection splits into three edges. *Effort/benefit is low*
  relative to §2–§3; treat it as optional.

---

## 6. Domain precedents

| study | bipartite object | projection used | relevance |
|---|---|---|---|
| Goh et al. 2007, PNAS (human diseasome) | disorders × disease genes (OMIM; 1,284 × 1,777) | both: disease network and gene network, edge weight = shared genes | canonical "bipartite → two projections" figure; raw-count weights, no null |
| Yıldırım et al. 2007, Nat Biotechnol (drug–target) | FDA drugs × protein targets | drug–drug, target–target | same template, pharmacology |
| Hidalgo et al. 2007, Science (product space) | countries × products (RCA ≥ 1) | product–product proximity (min of conditional co-export probabilities) | asymmetric conditional weighting ≈ plan's `w_cond`; core–periphery result |
| Hidalgo & Hausmann 2009; Tacchella et al. 2012 | countries × products | none; bipartite iterative ranking | §4 |
| Glass & Girvan 2015, PLoS Comput Biol | genes × GO terms | term–term, then communities vs GO branches | **closest precedent to H1c** |
| Wood et al. 2020, Open Biol ("Term Matrix") | genes × GO terms | term co-annotation matrix | co-annotation of mutually exclusive terms flags **curation errors**; 289 manual and >52,700 automatic annotations corrected (abstract figures). This is the GO analogue of our R5 violations and of a *signed* backbone's negative edges |
| Tsekenis et al. 2024, Sci Rep | inorganic compounds × elements | compound–compound, element–element; fitness (max-ent) null | maximum-entropy null on a chemistry bipartite graph; oxygen as hub |
| Fialkowski et al. 2005, Angew Chem | compounds → reactions (Beilstein) | directed one-mode | scale-free claims in chemistry networks (pre-Clauset–Shalizi–Newman; read critically) |
| Ertl 2017, J Cheminform | molecules × functional groups (3,080 groups in ChEMBL) | none | defines a functional-group vocabulary; context for the 396 ClassyFire substituents |
| Ernst et al. 2019, Metabolites (MolNetEnhancer) | ClassyFire labels mapped onto MS/MS molecular networks | none | ClassyFire used as **node annotation**, not analysed as a network |

**ClassyFire/ChemOnt as a network: what the search found.** Targeted searches (ClassyFire or
ChemOnt with "co-occurrence", "network", "projection", "alternative parents") found **no
network-science analysis of ClassyFire label co-occurrence or of the alternative-parent
structure**. The uses found are annotation layers on spectral networks (MolNetEnhancer; ConCISE,
Metabolites 2022), class matching between genome-mined and metabolome-mined compounds
(NPClassScore, PMC9869629), and R accessors for `alternative_parents` (classyfireR). This is
absence of evidence from a bounded search, not proof of novelty; state it that way in the report.
For ChEBI, graph statistics exist (node/edge counts, relation-type shares), but we found no
co-annotation network study.

---

## 7. Ranked recommendations

Columns: **(a)** what it would show; **(b)** feasibility with stored aggregates
(✔ now · ↻ needs a ~4-min rescan · S = 200k sample only · ✘ infeasible); **(c)** syllabus link
(Newman 2nd ed.); **(d)** effort; **plan** = relation to `phase2_plan.md`.

| # | analysis | (a) what it shows | (b) | (c) syllabus | (d) | plan |
|---|---|---|---|---|---|---|
| 1 | **BiCM-validated class projection** (grouped BiCM §2.3 + Poisson–binomial V-motif test §2.4 + BH-FDR, with a lift ≥ 2 effect-size floor) | which class pairs co-occur beyond what compound and class degrees force; gives `μ_cc'` for every edge | ✔ (`k_c`, `n_d`, `C`) | configuration model §12.1, bipartite CM §12.11.2; significance vs null | 1–2 d | **new**; replaces the heuristic Chung–Lu "size-preserving edge null" of §6.1 for the projection, and is the principled input to §3–4 |
| 2 | **Guimerà M_B modularity** on `C`, compared with NG-modularity Louvain and the backbone partition; NMI/ARI vs SuperClass | whether "communities" survive a bipartite null, i.e. H1c done correctly | ✔ | modularity §14.2; assortative mixing §7.7 | 0.5–1 d | **extends** planned Louvain/NMI (§5, H1c) |
| 3 | **Projection-artefact baseline**: compound-mode clique closed forms; NSW Eq. 90 expected clustering vs observed clustering of `C`; Newman identity check on `W` | how much of the projection's clustering and heavy tail is built in | ✔ | clustering §7.3; cliques §7.2.1; degree distributions §10.3–10.4 | 0.5 d | **planned** (closed forms; the docstring names `scripts/14_bipartite_analysis.py`, which is not yet in `scripts/`) |
| 4 | **Signed backbone** (under-representation tail of #1) | class pairs ChemOnt *avoids* co-assigning; compare with the 93 R5-violating pairs, as in Term Matrix | ✔ (same computation as #1) | null models; structural balance §7.5 (signed edges) | +0.5 d | **new** |
| 5 | **RCA → NODF + Fitness–Complexity** on `D` | whether functional-group repertoires are nested; ranking of substituent "sophistication" and class versatility | ✔ | bipartite networks §6.6; eigenvector-type centrality §7.1; degree distributions | 1 d | **new** |
| 6 | **BRIM / LPAwb+ on `D`** (3,851 × 396, weighted) | co-modules of classes and functional groups (which groups *define* which class blocks) | ✔ | modularity §14.2 | 1 d | **new** |
| 7 | **Weighting comparison**: raw `C` vs Newman `W` vs resource allocation `WRA` (directed) for centrality/PageRank | sensitivity of rankings to the weighting choice; `WRA` is the natural directed flow | ✔ for `C`, `W`; ↻ for `WRA` | weighted networks §6.3; centrality §7.1 | 0.5 d (+ rescan) | **extends** planned PageRank on `w_cond` |
| 8 | **Disparity-filter backbone** of `C`/`W` | local-null backbone; contrast with #1 (what the bipartite null adds) | ✔ | weighted networks; degree heterogeneity | 0.25 d | **new** (contrast only) |
| 9 | **Stratified hypergeometric SVN** (Tumminello) | cross-check of #1 under a different null | ↻ (`n_{d,c}` + per-degree `C^(d)`) | hypergeometric null; similarity §7.6.1 | 0.5 d + rescan | optional |
| 10 | **Two-mode clustering** (Opsahl) and Latapy et al. statistics | clustering without projection artefacts | S | clustering §7.3 | 0.5 d | optional |
| 11 | **Hyperedge triple test** | genuine 3-way chemotypes lost by projection | S | (beyond syllabus: hypergraphs §6.5) | 1–2 d | optional |
| 12 | BRIM/biLouvain on compounds × categories | full bipartite communities | ✘ at 73M; S only | modularity | — | skip |

**Power-law fits** of `k_c`, `n_d` and the projected degrees remain in the plan (`phase2_plan.md`
§5: CSN with likelihood-ratio tests). Note that `n_d` is a *hyperedge-size* distribution: it is
bell-shaped with a tail to d ≈ 92, so it is not a candidate for a power law.

**Suggested minimal package for the report** (≈ 4 working days): #3 → #1 (+#4 free) → #2 → #5.
Together these answer the professor's prompt directly. The projection (#3) is shown to be mostly
artefact. A statistically validated projection (#1) is built from the bipartite graph's own
degree sequences. Its communities are judged with a bipartite modularity (#2). One genuinely
bipartite analysis (#5) uses the second relation, functional groups.

---

## 8. References (all verified unless marked)

1. Newman M. E. J. (2018). *Networks*, 2nd ed. Oxford University Press. doi:10.1093/oso/9780198805090.001.0001. Sections used: 4.5, 6.5, 6.6–6.6.1, 7.1–7.3, 7.6.1, 12.11.2, 14.2 (table of contents checked).
2. Newman M. E. J. (2001). The structure of scientific collaboration networks. *PNAS* 98(2):404–409. doi:10.1073/pnas.98.2.404
3. Newman M. E. J. (2001). Scientific collaboration networks. II. Shortest paths, weighted networks, and centrality. *Phys. Rev. E* 64:016132. doi:10.1103/PhysRevE.64.016132 (weighting `1/(n_k − 1)`, checked in the text).
4. Newman M. E. J., Strogatz S. H., Watts D. J. (2001). Random graphs with arbitrary degree distributions and their applications. *Phys. Rev. E* 64:026118. doi:10.1103/PhysRevE.64.026118 (Eqs. 90–91 checked in arXiv:cond-mat/0007235).
5. Guillaume J.-L., Latapy M. (2004). Bipartite structure of all complex networks. *Inf. Process. Lett.* 90(5):215–221. doi:10.1016/j.ipl.2004.03.007
6. Latapy M., Magnien C., Del Vecchio N. (2008). Basic notions for the analysis of large two-mode networks. *Social Networks* 30(1):31–48. doi:10.1016/j.socnet.2007.04.006
7. Zhou T., Ren J., Medo M., Zhang Y.-C. (2007). Bipartite network projection and personal recommendation. *Phys. Rev. E* 76:046115. doi:10.1103/PhysRevE.76.046115 (resource-allocation formula checked in arXiv:0707.0540).
8. Opsahl T. (2013). Triadic closure in two-mode networks: redefining the global and local clustering coefficients. *Social Networks* 35(2):159–167. doi:10.1016/j.socnet.2011.07.001
9. Tumminello M., Miccichè S., Lillo F., Piilo J., Mantegna R. N. (2011). Statistically validated networks in bipartite complex systems. *PLoS ONE* 6(3):e17994. doi:10.1371/journal.pone.0017994
10. Saracco F., Di Clemente R., Gabrielli A., Squartini T. (2015). Randomizing bipartite networks: the case of the World Trade Web. *Sci. Rep.* 5:10595. doi:10.1038/srep10595
11. Saracco F., Straka M. J., Di Clemente R., Gabrielli A., Caldarelli G., Squartini T. (2017). Inferring monopartite projections of bipartite networks: an entropy-based approach. *New J. Phys.* 19:053022. doi:10.1088/1367-2630/aa6b38 (V-motifs, Poisson–binomial, FDR at t = 0.01 checked in the text).
12. Neal Z. P. (2014). The backbone of bipartite projections: inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors. *Social Networks* 39:84–97. doi:10.1016/j.socnet.2014.06.001
13. Neal Z. P. (2022). backbone: an R package to extract network backbones. *PLoS ONE* 17:e0269137. doi:10.1371/journal.pone.0269137 (as cited in the CRAN `sdsm` reference page; the paper itself was not opened).
14. Neal Z. P., Neal J. W. (2024). Stochastic Degree Sequence Model with Edge Constraints (SDSM-EC) for backbone extraction. In *Complex Networks & Their Applications XII*, SCI vol. 1141, pp. 127–136. doi:10.1007/978-3-031-53468-3_11
15. Serrano M. Á., Boguñá M., Vespignani A. (2009). Extracting the multiscale backbone of complex weighted networks. *PNAS* 106(16):6483–6488. doi:10.1073/pnas.0808904106
16. Gualdi S., Cimini G., Primicerio K., Di Clemente R., Challet D. (2016). Statistically validated network of portfolio overlaps and systemic risk. *Sci. Rep.* 6:39467. doi:10.1038/srep39467
17. Cimini G., Squartini T., Saracco F., Garlaschelli D., Gabrielli A., Caldarelli G. (2019). The statistical physics of real-world networks. *Nat. Rev. Phys.* 1(1):58–71. doi:10.1038/s42254-018-0002-6
18. Vallarano N., Bruno M., Marchese E., Trapani G., Saracco F., Cimini G., Zanon M., Squartini T. (2021). Fast and scalable likelihood maximization for exponential random graph models with local constraints. *Sci. Rep.* 11:15227. doi:10.1038/s41598-021-93830-4
19. Benjamini Y., Hochberg Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* 57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
20. Barber M. J. (2007). Modularity and community detection in bipartite networks. *Phys. Rev. E* 76:066102. doi:10.1103/PhysRevE.76.066102
21. Guimerà R., Sales-Pardo M., Amaral L. A. N. (2007). Module identification in bipartite and directed networks. *Phys. Rev. E* 76:036102. doi:10.1103/PhysRevE.76.036102 (Eq. 5 checked in PMC2262941).
22. Beckett S. J. (2016). Improved community detection in weighted bipartite networks. *R. Soc. Open Sci.* 3:140536. doi:10.1098/rsos.140536
23. Pesántez-Cabrera P., Kalyanaraman A. (2016). Detecting communities in biological bipartite networks. *Proc. ACM-BCB '16*, pp. 98–107. doi:10.1145/2975167.2975177
24. Almeida-Neto M., Guimarães P., Guimarães P. R. Jr, Loyola R. D., Ulrich W. (2008). A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. *Oikos* 117(8):1227–1239. doi:10.1111/j.0030-1299.2008.16644.x
25. Mariani M. S., Ren Z.-M., Bascompte J., Tessone C. J. (2019). Nestedness in complex networks: observation, emergence, and implications. *Phys. Rep.* 813:1–90. doi:10.1016/j.physrep.2019.04.001
26. Hidalgo C. A., Hausmann R. (2009). The building blocks of economic complexity. *PNAS* 106(26):10570–10575. doi:10.1073/pnas.0900943106
27. Tacchella A., Cristelli M., Caldarelli G., Gabrielli A., Pietronero L. (2012). A new metrics for countries' fitness and products' complexity. *Sci. Rep.* 2:723. doi:10.1038/srep00723
28. Hidalgo C. A., Klinger B., Barabási A.-L., Hausmann R. (2007). The product space conditions the development of nations. *Science* 317(5837):482–487. doi:10.1126/science.1144581
29. Battiston F., Cencetti G., Iacopini I., Latora V., Lucas M., Patania A., Young J.-G., Petri G. (2020). Networks beyond pairwise interactions: structure and dynamics. *Phys. Rep.* 874:1–92. doi:10.1016/j.physrep.2020.05.004
30. Goh K.-I., Cusick M. E., Valle D., Childs B., Vidal M., Barabási A.-L. (2007). The human disease network. *PNAS* 104(21):8685–8690. doi:10.1073/pnas.0701361104
31. Yıldırım M. A., Goh K.-I., Cusick M. E., Barabási A.-L., Vidal M. (2007). Drug–target network. *Nat. Biotechnol.* 25:1119–1126. doi:10.1038/nbt1338
32. Glass K., Girvan M. (2015). Finding new order in biological functions from the network structure of gene annotations. *PLoS Comput. Biol.* 11(11):e1004565. doi:10.1371/journal.pcbi.1004565
33. Wood V., Carbon S., Harris M. A., Lock A., Engel S. R., Hill D. P., Van Auken K., Attrill H., Feuermann M., Gaudet P., Lovering R. C., Poux S., Rutherford K. M., Mungall C. J. (2020). Term Matrix: a novel Gene Ontology annotation quality control system based on ontology term co-annotation patterns. *Open Biol.* 10:200149. doi:10.1098/rsob.200149
34. Tsekenis G., Cimini G., Kalafatis M., Giacometti A., Gili T., Caldarelli G. (2024). Network topology mapping of chemical compounds space. *Sci. Rep.* 14:5266. doi:10.1038/s41598-024-54594-9
35. Fialkowski M., Bishop K. J. M., Chubukov V. A., Campbell C. J., Grzybowski B. A. (2005). Architecture and evolution of organic chemistry. *Angew. Chem. Int. Ed.* 44(44):7263–7269. doi:10.1002/anie.200502272
36. Ertl P. (2017). An algorithm to identify functional groups in organic molecules. *J. Cheminform.* 9:36. doi:10.1186/s13321-017-0225-z
37. Ernst M. *et al.* (2019). MolNetEnhancer: enhanced molecular networks by integrating metabolome mining and annotation tools. *Metabolites* 9(7):144. doi:10.3390/metabo9070144 (DOI and title verified; **full author list [unverified]**).
38. Djoumbou Feunang Y. *et al.* (2016). ClassyFire: automated chemical classification with a comprehensive, computable taxonomy. *J. Cheminform.* 8:61. doi:10.1186/s13321-016-0174-y (baseline).

Mentioned but not given a full citation here: ConCISE (*Metabolites* 2022, mdpi.com/2218-1989/12/12/1275) and NPClassScore (PMC9869629). We saw their abstracts only; authors **[unverified]**.
