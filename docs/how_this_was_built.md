# How this project was built

A complete account of the work, in the order it happened — including the things that
went wrong, because several of them changed the results and all of them are material
for the course's mandatory Declaration of Tool Usage reflection.

Nothing here is reconstructed from intent. Every number quoted was measured, and where
a figure was later corrected, both the wrong and the right value are shown.

---

# Part I — Reference: what the data actually is

Read this before the narrative if you need to explain the project to someone. Everything
here is a fact about the inputs and outputs, not about the process.

## A. The domain in one page

**ClassyFire** is a program that reads a molecule's *structure* and assigns it to a
category. It uses no biology, no literature, no bioactivity — structure only. It works by
matching the molecule against a library of >9000 hand-written **SMARTS** patterns
(a substructure query language: `[$([#16]-1-[#6]=[#6]-[#6]=[#7]-1)]` is the paper's own
example for thiazoles), plus **Markush** structures for patterns SMARTS cannot express,
plus ~200 regular expressions over IUPAC names for cases like leukotrienes where no
single backbone works.

**ChemOnt** is the taxonomy it assigns into: a strict tree of chemical categories.

The paper's four algorithmic steps:

| Step | What happens |
|---|---|
| 1. Preprocessing | Input (SMILES / SDF / InChI / IUPAC name / FASTA) becomes a chemical object; physico-chemical properties computed via ChemAxon JChem |
| 2. Feature extraction | Superstructure search over the SMARTS/Markush library; logical rules; IUPAC-name parsing |
| 3. Category assignment + reduction | Each feature maps to a category; for any parent–child pair among the hits, **only the child is kept** (this is rule **R4**) |
| 4. Direct-parent selection | The surviving category with the **largest** structural feature wins, measured in non-hydrogen atoms; ties broken by ring counts, heteroatoms, fused rings |

### The five hierarchy levels

Named after the Linnaean scheme. `depth(root) = 0`, so:

| Depth | Level | Count in ChemOnt 2.1 | Example |
|---|---|---:|---|
| 1 | **Kingdom** | 2 | Organic compounds / Inorganic compounds |
| 2 | **SuperClass** | 31 (26 organic + 5 inorganic) | Lipids and lipid-like molecules |
| 3 | **Class** | 765 | Fatty Acyls |
| 4 | **SubClass** | 1,729 | Fatty acid esters |
| 5–11 | below SubClass | 2,297 | Acyl carnitines |

**Counts here are from the OBO** (765 / 1,729 / 2,297). The dataset's dictionary gives
**766 / 1,729 / 2,296**, because the 4 drifted parent edges (§B.4) move four nodes
between levels:

| Node | depth in OBO | depth in dictionary |
|---|---:|---:|
| `Sulfinylamines` | 4 | 3 |
| `Ureides` | 5 | 4 |
| `1,2-diacyl-3-O-beta-D-galactosyl-sn-glycerols` | 7 | 6 |
| `1-acyl-3-O-beta-D-galactosyl-sn-glycerols` | 7 | 6 |

Always say which artefact a level count came from. The two disagree, and the difference
is exactly these four nodes.

Kingdom is decided purely by molecular formula: organic = contains carbon, with a short
list of exceptions (cyanides, CO, CO₂, CS₂ …) that are called inorganic anyway.

### The four label types — this distinction matters more than anything else

Taken from the paper's Additional file 1, Table S1:

| Label | Definition | Is it a competing candidate? |
|---|---|---|
| **Direct parent** | "the category corresponding to the largest skeleton or most dominant feature". May sit at *any* depth, not necessarily 5 | This is *the* assignment |
| **Alternative parents** | "other categories that describe the compound and do not display a parent-child relationship to each other or to the direct parent" | **No — descriptors by design.** Separating the few that *are* competing is the whole Phase 2 problem |
| **Intermediate nodes** | "descendants of the subclass and ascendants of the direct parent" | No — they are lineage filler |
| **Substituents** | functional groups present, with anything already implied by the classification removed | No |

**Note the paper contradicts itself here.** The Results section says alternative parents
have no *ancestor–descendant* relation to the direct parent; Table S1 says the weaker
*parent-child*. Phase 1c tests both readings, and the violations breach both.

---

## B. The source data, field by field

### B.1 The compound file — `classyfire_dedup_inchikey_smiles.enriched.tsv.zst`

2.0 GB compressed, **24.4 GB** decompressed, **73,105,281** rows, tab-separated, one
header line. Six columns:

| # | Column | Type | Notes |
|---|---|---|---|
| 1 | `inchikey` | 27-char string | The primary key. **Not** tautomer-normalised |
| 2 | `cid` | int or empty | PubChem CID; present on 68,693,295 rows |
| 3 | `zinc_id` | string or empty | ZINC20 ID; present on 30,187,323 rows |
| 4 | `smiles` | string | Structure |
| 5 | `chemont_tree_json` | JSON array of 5 | `[kingdom, superclass, class, subclass, direct_parent]` as **numeric** IDs |
| 6 | `chemont_other_json` | JSON object | `intermediate_nodes`, `alternative_parents`, `substituents`, `mapped_features`, `geometric_descriptor` |

A real row:

```
inchikey            RDHQFKQIGNGIED-UHFFFAOYSA-N
cid                 1
zinc_id             (empty)
smiles              CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C
chemont_tree_json   [0,12,3909,324,1095]
chemont_other_json  {"intermediate_nodes":[324],
                     "alternative_parents":[346,1205,1238,1831,2449,3865,3919,3940,4150,4225,4557],
                     "substituents":[...], "mapped_features":[...]}
```

Resolved through the dictionary, that row reads:

| Slot | ID | Name |
|---|---:|---|
| kingdom | 0 | Organic compounds |
| superclass | 12 | Lipids and lipid-like molecules |
| class | 3909 | Fatty Acyls |
| subclass | 324 | Fatty acid esters |
| **direct parent** | 1095 | **Acyl carnitines** |
| intermediate node | 324 | Fatty acid esters |
| alternative parents | 346, 1205, 1238, 1831 … | Dicarboxylic acids and derivatives, Carboxylic acids, Carboxylic acid esters, Carbonyl compounds … |

That is acetylcarnitine, PubChem CID 1.

### B.2 The repeat-padding convention — a schema fact, not documented upstream

The five slots are **not** five independent depths. When a compound's direct parent is
shallower than a slot, the trailing slots **repeat the deepest node**. Measured over all
73.1M rows:

| Pattern | Rows | Meaning |
|---|---:|---|
| `xxxxx` | 39,276,855 | direct parent at depth 5 or below — all five slots distinct |
| `xxxx=` | 22,340,569 | direct parent at depth 4 — slot 5 repeats slot 4 |
| `xxx==` | 11,368,062 | direct parent at depth 3 — slots 4 and 5 repeat slot 3 |
| `xxNNx` | 119,735 | class and subclass unresolved (null) |
| `xNNNx` | 53 | superclass, class, subclass all unresolved |

Missing this convention is not academic: reading the slots as literal depths produced
**7,724 false violations per 50,000 rows** in the first implementation. It is now the
first L1 unit test.

### B.3 The dictionary — `chemont_dictionary.tsv`

4,826 lines = header + 4,824 categories + the root. Five columns:

```
numeric_id  chemont_id          name                parent_numeric_id  parent_name
0           CHEMONTID:0000000   Organic compounds   9999999            Chemical entities
1           CHEMONTID:0000001   Inorganic compounds 9999999            Chemical entities
```

The compact `numeric_id` is what appears in the compound file; `chemont_id` is the
canonical published identifier. The root is `9999999` with a `null` parent.

### B.4 The taxonomy — `ChemOnt_2_1.obo`

The authors' own release, OBO format, `data-version: 2.1`, dated 26 Aug 2016,
`sha256 8616a6ec…51fe22`. 4,825 `[Term]` stanzas = 1 root + 4,824 categories:

```
[Term]
id: CHEMONTID:0000000
name: Organic compounds
def: "Compounds that contain at least carbon atom, excluding isocyanide/cyanide …" []
synonym: "an organic compound" EXACT CHEMONT_TERM []
synonym: "Organic Chemicals" RELATED MeSH_TERM [MESH:D02]
synonym: "organic molecule" EXACT ChEBI_TERM [CHEBI:72695]
is_a: CHEMONTID:9999999 ! Chemical entities
```

Each synonym carries a **scope** (EXACT / NARROW / BROAD / RELATED), a **source tag**
(ChEBI_TERM, MeSH_TERM, LIPIDMAPS_TERM, IUPAC_TERM, UniProt_TERM, CHEMONT_TERM) and a
bracketed **external ID list**. A single line can carry several IDs — which is why
"number of mappings" has three defensible readings — lines, distinct (category, ID)
pairs, or distinct IDs. MeSH gives 966 / 1,165 / 514 respectively, a 20 % spread, so any
mapping count must state which reading it uses.

**The dictionary and the OBO disagree on 4 parent edges**: `Sulfinylamines`, `Ureides`,
`1,2-diacyl-3-O-beta-D-galactosyl-sn-glycerols`, `1-acyl-3-O-beta-D-galactosyl-sn-glycerols`.
Same 4,825 terms, same names — four different parents. That is why every violation count
is reported alongside how many of them touch a drifted edge (73 of 198,959).

### B.5 The paper's supplementary files

| File | Format | Size | Contents |
|---|---|---|---|
| Additional file 1 | .docx | 103 KB | Mathematical definition of ChemOnt; **Table S1** (the label definitions above) |
| Additional file 2 | .docx | 4.9 MB | Figures S1–S3, SMARTS/Markush examples |
| Additional file 3 | .csv | **118 MB** | ChEBI annotation dump: 2,687,046 rows over 84,378 compounds |
| Additional file 4 | .xlsx | 1.9 MB | **The 800-compound test set.** Sheet 1 = CID + SMILES; Sheet 2 = 22,043 rows of assignments with FALSE POSITIVE / FALSE NEGATIVE columns marked by the expert panel |
| Additional file 5 | .xlsx | 87 KB | The 20-compound ChEBI comparison, with an `=AVERAGE(B2:B21)` row the authors left in |
| Additional file 6 | .pdf | 4.3 MB | Figure S4, the text-search example |

Files 3–5 are what made the evaluation claims reproducible at all.

---

## C. What the pipeline produces

### C.1 Sufficient statistics — the design that makes re-analysis free

The scan writes three small files that are **sufficient** to recompute every Phase-1
invariant without touching the 2 GB source again:

| File | Schema | Rows | Purpose |
|---|---|---:|---|
| `path_counts.tsv` | `kingdom, superclass, class, subclass, direct_parent, rows` | 3,631 | every distinct label path with its exact row count |
| `dp_alt_pairs.tsv` | `direct_parent, alternative_parent, rows` | 1,047,341 | the weighted class-pair structure; **this is the Phase-2 overlay** |
| `inode_sets.tsv` | `subclass, direct_parent, intermediate_nodes, rows` | 5,178 | distinct intermediate-node sets |

Example rows:

```
path_counts.tsv    0  264  265  13  2309   1878192
dp_alt_pairs.tsv   2309  3940                1878047
inode_sets.tsv     13  2309  60,347          1838568
```

The compression ratio is the point: 73,105,281 rows collapse to 3,631 distinct paths
(20,000× fewer) because most compounds share a classification.

### C.2 The point cloud — `points.bin`

**292,421,124 bytes = 73,105,281 points × 4 bytes.** No header, no delimiter: a flat
array of `(int16 x, int16 y)` little-endian pairs. World coordinates are recovered as
`x / scale` where `scale` is published in `layout.json` (`34,520.51` for the current
build).

Points are **grouped by class**, and `layout.json` records each class's `off` (point
index) and `cnt`, so the renderer issues one GPU draw call per *visible* class and never
draws the rest.

Within a class, points are placed on a **Vogel (sunflower) spiral**:

```
θ = i · 2.39996…      (the golden angle, π(3−√5))
r = R · √((i+0.5)/n)  (√ keeps areal density uniform)
```

The halo radius is `R = K·√(compounds)` with `K = 6.867×10⁻⁵`, chosen so on-screen density
is **algebraically independent of class size**: `n / (π(K√n·z)²) = 1/(πK²z²)`. Measured, a
1,878,192-compound class and a 601-compound class both render at exactly 30 points/pixel.

**Honest limitation:** position *within* a class carries no meaning — there is no 2D
structural embedding of molecules here. Position *between* classes is the taxonomy layout.

### C.3 `layout.json` — the frontend's single input

3.6 MB. One record per class:

| Key | Meaning | Key | Meaning |
|---|---|---|---|
| `i` | numeric ID | `dp` | compounds with this as direct parent |
| `n` | name | `st` | compounds in this subtree |
| `c` | CHEMONTID | `ip` | compounds with it anywhere in the path |
| `p` | parent ID (−1 = kingdom) | `e` | alternative-parent entropy, bits |
| `d` | depth | `px, py` | layout position, normalised |
| `k` | kingdom (0 organic, 1 inorganic) | `nr` | node radius |
| `s, sc` | superclass ID and name | `hr` | halo radius (point culling) |
| `x` | example compounds | `off, cnt` | byte range in `points.bin` |

Plus `tree` (4,822 dendrogram links), `links` (14,000 confusability edges), `bundles`
(5,200 bundled arcs), and `meta` (counts, `scale`, content `version`).

---

## D. Algorithms, stated precisely

**Fruchterman–Reingold** (the class layout, and the abandoned attempts). Repulsion
`k²/d` between all pairs, attraction `d²/k` along edges, displacement clamped by a
"temperature" that cools each iteration. The clamp is not optional: without it the
displacement diverges, so the generator asserts the layout is finite before writing.

**Radial tree layout** (what shipped). Leaves get an angular slot proportional to
`log₁₀(1+compounds)+0.65`; an internal node spans its children's slots; radius is a
function of depth (`RING = linspace(0,1,12)^0.78`). Zero edge crossings by construction —
which is why it beat the force layout.

**Hierarchical edge bundling** (Holten). For a cross-link `a→b`, take the path
`a → LCA(a,b) → b` through the tree as control points, then relax toward the straight
chord by `β = 0.86`:

```
ctrl = β · path + (1−β) · straight
```

Rendered as quadratic Béziers. This is why the arcs flow along the hierarchy instead of
cutting through the middle.

**Squarified treemap** (attempt 4). Recursively split the remaining rectangle, choosing
the split that minimises worst-case aspect ratio.

**Shannon entropy** (the ambiguity score). For class `d` with alternative-parent counts
`n(d,a)`:

```
H(d) = −Σ_a p(a|d) · log₂ p(a|d),   p(a|d) = n(d,a) / Σ_x n(d,x)
```

High `H` = the system reaches for many different alternatives = ambiguous.

**Permutation test with size matching.** Draw 20,000 random class sets of the same size
as the flagged set, compute the mean entropy of each, and locate the observed value in
that distribution. Because entropy correlates with class size at `r = 0.788`, controls
are drawn from the *same size decile*, which is what turns an uncontrolled `p = 0.0001`
into an honest `p = 0.0101`.

---

## E. The scripts, and what each one costs

| Script | Reads | Writes | Runtime |
|---|---|---|---|
| `chemont.py` | — | — | library: tree model, ancestor/depth queries, the R5/R6/PATH checks |
| `paths.py` | — | — | library: resolves large inputs via `$ATLAS_DATA` → `data/` → dev layout |
| `01_verify_taxonomy.py` | `ChemOnt_2_1.obo` | `taxonomy_claims.json`, `taxonomy_nodes.tsv` | <1 s |
| `02_scan_dataset.py` | the 2 GB `.zst` | `scan.json`, 3 sufficient-statistic TSVs | **6 min 35 s** |
| `03_verify_evaluation.py` | Additional files 3, 4, 5 | `evaluation_claims.json` | ~20 s |
| `05_build_graph_data.py` | the statistics | `graph_data.json` | ~10 s |
| `06_build_pointcloud.py` | `graph_data.json` | (superseded treemap layout) | ~3 min |
| `07_build_network.py` | `graph_data.json` | `network_layout.json`, preview PNG | ~90 s |
| `08_build_network_points.py` | `network_layout.json` | `points.bin` (292 MB), `layout.json` | ~2 min |
| `test_checks.py` | fixtures + sample | — | ~2 s, 31 assertions |

Environment: Python 3.11+ with numpy, matplotlib, openpyxl; the `zstd` CLI (the pipeline
shells out rather than depending on a Python binding); Node 20 and Vite 5 for the
frontend. No cheminformatics toolkit is needed — nothing here parses a molecule.

---

# Part II — Method and results

## 1. The baseline paper and what was checked against it

Djoumbou Feunang, Y. *et al.* (2016) *ClassyFire: automated chemical classification with a
comprehensive, computable taxonomy.* **J Cheminform 8, 61.** doi:10.1186/s13321-016-0174-y

Every quantitative claim in the paper was extracted into a ledger, `docs/paper_claims.md`
— **43 claims**, each quoted with its section of origin and tagged by whether it is
reachable from artefacts obtainable outside the authors' group:

| Group | Claims | Subject |
|---|---:|---|
| C1–C12 | 12 | ChemOnt's structure: category counts, depth, the tree property |
| D1–D8 | 8 | The dictionary: synonyms, ChEBI / LIPID MAPS / MeSH mappings |
| R1–R6 | 6 | The classification algorithm and its stated invariants |
| E1–E9 | 9 | The evaluation: test set, false positives, scores, ChEBI comparison |
| S1–S8 | 8 | Deployment scale: 77M compounds, runtimes, database sizes |

Without such a ledger "we reproduced the paper" is unfalsifiable. With it, every claim has
a verdict and a reason.

### Input verification

No analysis was run on a file whose origin had not been established.

| File | Check | Result |
|---|---|---|
| `chemont_dictionary.tsv` | md5 vs Zenodo API | `312aa8b2…` exact |
| `ChemOnt_2_1.obo.zip` | md5 vs Zenodo API | `0616fc94…` exact |
| `…head20.tsv` | md5 vs Zenodo API | `690bf54c…` exact |
| `vocabulary.json` | md5 vs Zenodo API | `55667d9a…` exact |
| `classyfire_…enriched.tsv.zst` (2.0 GB) | md5 vs Zenodo API | `070b2936…` exact |
| the same archive | `zstd -t` | passes, 24.4 GB decompressed |
| `ChemOnt_2_1.obo` | sha256 of the file unzipped from the verified zip | `8616a6ec…51fe22` |

The OBO's digest is asserted at runtime by `scripts/01_verify_taxonomy.py`, so the
replication refuses to run against a substituted file. Note the OBO reached us via
Zenodo's redistribution rather than from `classyfire.wishartlab.com` directly — a verified
chain, but one link longer than ideal.

---

## 2. Phase 1a — the taxonomy structure

`scripts/01_verify_taxonomy.py` parses the OBO and checks 26 claims, tiered so that
near-vacuous shape checks cannot inflate the score.

**Result: 8 of 11 structural claims reproduce exactly** — the tree property (every
non-root node has exactly one `is_a` parent, no cycles, no dangling references, exactly
one parentless node), 2 kingdoms, 31 superclasses = 26 organic + 5 inorganic, 1,729
subclasses, the 4,146 / 678 organic-inorganic partition, maximum depth 11, and the paper's
own worked example (`Phenylpropanoids and polyketides`: 34 direct children, 273
descendants).

Three discrepancies are findings about the paper, not about the reproduction:

- **The paper's category arithmetic does not close.** It states 4,825 categories "in
  addition to the root", split 4,146 organic + 678 inorganic — but 4,146 + 678 = **4,824**.
  Both halves match the OBO exactly, which pins the total independently, so the headline
  4,825 counts the root as a category. The Zenodo release notes, written by an unrelated
  party, also say 4,824.
- **The paper contradicts itself.** Methods says 678 inorganic categories; the Conclusion
  says 674. The OBO says 678.
- **Every cross-ontology mapping count has drifted upward** — ChEBI 6,014→6,034, LIPID MAPS
  789→829, MeSH 844→865, synonyms 9,012→9,024. The OBO is dated 26 Aug 2016 against a
  4 Nov 2016 publication and the paper says the MeSH mapping was ongoing, so this is
  version drift rather than error.

One category, `CHEMONTID:0001218` (Inorganic isocyanides), carries no ChEBI synonym at
all, against the paper's "Each ClassyFire category has one or more mapped ChEBI terms."

---

## 3. Phase 1d — the published evaluation

Additional file 4 contains the 800-compound test set *and* the seven experts' markings of
every error, 22,043 rows. That makes the evaluation claims checkable.

| Claim | Paper | Observed | |
|---|---|---|---|
| E1 test compounds | 800 | 800 | exact |
| E5 false positives | 17 | 17 | exact |
| E5 false negatives | 13 | 13 | exact |
| E3 assignments/compound | 26.38 | 26.38 | exact |
| E7 **maximum score** | 7067.24 | **7066.93** | 0.004 % residual |
| E9 (six ChEBI metrics) | ~31, ~27, ~45, ~33, ~14, ~94 %, 43.6 % | 30.85, 26.90, 45.40, 32.90, 14.30, 94.5 %, 43.5 % | all match |
| E3 total assignments | 21,102 | 21,101 | −1 |
| E4 distinct categories | 1,308 | 1,305 | −3 |
| E6 compounds per category | 2.6 | **16.17** | irreconcilable |

**E7's scoring scheme was recovered from prose.** The paper gives a score of
7067.04 / 7067.24 but never states the formula. Taking its description literally — "each
category was assigned a normalized weight based on its number of occurrences among the 800
chemical entities", so errors on populated categories are penalised more — the weight
`w(c) = occ(c)/800` reproduces the stated maximum as `Σ occ(c)²/800 = 7066.93`, a residual
of 0.004 %. The penalty half of the formula was not recovered and is reported as
unrecovered.

Two further errors in the paper:

- **E6 is wrong.** "Each category was assigned to an average of 2.6 compounds" is
  unreproducible and contradicted by the paper's own figures: 21,102 ÷ 1,308 = 16.13. The
  observed median is 2 and the mode 1 (473 of 1,305 categories occur exactly once), so 2.6
  resembles a typical-value statistic mislabelled as a mean.
- **E7's percentage drops a digit.** 7067.04 / 7067.24 = 99.997 %, printed as "99.97 %".
- A manuscript typo: the text says CID 46936568 was missed as a "*purine* nucleotide
  sugar"; the supplement says "*Pyrimidine*", which is correct for a cytidine derivative.

**Result: 14 of 22 evaluation claims reproduced.**

---

## 4. Phase 1b/1c — one pass over 73,105,281 rows

### Design

Applying the paper's rules per row is prohibitively slow: raw parsing alone runs at
0.64 M rows/s, and rule evaluation dominates everything else. The pipeline therefore
accumulates **bounded sufficient statistics** in the hot loop and evaluates the rules
afterwards over distinct keys only:

| Statistic | Distinct keys | vs 73.1M rows |
|---|---:|---|
| label paths | 3,631 | 20,000× fewer |
| intermediate-node sets | 5,178 | 14,000× fewer |
| (direct parent, alternative parent) pairs | 1,047,341 | 70× fewer |

**73,105,281 rows in 6 min 35 s, 184,087 rows/s.** The dumped statistics are sufficient to
recompute the invariants under any rule reading, so no re-analysis needs another pass over
the 2 GB source.

### Phase 1b — provenance, not paper reproduction

The Zenodo file is a **2026 third-party, InChIKey-deduplicated aggregation**, not the 77M
population the 2016 paper describes. Verifying it checks the release notes, not the paper —
a distinction the report must keep.

Six of seven figures reproduce exactly: rows, ZINC IDs, both-IDs, unresolved slots, classes
seen, classes total. One does not: **PubChem CID 68,693,295 observed against 68,693,298
stated, −3**, reproducibly.

### Phase 1c — the paper's invariants at full scale

The paper states PATH, R4, R5 and R6 as always true. They are not:

| Invariant | Violating rows | Share |
|---|---:|---:|
| PATH — label is a root-path of its direct parent | 198,959 | 0.272 % |
| — with a complete, non-null path | **79,171** | 0.108 % |
| R5 — alternative parents unrelated to the direct parent | 596,367 * | 0.816 % |
| R6 — intermediate nodes | 85,036 | 0.116 % |

\* an upper bound on distinct rows; over a 200,000-row sample, 1,661 pair instances against
1,657 distinct rows (+0.24 %).

Three controls make these robust rather than artefacts:

1. **Taxonomy drift does not explain them.** Only **73** of the 198,959 path violations
   touch any of the 4 edges where the release and ChemOnt 2.1 disagree.
2. **R5 fails under both of the paper's own conflicting definitions.** The Results section
   says alternative parents have no *ancestor–descendant* relation to the direct parent;
   Table S1 says the weaker *parent-child*. **69 of the 93** violating pairs are direct
   parent-child, so they breach both. The dominant case, `Aryl thioethers` →
   `Alkylarylthioethers` (551,743 rows), also breaks R4: category reduction should have
   promoted the more specific child.
3. **It contradicts the release's own headline.** Zenodo states "the canonical five-tier
   hierarchical path is now correct for every row"; 79,171 rows carry a complete, non-null
   path that is still not a root-path.

---

## 5. Verification

`scripts/test_checks.py`, 31 assertions in three tiers:

- **L1 unit** — hand-built toy taxonomies with known answers, no data required.
- **L2 fixture** — the real ChemOnt against facts taken from the paper.
- **L3 sample** — the streaming scan cross-checked against an **independent naive
  re-implementation** of the same rules over 20,000 rows.

A skipped tier exits 2 rather than 0, so an environment missing the large inputs cannot
report a pass on L1 alone.

---

## 6. The visualization

Two design constraints, both measured rather than assumed.

**Density.** 73.1M points on a ~2M-pixel screen puts 550–1,500 points on every occupied
pixel; an early layout measured 1,494 points on the median occupied pixel with 77 % of them
past saturation. Point rendering is therefore only meaningful past a zoom threshold, and
the halo radius scales as `K√n` so that on-screen density is algebraically independent of
class size — a 1,878,192-compound class and a 601-compound class both render at 30
points/pixel.

**Topology.** 4,822 of the edges form a strict tree, and force-directed layouts produce a
hairball on tree-structured data with high-degree hubs. A **radial dendrogram** has an exact
layout with no edge crossings — radius from depth, angle from subtree extent — and the
confusability links are **bundled along the hierarchy** (Holten, `β = 0.86`), which is the
standard technique for a hierarchy plus cross-links.

The frontend renders two layers on one camera: a WebGL2 point cloud of all 73,105,281
compounds streamed into a pre-allocated 292 MB GPU buffer, and a 2D overlay carrying the
4,824 classes, 4,822 dendrogram links and 5,200 bundled arcs. Because the cloud is 292 MB
on an 8 GB machine it is never held twice: the binary is streamed in and each chunk dropped
after upload, and the renderer issues one draw call per *visible* class using per-class byte
ranges.

---

## 7. Shipping

**Repository.** 30 commits, ordered by how the project builds up rather than
chronologically. 71 tracked files, kept small because `.gitignore` excludes the 292 MB
point cloud, the 118 MB supplementary CSV, the 12 MB pairs file and all generated assets.

**Deployment.** Built locally (server and build machine are both aarch64); only the
compiled app (176 KB), `layout.json` (3.5 MB) and `points.bin` (292 MB) ship. Both data
files are verified by sha256 on both ends.

**Caching.** The cloud is immutable, so its URL carries a content hash
(`points.bin?v=…`) and is served `max-age=31536000, immutable` — a returning visitor never
re-fetches 292 MB. `layout.json` is deliberately not immutable, since it carries that
version; it revalidates by ETag. Range requests stay enabled so an interrupted stream
resumes, gzip is off for the binary (it compresses only 14 %) and on for the JSON
(3.67 MB → 973 KB).

---

## 8. What is not reproduced

- **ClassyFire was never run.** The paper is fundamentally about a classification program —
  four algorithmic steps, >9000 SMARTS patterns, Markush structures, ~200 IUPAC regexes.
  No molecule was classified here. What is reproduced is the taxonomy *artefact* and the
  *bookkeeping* of the evaluation.
- **R1–R3** (the rule base) — never released.
- **E2** (249.9 s on their hardware) and **S1–S6, S8** (77M compounds, the 424-hour PubChem
  run, 550 ms/compound, 13,100 descriptions) — require their servers.
- **The R5 figure is an upper bound.** 596,367 counts pair instances, not distinct rows;
  the true distinct-row count is roughly 1,400 lower.
- **R5's alt-vs-alt clause is not evaluated** by the scan, because the sufficient statistics
  do not retain per-row alternative-parent sets. Over 200,000 rows, zero rows violate that
  clause without also violating dp-vs-alt, so no reported figure depends on it.

**Scorecard: ~25 of ~41 attemptable claims reproduced, 16 measured and mismatched with
explanations, ~11 structurally impossible.**

---

*Errors encountered during development, and how they were resolved, are recorded separately
in [`process_notes.md`](process_notes.md) — that file is the evidence for the course's
Declaration of Tool Usage, not part of these results.*
