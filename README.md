# Hierarchical Topology of Chemical Space — reproducing ClassyFire/ChemOnt

CLL798D / CHL7903 Network Science · Group 2 · Pathway 1 (reproduce, then extend)
Prem Bhugra (2022CH71038) · Sakhare Yash Balram (2022CH71496) · Sudarshan Kumar Oraon (2022CH71511)

**Baseline paper.** Djoumbou Feunang, Y. *et al.* (2016) *ClassyFire: automated chemical
classification with a comprehensive, computable taxonomy.* J Cheminform **8**, 61.
doi:[10.1186/s13321-016-0174-y](https://doi.org/10.1186/s13321-016-0174-y)

**Dataset.** InChIKey-Deduplicated ClassyFire/ChemOnt Label Collection, v3, Zenodo
DOI:[10.5281/zenodo.20472700](https://zenodo.org/records/20472700), CC BY 4.0.

## Data sources, licensing and attribution

This repository **redistributes** third-party material. Everything below is other
people's work, reused under the licence stated, and none of it is ours.

| Material in this repo | Source | Licence | Verified |
|---|---|---|---|
| `data/ChemOnt_2_1.obo` | ChemOnt 2.1 taxonomy by Djoumbou Feunang *et al.*, obtained from the Zenodo record below | CC BY 4.0 (as redistributed by that record) | md5 of the source zip against the Zenodo API; `sha256 8616a6ec…51fe22` asserted at runtime |
| `data/chemont_dictionary.tsv`, `data/vocabulary.json`, `…head20.tsv` | Zenodo `10.5281/zenodo.20472700` | CC BY 4.0 | md5 against the Zenodo API |
| `data/paper_supplementary/*.xlsx` | Additional files 4 and 5 of the baseline paper, Springer Nature | CC BY 4.0 — the article's own licence, [creativecommons.org/licenses/by/4.0](http://creativecommons.org/licenses/by/4.0/), confirmed via Crossref | downloaded from the publisher |
| the 2 GB compound file *(not committed)* | Zenodo `10.5281/zenodo.20472700` | CC BY 4.0 | md5 against the Zenodo API; `zstd -t` |

**Not redistributed here**, and deliberately so: the 2.0 GB compound archive, the 118 MB
ChEBI supplement (Additional file 3), and every generated artefact. `scripts/paths.py`
resolves the large inputs at run time and names the DOI if they are missing.

**What is ours**, and the only thing you may reuse without consulting the sources above:
the code in `scripts/` and `frontend/src/`, the analysis in `results/`, and the documents
in `docs/`.

**ClassyFire itself was never run here.** No molecule was classified by this project; we
reproduce the released taxonomy and the published evaluation bookkeeping. See
[`docs/how_this_was_built.md`](docs/how_this_was_built.md) §11 for exactly what is and is
not reproduced.

**Citing the originals.** Any use of this work should cite the baseline paper and the
dataset, not this repository:

```bibtex
@article{djoumbou2016classyfire,
  title   = {ClassyFire: automated chemical classification with a comprehensive,
             computable taxonomy},
  author  = {Djoumbou Feunang, Yannick and Eisner, Roman and Knox, Craig and
             Chepelev, Leonid and Hastings, Janna and Owen, Gareth and Fahy, Eoin and
             Steinbeck, Christoph and Subramanian, Shankar and Bolton, Evan and
             Greiner, Russell and Wishart, David S.},
  journal = {Journal of Cheminformatics},
  volume  = {8}, pages = {61}, year = {2016},
  doi     = {10.1186/s13321-016-0174-y}
}
```

## Phase structure

| Phase | Scope | Status |
|---|---|---|
| **1a** | Reproduce the paper's ChemOnt **taxonomy structure** claims from the released OBO | done, self-audited — `scripts/01_verify_taxonomy.py` |
| **1b** | Verify the Zenodo release's **scale and coverage** claims by one streaming pass | done — `scripts/02_scan_dataset.py` |
| **1c** | Test the paper's stated invariants (PATH / R5 / R6) against all 73.1M labels | done — same pass |
| **1d** | Reproduce the paper's **evaluation** claims from its own supplementary files | done — `scripts/03_verify_evaluation.py` |
| **2** | Confusability overlay: are the paper's expert-flagged errors where ClassyFire hedges? | done, pre-registered — `scripts/10`–`13c`, [`docs/phase2_results.md`](docs/phase2_results.md) |
| **E** | Bipartite extension: compound × class and compound × functional-group networks, their one-mode projections, statistically validated backbones | done — `scripts/09`, `14`–`16`, [`docs/bipartite_extension.md`](docs/bipartite_extension.md) |

Phase 1 makes no claim that is not stated in the baseline paper or the dataset's own
release notes. Every target value is quoted, with its section of origin, in
[`docs/paper_claims.md`](docs/paper_claims.md).

## Layout

```
data/     ChemOnt_2_1.obo          taxonomy as released by the ClassyFire authors
docs/     paper_claims.md          ledger of every quantitative claim in the paper
docs/     phase1a_audit.md        self-audit of the replication logic (6 flaws, fixed)
scripts/  01_verify_taxonomy.py    Phase 1a
results/  taxonomy_claims.json     machine-readable verdicts
          taxonomy_nodes.tsv       one row per ChemOnt node: parent, depth, kingdom, fan-out
```

## Reproducing Phase 1a

```bash
python3 scripts/01_verify_taxonomy.py
```

Standard library only. Runs in under a second.

## Phase 1a headline result (v2, post-audit)

The ChemOnt 2.1 OBO holds **4825 `[Term]` stanzas = 1 root + 4824 categories**, and is a
**strict tree** (every non-root node has exactly one `is_a` parent, no cycles, no dangling
references, exactly one parentless node), with **maximum depth 11** — as the paper asserts.

Claims are tiered, because a file of roughly the right shape satisfies ">4800 categories"
regardless of whether the replication worked:

| tier | reproduced exactly |
|---|---|
| **STRUCTURAL** — single-definition counts the paper's argument rests on | **8 / 11** |
| **BOOKKEEPING** — counts whose definition the paper leaves ambiguous | **0 / 5** |
| **TRIVIAL** — near-vacuous shape checks | 3 / 3 |

Excluded from scoring: one claim that holds only under rounding (C9, mean depth 4.636 vs
"five levels per node") and one ratio derived from an already-failing count (D3-avg).

**What reproduced exactly:** strict tree; 2 kingdoms; 31 SuperClasses = 26 organic + 5
inorganic; 1729 SubClasses; the 4146/678 organic-inorganic partition; maximum depth 11;
the "Phenylpropanoids and polyketides" worked example at 34 direct children and 273
descendants; and a text definition on every one of the 4824 categories.

**What did not, in three groups:**

1. **The paper's category arithmetic does not close.** It claims 4825 categories "in
   addition to the root", split 4146 organic + 678 inorganic — but 4146 + 678 = 4824.
   The observed partition matches *both* halves exactly, which pins the total at 4824
   independently, so the headline 4825 counts the root as a category. Its level tally is
   separately one short twice: 764 Classes and 2296 below-SubClass, against 765 and 2297
   observed. With the observed figures the tree closes exactly: 2 + 31 + 765 + 1729 + 2297
   = 4824. This +1 is **unresolved** — settling drift vs manuscript slip requires the
   paper's own Additional file 1, Table S1.
2. **The paper contradicts itself once.** Methods says 678 inorganic categories, the
   Conclusion says 674. The OBO gives 678. This is one measurement, reported once.
3. **Every cross-ontology mapping count has drifted upward, consistently and in one
   direction** — ChEBI 6014→6034, LIPID MAPS 789→829 categories over 307→309 terms,
   MeSH 844→865 categories over 945→966 mappings, synonyms 9012→9024. The paper describes
   these mappings as ongoing ("This MeSH mapping will likely continue for another year or
   two"), and the OBO is dated 26 Aug 2016 against a 4 Nov 2016 publication date, so the
   released file is a later snapshot than the manuscript was written against. One category,
   `CHEMONTID:0001218` (Inorganic isocyanides), carries no ChEBI synonym at all, against the
   paper's "Each ClassyFire category has one or more mapped ChEBI terms."

Nothing here undermines the paper's substantive result. The taxonomy is exactly the object
the paper describes; the discrepancies are bookkeeping in the manuscript and version drift
in the mapping tables, and both are quantified rather than smoothed over.

The replication logic was itself audited and corrected — six flaws found, two of which
changed reported verdicts. See [`docs/phase1a_audit.md`](docs/phase1a_audit.md).


---

## Phase 1b — provenance (Zenodo release notes, not paper claims)

73,105,281 rows read in 395 s (185k rows/s). **Six of seven figures reproduce exactly.**

| figure | stated | observed | delta |
|---|---:|---:|---:|
| rows | 73,105,281 | 73,105,281 | 0 |
| with PubChem CID | 68,693,298 | 68,693,295 | **−3** |
| with ZINC20 ID | 30,187,323 | 30,187,323 | 0 |
| with both | 25,775,344 | 25,775,344 | 0 |
| unresolved slot rows | 119,795 | 119,795 | 0 |
| classes seen | 4,120 | 4,120 | 0 |
| classes total | 4,824 | 4,824 | 0 |

The −3 on PubChem CID is reproducible and matches an earlier independent scan exactly.
Note also that the release's own denominator is **4,824**, not the paper's 4,825 —
independent corroboration of the Phase 1a finding.

## Phase 1c — the paper's invariants, tested on all 73.1M labels

The paper asserts these as always true. They are not.

| invariant | violating rows | share |
|---|---:|---:|
| PATH — label is a root-path of its direct parent | 198,959 | 0.272 % |
| — of which genuinely wrong ancestors (not just null slots) | **79,171** | 0.108 % |
| R5 — alternative parents unrelated to the direct parent | 596,367 * | 0.816 % |
| — under Table S1's weaker *parent-child* reading | 589,073 * | 0.806 % |
| R6 — intermediate nodes (inclusive reading) | 85,036 | 0.116 % |
| R6 — strict Table S1 reading | 9,589 | 0.013 % |

\* upper bound on distinct rows. Rows carrying two violating pairs do occur: over a
200,000-row sample, 1,661 pair instances against 1,657 distinct rows (+0.24 %), so the
true distinct-row count is roughly 1,400 below the figure above.

Three points make these findings robust rather than artefacts:

1. **Taxonomy drift does not explain them.** Only **73** of the 198,959 path violations
   touch any of the 4 edges where the release and ChemOnt 2.1 disagree.
2. **R5 fails under both of the paper's own conflicting definitions.** The main text says
   alternative parents have no *ancestor-descendant* relation to the direct parent;
   Table S1 says the weaker *parent-child*. **69 of the 93** violating pairs are direct
   parent-child, so they breach both. The dominant case is
   `Aryl thioethers` → `Alkylarylthioethers` (551,743 rows), which also breaches R4's
   category-reduction rule: the more specific child should have become the direct parent.
3. **It contradicts the release's own headline.** The Zenodo notes state "the canonical
   five-tier hierarchical path is now correct for every row"; **79,171 rows carry a
   complete, non-null path that is still not a root-path** of their direct parent. That
figure is emitted directly by the scan (`path_violation_rows_fully_populated`); deriving
it as 198,959 − 119,795 mixes two populations, because the 7 rows with a null direct
parent are unresolved but never counted as path violations.

## Phase 1d — the paper's evaluation claims (14/22)

Reproduced exactly from the authors' supplementary files: the 800-compound test set,
**17 false positives and 13 false negatives**, 26.38 assignments per compound, and all six
ChEBI-comparison metrics (~31, ~27, ~45, ~33, ~14, 94.5 %, 43.5 %).

**E7's scoring scheme was recovered from prose.** The paper gives a score of 7067.04/7067.24
but never the formula. Weighting each category by its occurrence count, `w(c) = occ(c)/800`,
gives a maximum of `sum_c occ(c)^2/800 = 7066.93` — a **0.004 % residual**. The penalty half
is not recovered, and is reported as unrecovered.

Further errors found in the paper: **E6** ("average of 2.6 compounds per category") is
unreproducible and contradicted by the paper's own 21,102 / 1,308 = 16.13; **E7's
percentage** drops a digit (7067.04/7067.24 = 99.997 %, printed as 99.97 %); and the text
says CID 46936568 was missed as a "*purine* nucleotide sugar" where its own supplement says
"*Pyrimidine*" — the supplement is right for a cytidine derivative.

## Re-analysis without rescanning

The scan writes **bounded sufficient statistics**, so no change of rule, definition or
taxonomy version requires another pass over the 2 GB file:

| file | contents | size |
|---|---|---|
| `results/path_counts.tsv` | all 3,631 distinct label paths with row counts | 85 KB |
| `results/dp_alt_pairs.tsv` | all 1,047,341 (direct parent, alternative parent) pairs | 13 MB |
| `results/inode_sets.tsv` | all 5,178 distinct intermediate-node sets | 89 KB |

These are sufficient to recompute every Phase-1 invariant exactly. `--pairs` additionally
emits alternative-parent co-occurrence pairs, the Phase-2 class-class network input.

## Testing

```bash
python3 scripts/test_checks.py      # 31 tests, ~2 s, no dependencies
make test PY=.venv/bin/python       # + test_phase2.py (26) and test_bipartite.py (16)
```

L1 unit tests on toy taxonomies (each a regression for a bug actually made), L2 fixtures
against the real ChemOnt, L3 agreement between the streaming scan and an independent naive
re-implementation over the sample. The extension suites test every estimator on toys with a
known answer (BiCM against exhaustive enumeration, closed-form clique projections against
networkx, a planted over-expressed pair, a perfectly nested matrix) and assert that the
Phase-1 outputs are byte-identical before and after.

---

## Extension — the release as a bipartite network

ChemOnt shows each compound in one box, but ClassyFire also attaches ~12.6 *alternative
parents* and ~9.3 functional groups to it. The extension reads the release as what it is:
an affiliation network of 73,105,281 compounds × 3,851 classes (995,524,038 incidences) and
73.1M compounds × 396 functional groups (679,126,726 incidences).

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
make extension PY=.venv/bin/python   # ~4 min scan + ~10 min analysis, then copies the UI data
```

`scripts/09_scan_bipartite.py` is the only step that reads the 2 GB file (one parallel pass,
~4 min). It writes the class-mode projection `C = BᵀB`, Newman's weighted projection, the
direct-parent × substituent matrix and both degree sequences, and asserts the identities that
tie them together (`Σ k_c = Σ L·n_L`, `Σ_{c<c'} C = Σ C(L,2)·n_L`). Its totals reconcile
exactly with Phase 1b: 73,105,281 rows, 922,418,764 alternative-parent labels.

**One-mode projections, and why the naive ones mislead.**
- Per hierarchy level the compound projection is a disjoint union of cliques, so its
  clustering is 1 by construction; its density, from Kingdom to direct parent, falls from
  0.998 to 0.0053, exactly, without building a ~10¹⁵-edge graph.
- Joined on any shared label, 99.2 % of compound pairs are linked; removing the 3 hub labels
  (Hydrocarbon derivatives, Organic oxides, Azacyclic compounds) drops that to 65 %.
- 24.3M distinct label sets: 79.8 % of compounds are indistinguishable from another in the
  projection.

**Validated class backbone.** The raw class projection has 1,393,402 edges (density 0.188).
Against the bipartite configuration model (BiCM, Saracco et al. 2015/2017), solved exactly
over compound-degree classes, with BH-FDR 0.01 and an effect-size floor of observed ≥ 2 ×
expected, **174,525 edges survive**. Their clustering is 2.8× a degree-preserving null.
Louvain communities on that backbone agree only moderately with the SuperClasses (NMI 0.41,
ARI 0.32). On bipartite modularity (Guimerà et al. 2007), which scores partitions on the
bipartite data itself rather than on the projection, they score 0.074 against 0.042 for the
official SuperClass partition. The strongest surviving edges are near-synonyms (Aryl chlorides
– Organochlorides, Anisoles – Methoxybenzenes): genuine redundancy in the alternative-parent
layer.

**Class × functional group.** The taxonomy orders functional-group chemistry only weakly
(Mantel r = −0.20 between profile similarity and tree distance, p = 0.001). The matrix is no
more nested than its margin-preserving null. Reactive groups cut across the skeleton classes:
nitro groups occur in 1,960 classes. A class-level hazard screen would miss most compounds
that carry them.

## Phase 2 — the confusability overlay (pre-registered)

The plan is [`docs/phase2_plan.md`](docs/phase2_plan.md). The pre-registration
[`docs/phase2_prereg.md`](docs/phase2_prereg.md) was frozen by sha256 before the single
primary test ran; both hashes are recorded in `results/phase2/validation.json`.

| hypothesis | result | verdict |
|---|---|---|
| **H1a**: the 21 expert-flagged classes are more ambiguous than size-matched controls | p = 0.453, effect +0.002 bits, 95 % CI [−0.23, 0.33] | **not supported** |
| **H1b**: overlay clustering and modularity exceed a degree-preserving null | z = 182 (weighted clustering), z = 134 (modularity) | **supported** |
| **H1c**: overlay communities are not a relabelling of the SuperClasses | NMI 0.64, but every community lies inside one SuperClass by construction | **inconclusive** |

The negative H1a is the substantive finding. The unfiltered pilot reproduces (p = 0.0102 against
the brief's 0.0101). But the pre-registered tree-proximity filter leaves 86.7 % of classes,
including both positive controls, with no competitor at all. A tree-free PMI score gives
p = 0.0001, which is reported as exploratory because it was not pre-registered. None of the
36 sensitivity-sweep cells reaches p < 0.05. **Monte Carlo error propagation**: swapping a
compound's direct parent for a statistically likely competitor changes its SuperClass in
54.8 % of draws, against 85.9 % for a uniformly drawn wrong class. Errors guided by the
classifier's own hedging stay local.

## Interactive view

`make dev`, then choose **Extension**. The map draws either the Phase-2 competing overlay or
the validated BiCM backbone over the radial taxonomy, coloured by taxonomy, by either
community partition or by ambiguity, with the expert-flagged classes ringed. The side panel
has the results, charts and method notes for both analyses. See
[`README_FRONTEND.md`](README_FRONTEND.md). The literature survey behind the method choices
is [`docs/extension_research.md`](docs/extension_research.md).
