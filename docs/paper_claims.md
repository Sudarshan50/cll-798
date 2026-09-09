# Baseline paper — ledger of reproducible quantitative claims

Djoumbou Feunang, Y., Eisner, R., Knox, C., Chepelev, L., Hastings, J., Owen, G., Fahy, E.,
Steinbeck, C., Subramanian, S., Bolton, E., Greiner, R., Wishart, D. S. (2016).
*ClassyFire: automated chemical classification with a comprehensive, computable taxonomy.*
Journal of Cheminformatics 8, 61. doi:10.1186/s13321-016-0174-y

Every claim below is quoted from the paper. Each is tagged with whether it is
reproducible from artefacts we hold, and if so from which one.

Artefacts held locally (all verified against publisher / Zenodo checksums):
- `data/paper_supplementary/` — **the paper's own six Additional files**, downloaded from
  the publisher. These make E1, E3–E9 and S7 reproducible; an earlier revision of this
  ledger wrongly recorded them as out of reach.
- `data/ChemOnt_2_1.obo` — the ChemOnt taxonomy in OBO format, `data-version: 2.1`,
  `date: 26:08:2016`, `auto-generated-by: OBO-Edit 2.3.1`. This is the file the paper
  describes in Methods/Component 1 ("The resulting OBO file was generated with OBO-Edit
  [34], and can be downloaded from the ClassyFire website").
- (Phase 1b) the Zenodo v3 InChIKey-deduplicated ClassyFire label collection.

Legend for **Reproducible?**
- **OBO** — checkable exactly from the OBO file.
- **DATA** — checkable from the Zenodo compound-level release (Phase 1b).
- **SUP** — reproducible from the paper's Additional files (3 = ChEBI dump, 4 = the
  800-compound test set with the expert panel's markings, 5 = the ChEBI comparison).
  Verified by `scripts/03_verify_evaluation.py`.
- **NO** — genuinely out of reach: depends on the closed ClassyFire codebase, the ChemAxon
  JChem licence, the private MySQL store, or their hardware. Not attempted, and not
  attemptable by anyone outside the authors' group.

**On E7.** The paper states a score of 7067.04 out of 7067.24 but never gives the scoring
formula. Taking its description literally — "each category was assigned a normalized weight
based on its number of occurrences among the 800 chemical entities", so that errors on
populated categories are penalised more — the weight w(c) = occ(c)/800 reproduces the stated
maximum as sum_c occ(c)^2/800 = 7066.93, a residual of 0.004 %. That is reported as a
*reconstruction with its residual*, never as the authors' method. The penalty half of the
formula is not recovered.

---

## C. Taxonomy structure  (Methods → Component 1; Conclusion)

| # | Claim (paper's words) | Value | Reproducible? |
|---|---|---|---|
| C1 | "a taxonomy consisting of >4800 different categories" (Abstract) | > 4800 | OBO |
| C2 | "this extensive chemical taxonomy contains a total of 4825 chemical categories of organic (4146) and inorganic (678) compounds, in addition to the root category (Chemical entities)" | 4825 = 4146 + 678, plus root | OBO |
| C3 | "The top level (Kingdom) partitions chemicals into two disjoint categories: organic compounds versus inorganic compounds" | 2 kingdoms | OBO |
| C4 | "SuperClasses (which includes 26 organic and 5 inorganic categories)" | 31 = 26 + 5 | OBO |
| C5 | "the Class level, which now includes 764 nodes" | 764 | OBO |
| C6 | "There are 1729 SubClasses in the current taxonomy." | 1729 | OBO |
| C7 | "there are 2296 additional categories below the SubClass level covering taxonomic levels 5–11" | 2296 | OBO |
| C8 | "can be represented as a tree with a maximum depth of 11 levels" | max depth 11 | OBO |
| C9 | "and an average depth of five levels per node" | mean depth ≈ 5 | OBO |
| C10 | "The main relationship type connecting these different categories is the 'is_a' relationship... arranged in a tree structure" — i.e. every node has at most one parent, no cycles, one root | strict tree | OBO |
| C11 | "we have also developed a comprehensive taxonomy for inorganic compounds consisting of 674 categories" (Conclusion) | 674 | OBO — **note: conflicts with C2's 678** |
| C12 | "The 'Phenylpropanoids and polyketides' category currently has 34 direct children and a total of 273 descendant categories" | 34 / 273 | OBO |

## D. Dictionary, synonyms and cross-ontology mappings  (Methods → Component 2; Mapping section)

| # | Claim | Value | Reproducible? |
|---|---|---|---|
| D1 | "A total of 9012 English synonyms were added to the ChemOnt terminology data set." | 9012 | OBO |
| D2 | "over 9000 synonyms" (Conclusion) | > 9000 | OBO |
| D3 | "A total of 6014 category mappings were created, with an average of 1.24 ChEBI synonyms per category." | 6014 / 1.24 | OBO |
| D4 | "Each ClassyFire category has one or more mapped ChEBI terms." | 100 % ChEBI coverage | OBO |
| D5 | "A total of 789 ClassyFire categories were mapped to one of 307 LIPID MAPS terms each." | 789 / 307 | OBO |
| D6 | "844 ClassyFire categories have been mapped to at least one corresponding MeSH term, accounting for a total of 945 mappings" | 844 / 945 | OBO |
| D7 | "each term was formally defined using a precise, yet easily understood text description" (4825 terms) | every term has a `def:` | OBO |
| D8 | "For any ChemOnt term, a synonym can have the identical meaning (exact scope), a more specific meaning (narrow scope), or a less specific meaning (broad scope)... otherwise related" | scopes ⊆ {EXACT, NARROW, BROAD, RELATED} | OBO |

## R. Rule base and algorithm  (Methods → Component 3; classification process)

| # | Claim | Value | Reproducible? |
|---|---|---|---|
| R1 | "Converting the 4825 definitions in our Chemical Classification Dictionary led to the creation of >9000 SMARTS strings." | > 9000 | NO — rule base not released |
| R2 | "The Markush patterns used by ClassyFire constitute only about 4% of the set of patterns" | ~4 % | NO |
| R3 | "a module, which uses a set of ~200 regular expressions" for IUPAC-name parsing | ~200 | NO |
| R4 | "if there is a parent–child relationship... only the child node is retained" (category reduction) | algorithmic | NO (but testable as a *consistency property* of released labels — Phase 1c) |
| R5 | "Alternative parents are categories that describe the compound but do not have an ancestor–descendant relationship with each other or with the Direct Parent." | algorithmic | DATA — testable against the released alternative-parent lists |
| R6 | "Intermediate Nodes... are descendants of a subclass (any category with a depth of 4), but have a depth lower than the direct parent." | definition | OBO + DATA |

## E. Evaluation  (Training and evaluation; Results)

| # | Claim | Value | Reproducible? |
|---|---|---|---|
| E1 | test set of 800 unique structures from DrugBank / LIPID MAPS / HMDB / T3DB | 800 | **SUP4** |
| E2 | "The classification process took 249.9 s on... 4 CPU CentOS nodes, with 3.6 GB of RAM... maximum of 16 threads." | 249.9 s | NO |
| E3 | "A total of 21,102 category assignments were made, for an average of 26.38 assignments per compound." | 21102 / 26.38 | **SUP** (21102/800 = 26.3775 — self-consistent) |
| E4 | "ClassyFire assigned a total of 1308 distinct Categories." | 1308 | **SUP** |
| E5 | "a total of 17 false positives (out of 21,102 assignments)"; "13 missing assignments (false negatives)" | 17 / 13 | **SUP** |
| E6 | "Each category was assigned to an average of 2.6 compounds." | 2.6 | **SUP** (21102/1308 = 16.1; 800·26.38/1308 = 16.1 — **does not equal 2.6**; flagged) |
| E7 | "ClassyFire obtained score of 7067.04, or 99.97% of a maximum score of 7067.24." | 99.97 % | **SUP** (7067.04/7067.24 = 99.997 % — **paper states 99.97 %**; flagged) |
| E8 | "precision of 99.8% and a recall of 99.9%" | — | **SUP** |
| E9 | ChEBI comparison on 20 compounds: ~33 ChEBI classes/compound, ClassyFire ~31 categories, mapping → ~27 terms, extended → ~45, ~14 missing from ChEBI, >2 ChEBI terms unreachable, "reproduce ~94% of the ChEBI annotations... increase the number of annotations by another 43.6%" | — | **SUP** |

## S. Scale of deployment  (Abstract; Results; Use cases)

| # | Claim | Value | Reproducible? |
|---|---|---|---|
| S1 | "ClassyFire has been used to annotate over 77 million compounds" | > 77 M | DATA (partial — see note) |
| S2 | "more than 70 million compounds to date" stored in the local MySQL DB | > 70 M | DATA (partial) |
| S3 | ">60,000,000 molecules in PubChem" classified | > 60 M | DATA (partial) |
| S4 | PubChem run "completed in 424 h for an average of 550 ms per compound" | 424 h / 550 ms | NO |
| S5 | "<50 ms" cached lookup, "average of 540 ms" for a novel compound | — | NO |
| S6 | "only 0.12% of the >91,000,000 compounds (as of June 2016)" in PubChem carry a MeSH class | 0.12 % | NO |
| S7 | ">6000 molecules in DrugBank, >25,000 in LIPID MAPS, >42,000 in HMDB, >43,000 in ChEBI" | — | **SUP3** (ChEBI part only) |
| S8 | "over 13,100 meaningful, 20–50 word descriptions" generated | 13100 | NO |

**Note on S1–S3.** The paper's 77 M is a count of *structures processed by the ClassyFire
server as of June 2016*, not a released file. The Zenodo v3 collection we verify in Phase 1b
is a 2025 third-party, InChIKey-deduplicated redistribution of ClassyFire labels. It is
therefore evidence *about* the same classification effort at a later date and after
deduplication — it is not the same population, and we do not claim it reproduces S1.
This distinction is stated explicitly wherever the numbers are compared.
