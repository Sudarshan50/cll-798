# Process notes — errors observed and how they were resolved

Split out of `how_this_was_built.md`, which is now a technical account. This file exists
because the course's Declaration of Tool Usage asks what errors, limitations and gaps were
observed when using automated tools, and how the outcome was refined. It is the evidence
for that answer, not part of the project's results.

# Part II — The narrative: how it was actually built

## 0. Reading the brief

Two PDFs defined the task: the group's approved abstract and the course policy.

The machine had no `pdftotext` and no Python PDF library, and the system Python refused
`pip install` (PEP 668). Resolution: a throwaway virtualenv in the scratch directory
with `pypdf`, used only to extract text. This is worth recording because it is the first
of several places where the obvious tool was not available and a substitute was needed.

What the policy required, and what it means in practice:

| Requirement | Consequence for this project |
|---|---|
| Pathway 1: reproduce **and extend** a baseline paper | Reproduction alone does not satisfy the brief |
| "replicating its core statistical properties (fully or partially)" | Partial reproduction is acceptable if honest about scope |
| All assumptions "MUST be clearly stated and justified" | Every definitional choice needs to be written down |
| Declaration of Tool Usage, including LLM use | This document exists partly to serve that |

The abstract named the baseline paper (Djoumbou Feunang et al. 2016, ClassyFire/ChemOnt,
*J Cheminform* **8**:61) and the dataset (Zenodo `10.5281/zenodo.20472700`).

**A finding before any work started.** A sibling directory held a substantial earlier
attempt at the same project — a completed 73.1M-row scan, seven scripts, figures, a
React dashboard and a `report.pdf`. This was surfaced immediately rather than silently
duplicated or silently overwritten. The decision taken was a clean rebuild in a fresh
directory, reusing only verified input files.

---

## 1. Getting the paper, and turning it into something checkable

The paper is open access. It was fetched from PubMed Central (PMC5096306) and converted
to text with `html2text`.

Reading it once is not enough to reproduce it. The paper's claims are scattered across
Abstract, Methods, Results, Use cases and Conclusion, and they are of very different
kinds: some are properties of a released file, some are counts from a private database,
some are the output of software nobody else can run.

So every quantitative claim was extracted into a ledger, `docs/paper_claims.md` —
**43 claims**, each quoted verbatim with its section of origin, and each tagged by
whether it is reachable from artefacts obtainable outside the authors' group:

| Group | Claims | What they are about |
|---|---|---|
| C1–C12 | 12 | ChemOnt's structure: category counts, depth, the tree property |
| D1–D8 | 8 | The dictionary: synonyms, ChEBI/LIPID MAPS/MeSH mappings |
| R1–R6 | 6 | The classification algorithm and its stated invariants |
| E1–E9 | 9 | The evaluation: test set, false positives, scores, ChEBI comparison |
| S1–S8 | 8 | Deployment scale: 77M compounds, runtimes, database sizes |

This ledger is the single most useful artefact produced. Without it, "we reproduced the
paper" is unfalsifiable. With it, every claim has a verdict and a reason.

---

## 2. Verifying the inputs before trusting them

No analysis was run on a file whose origin had not been established.

| File | Check | Result |
|---|---|---|
| `chemont_dictionary.tsv` | md5 vs Zenodo API | `312aa8b2…` exact |
| `ChemOnt_2_1.obo.zip` | md5 vs Zenodo API | `0616fc94…` exact |
| `…head20.tsv` | md5 vs Zenodo API | `690bf54c…` exact |
| `vocabulary.json` | md5 vs Zenodo API | `55667d9a…` exact |
| `classyfire_…enriched.tsv.zst` (2.0 GB) | md5 vs Zenodo API | `070b2936…` exact |
| the same archive | `zstd -t` integrity | passes, 24.4 GB decompressed |
| `ChemOnt_2_1.obo` | sha256 of the file unzipped from the verified zip | `8616a6ec…51fe22`, chain closed |

The OBO's digest is now asserted at runtime by `scripts/01_verify_taxonomy.py`, so the
replication refuses to run against a substituted file.

**A provenance subtlety worth stating:** the OBO reached us via Zenodo's redistribution,
not from `classyfire.wishartlab.com` directly. That is a verified chain, but it is one
link longer than ideal, and the report should say so.

Later, all six of the paper's supplementary files were downloaded from the publisher
(~131 MB). That turned out to matter enormously — see §5.

---

## 3. Phase 1a — reproducing the taxonomy, then auditing the reproduction

`scripts/01_verify_taxonomy.py` parses the OBO and checks 26 claims.

**What reproduced exactly:** the tree property (every non-root node has exactly one
`is_a` parent, no cycles, no dangling references, exactly one parentless node), 2
kingdoms, 31 superclasses = 26 organic + 5 inorganic, 1729 subclasses, the 4146/678
organic-inorganic partition, maximum depth 11, and the paper's own worked example
(`Phenylpropanoids and polyketides`: 34 direct children, 273 descendants).

**What did not, and why it is a real finding:**

- **The paper's category arithmetic does not close.** It states 4825 categories "in
  addition to the root", split 4146 organic + 678 inorganic. But 4146 + 678 = **4824**.
  Both halves match the OBO exactly, which pins the total independently — so the
  headline 4825 counts the root as a category. The Zenodo release notes, written by an
  unrelated party, also say 4,824.
- **The paper contradicts itself.** Methods says 678 inorganic categories; the
  Conclusion says 674. The OBO says 678.
- **Every cross-ontology mapping count has drifted upward** — ChEBI 6014→6034, LIPID
  MAPS 789→829, MeSH 844→865, synonyms 9012→9024. The OBO is dated 26 Aug 2016 against
  a 4 Nov 2016 publication, and the paper says the MeSH mapping was ongoing. This is
  version drift, not error, and is reported as such.
- One category, `CHEMONTID:0001218` (Inorganic isocyanides), has **no ChEBI synonym at
  all**, against the paper's "Each ClassyFire category has one or more mapped ChEBI
  terms."

### 3.1 Auditing my own replication — six flaws, two of which changed the verdicts

The first run reported "15 of 26 claims reproduced". That number was inflated. Auditing
the checking logic found:

| # | Flaw | Effect |
|---|---|---|
| F1 | The inorganic count was scored twice — once as a pass (678), once as a fail (674) | Inflated both numerator and denominator |
| F2 | A ratio (1.24 ChEBI synonyms/category) scored as an independent pass while its own numerator had failed | A ratio cannot reproduce when its count does not |
| F3 | Two pass criteria were tolerances I invented (`abs(mean − 5) < 0.5`) | Made a claim pass that does not hold |
| F4 | "Mappings" counted as synonym *lines*, but one line can carry several IDs | MeSH: 966 lines vs 1165 pairs vs 514 IDs — a 20% spread |
| F5 | Sibling checks used inconsistent units (text vs lines vs IDs) | Not comparable |
| F6 | No provenance on the input file; the root's own synonym silently excluded | The whole phase rested on an unverified file |

Plus F7 (near-vacuous claims like ">4800 categories" scored equal to structural ones)
and F8 (latent: depth traversal assumed single-parenthood without asserting it).

**Corrected result: 8 of 11 structural claims, 0 of 5 bookkeeping, 3 of 3 trivial.**
Claims are now tiered so a shape check cannot pad the score, and a claim that holds only
under rounding is reported as CONSISTENT, never as a reproduction.

The audit is preserved in `docs/phase1a_audit.md`.

---

## 4. The test suite, written before it was convenient

`scripts/test_checks.py`, three tiers:

- **L1 unit** — hand-built toy taxonomies with known answers, no data required. Every
  test is a regression for a bug that actually occurred.
- **L2 fixture** — the real ChemOnt against facts taken from the paper.
- **L3 sample** — the streaming scan cross-checked against an **independent naive
  re-implementation** of the same rules over 20,000 rows.

Writing it immediately found **two bugs in the tests themselves** (a null direct parent
returns `None`/unclassifiable, not `False`; and one toy case violated a second R5
sub-rule I had forgotten to expect). The checks were right; the tests were wrong.

Much later, an audit found two tests that **could not fail**: `ok(..., True)` asserting
a literal, and a `>=` comparison between two different quantities. Both replaced.

---

## 5. Phase 1d — the evaluation claims, which I had wrongly written off

The ledger originally marked E1–E9 "not reproducible — depends on expert panels and the
closed codebase". **That was wrong.** Downloading the paper's supplementary files showed
Additional file 4 contains the actual 800-compound test set *and* the seven experts'
markings of every error — 22,043 rows.

| Claim | Paper | Observed | |
|---|---|---|---|
| E1 test compounds | 800 | 800 | exact |
| E5 false positives | 17 | 17 | exact |
| E5 false negatives | 13 | 13 | exact |
| E3 assignments/compound | 26.38 | 26.38 | exact |
| E7 **maximum score** | 7067.24 | **7066.93** | 0.004% residual |
| E9 (six ChEBI metrics) | ~31, ~27, ~45, ~33, ~14, ~94%, 43.6% | 30.85, 26.90, 45.40, 32.90, 14.30, 94.5%, 43.5% | all match |
| E3 total assignments | 21,102 | 21,101 | −1 |
| E4 distinct categories | 1,308 | 1,305 | −3 |
| E6 compounds per category | 2.6 | **16.17** | irreconcilable |

**E7's scoring scheme was recovered from prose.** The paper gives a score of
7067.04/7067.24 but never states the formula. Taking its description literally —
"each category was assigned a normalized weight based on its number of occurrences among
the 800 chemical entities", so that errors on populated categories are penalised more —
the weight `w(c) = occ(c)/800` reproduces the stated maximum as `Σ occ(c)²/800 = 7066.93`,
a residual of 0.004%. The penalty half was **not** recovered, and is reported as
unrecovered rather than fudged.

**Two further errors in the paper:**
- **E6 is wrong.** "Each category was assigned to an average of 2.6 compounds" is
  unreproducible and contradicted by the paper's own figures: 21,102 ÷ 1,308 = 16.13.
  The observed median is 2 and the mode 1 (473 of 1,305 categories occur exactly once),
  so 2.6 looks like a typical-value statistic mislabelled as a mean.
- **E7's percentage drops a digit.** 7067.04 / 7067.24 = 99.997%, printed as "99.97%".
- A manuscript typo: the text says CID 46936568 was missed as a "*purine* nucleotide
  sugar"; its own supplement says "*Pyrimidine*". For a cytidine derivative the
  supplement is right.

**A filtering bug of mine along the way:** the first run reported 0 false negatives.
False-negative rows have no `ChemOntID` by definition — they are *missing* assignments —
and my filter required one. Fixed, giving the exact 13.

Result: **14 of 22 evaluation claims reproduced.**

---

## 6. Phase 1b/1c — one pass over 73,105,281 rows

### 6.1 The first design, and why it failed

The obvious approach: stream the 2 GB zstd file and apply the paper's rules to each row.
Ancestor tests were done with bit masks.

It ran for **109 CPU-minutes without finishing**. Rather than wait, throughput was
measured directly: raw parsing alone runs at **0.64 M rows/s** (1.9 minutes for the whole
file). So essentially 100% of the time was rule evaluation, not I/O.

Attempting to cache per-row results made it worse: distinct alternative-parent tuples
were projected to reach **~31 million keys**, several GB of dictionary.

### 6.2 The redesign that worked

Move the rules out of the hot loop. Accumulate **bounded sufficient statistics** and
evaluate over distinct keys afterwards:

| Statistic | Distinct keys | vs 73.1M rows |
|---|---|---|
| label paths | 3,631 | 20,000× fewer |
| intermediate-node sets | 5,178 | 14,000× fewer |
| (direct parent, alternative parent) pairs | 1,047,341 | 70× fewer |

**Result: 73,105,281 rows in 6 min 35 s, 184,087 rows/s** — and the dumped statistics
mean no re-analysis ever needs another pass over the 2 GB file. Change a rule definition,
rerun the post-pass, done.

### 6.3 Phase 1b — provenance (not paper reproduction)

An important distinction the report must keep: the Zenodo file is a **2026 third-party,
InChIKey-deduplicated aggregation**, not the 77M-compound population the 2016 paper
describes. Verifying it checks the *release notes*, not the paper.

Six of seven figures reproduce exactly: rows, ZINC IDs, both-IDs, unresolved slots,
classes seen, classes total. One does not: **PubChem CID 68,693,295 observed vs
68,693,298 stated, −3** — reproducible, and matching an earlier independent scan.

### 6.4 Phase 1c — testing the paper's own invariants at full scale

The paper states PATH, R4, R5 and R6 as always true. They are not:

| Invariant | Violating rows | Share |
|---|---:|---:|
| PATH — label is a root-path of its direct parent | 198,959 | 0.272% |
| — with a complete, non-null path (the real finding) | **79,171** | 0.108% |
| R5 — alternative parents unrelated to the direct parent | 596,367 * | 0.816% |
| R6 — intermediate nodes | 85,036 | 0.116% |

\* an upper bound on distinct rows; measured over a 200,000-row sample, 1,661 pair
instances against 1,657 distinct rows (+0.24%).

Three things make these robust rather than artefacts:

1. **Taxonomy drift does not explain them.** Only **73** of the 198,959 path violations
   touch any of the 4 edges where the release and ChemOnt 2.1 disagree.
2. **R5 fails under both of the paper's own conflicting definitions.** The main text says
   alternative parents have no *ancestor–descendant* relation to the direct parent;
   Table S1 says the weaker *parent-child*. **69 of the 93** violating pairs are direct
   parent-child, so they breach both. The dominant case, `Aryl thioethers` →
   `Alkylarylthioethers` (551,743 rows), also breaks R4: category reduction should have
   promoted the more specific child.
3. **It contradicts the release's own headline.** Zenodo states "the canonical five-tier
   hierarchical path is now correct for every row"; 79,171 rows carry a complete,
   non-null path that is still not a root-path.

---

## 7. The visualization — five attempts, four of them wrong

This was the longest and least linear part of the work. It is recorded in full because
the failures were informative and the measurements are reusable.

### Attempt 1 — a published artifact
4,824 class nodes, force-directed, canvas + d3. Worked, and is still published. But it
did not show compounds.

### Attempt 2 — 73.1M points, first try
A Vite + React + WebGL2 app streaming 292 MB of int16 coordinates. Then a check of the
actual bytes on disk:

```
first 4 points: [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
px range: nan..nan
```

**The force layout had diverged to NaN**, collapsing all 73.1M points onto one pixel.
The compound counts were perfectly correct, so nothing else would have caught it — it
would have shipped as a single dot. Replaced with damped Fruchterman–Reingold plus three
assertions that now run every build: layout is finite, layout has not collapsed, and a
re-read of the written file confirms <2% of points at the origin.

### Attempt 3 — circle "islands", and the measurement that killed the whole approach
Grouping classes into 31 superclass islands looked reasonable. Measuring the rendered
density did not:

```
occupied pixels: 30,703 of 810,000 (3.8%)
counts/pixel: median 1,494   p99 14,606   max 28,084
77% of occupied pixels exceed 255
```

**96% of the canvas empty, and every occupied pixel saturated to flat white.** No
arrangement fixes this: 73.1M points on ~2M pixels is a density problem, not a layout
problem.

### Attempt 4 — squarified treemap
Rectangles tile without gaps, raising computed disc coverage to 51%. But measured
occupancy stayed at 3.0%. Investigation found **a dead radius block from the previous
layout still running afterwards and overwriting the treemap's radii** — 155 classes with
>100k compounds pinned at the 0.0016 floor. A real bug, found only because the computed
and measured numbers disagreed and the discrepancy was chased instead of ignored.

### Attempt 5 — the one that worked: radial dendrogram + hierarchical edge bundling
A force layout was tried on the 4,824 classes, rendered to PNG, and **looked at**. It was
a textbook hairball: the de-overlap pass had flung every node to the rim and the tree's
hub structure filled the middle with a white starburst.

That picture gave the insight: **this data is a tree** (4,822 edges for 4,824 nodes), and
force layouts are the wrong tool for trees. A radial tree has an exact layout with no
edge crossings — radius = depth, angle = subtree extent — and the confusability links can
be **bundled along the hierarchy** (Holten), which is the standard technique for exactly
this shape of data.

Three further bugs were caught by rendering and inspecting each iteration:

- **int16 clipping.** The strip ran to ±2.89 world units while the encoder used a fixed
  ÷32000 divisor built for ±1. `x range -32767..32767` meant the horizontal spread was
  being flattened back into a band. Scale now derived from the true extent, with an
  overflow guard.
- **Islands not centred with their contents.** Halo centres reached `cy=+1.56` while
  nodes only reached ±0.91 — the halos were scaled but never shifted. Same transform now
  applied to both.
- **Solid white fan wedges.** A parent with 259 children stacked 259 near-parallel spokes
  into a block. Link opacity is now scaled by sibling count so a dense fan sums to the
  same ink as a sparse one.

**The lesson, stated plainly:** four attempts were shipped or nearly shipped on reasoning
alone. Every one of them was fixed only after rendering the result and looking at it, or
measuring the actual pixels. Visual work cannot be validated by reasoning about the code.

---

## 8. Frontend correctness — bugs found by measurement, not inspection

| Bug | How it was found | Impact |
|---|---|---|
| Point culling used `nr * 4.2` instead of the generator's real halo `hr` | Compared the two formulas over all classes | **1,031 of 3,617 classes (28.5%)** had halos up to **19.6× larger** than the cull box; their compounds vanished when the node centre left the viewport |
| `stream()` yielded once per network chunk | Arithmetic on chunk count × clamped `setTimeout` | **~18 s of pure delay** and 4,462 React re-renders at 64 KB chunks; 71 s at 16 KB |
| Label collision boxes ignored rotation | Read the code against the rendering | Rejected labels that fit; permitted overlaps at diagonal angles |
| `measureText` used weight 500 while drawing 600 | Code read | Reserved box too short for the emphasised label |
| Labels culled by node **centre**, nodes by **radius** | Measured node radius at high zoom | At 99k zoom a disc is **1,020 px** in radius; its centre can be 900 px off-screen while it fills the view. Labels were skipped for plainly visible nodes |
| Label offset scaled with radius | Same measurement | At 99k zoom the label sat **1,029 px** from the node — off-screen |
| Every node is "big" past 50k zoom | Counted nodes over 30 px by zoom | 4,824 of 4,824 at 50k — so the "big node" clamp piled every off-screen label onto the viewport edge, colliding again |
| Zoom readout frozen | Read the code after a user report | `view` is a `useRef`; mutating it never re-renders. The map *was* zooming; the number was captured once at mount |
| **`schedule` identity churn re-streamed 292 MB** | Found during the SOLID refactor | `paint` was keyed on `labels`/`mode`, so toggling either changed `schedule`, re-ran the init effect, and **aborted the load and re-streamed the whole cloud from zero** |

The last one is the worst bug in the project and would have fired on every single mode
switch.

---

## 9. Parallel audits

Three workers were run concurrently on non-overlapping file scopes so they could not
collide: a SOLID refactor of `frontend/src/`, a deployment setup creating only new files,
and a read-only audit of the Python pipeline.

### The Python audit — six confirmed defects, two of which changed published numbers

- **`check_r6` compared depths where the paper requires ancestry.** Both the Results
  section and Table S1 say intermediate nodes must be *ascendants of the direct parent*;
  the code tested `depth[n] <= depth[dp]`, silently accepting **7,970 rows** carrying a
  node at the right depth on the wrong lineage. **Correct: 85,036, not 84,333.**
- **The 79,164 figure mixed two populations.** It was derived as 198,959 − 119,795, but
  only 119,788 of those are path violations: the 7 rows with a null direct parent are
  unresolved yet never reach `check_path`. **Correct: 79,171**, now emitted directly by
  the scan rather than subtracted in prose.
- **The ancestor bitmasks were shifted by the raw node id**, whose maximum is the root's
  9,999,999 — so every mask was a **10-million-bit integer**: 302 µs per R5 check and
  **6.4 GB** of int objects. *This was the real cause of the 109-minute scan* that had
  been attributed to "rule evaluation" in general. Replaced with set intersections.
- Two tests that could not fail (§4), a docstring claim the code does not implement,
  an int16 margin down to 0.35%, and a double normalisation.

The full 73.1M scan was re-run after the fixes; both corrected figures came back exactly
as predicted.

### The deployment audit
Found that `test_checks.py` and `05_build_graph_data.py` read from a sibling directory
**outside the repo**. In CI that directory is absent, so L2/L3 **skipped silently and the
suite still exited 0** — a green build testing 19 of 31 assertions. Paths now resolve via
`scripts/paths.py` (`$ATLAS_DATA` → in-repo `data/` → dev layout), and a skipped tier
exits 2.

---

## 10. Shipping

**Repository.** 28 commits (at the time of writing), ordered by how the project builds
up rather than chronologically — taxonomy model → Phase 1a → docs → tests → the scan → evaluation →
visualization data → frontend → refactor → fixes → deploy → docs. 76 tracked files,
kept small because `.gitignore` (verified with `git check-ignore`) excludes the 292 MB point cloud,
the 118 MB supplementary CSV, the 12 MB pairs file and all generated assets.

**Deployment.** Built locally (the server is aarch64; so is the build machine), and only
three things shipped: the compiled app (**176 KB**), `layout.json` (3.5 MB) and
`points.bin` (292 MB). No sources, no `node_modules`, no dataset, no git history.

Both data files were verified by **sha256 on both ends**. That mattered: macOS ships
`openrsync`, which reported "sent 88 bytes… total size is 735" while transferring
*nothing*. The checksum caught it; transfer switched to a plain pipe.

**Caching.** The cloud is immutable, so its URL carries a content hash
(`points.bin?v=50d4044cf9a1`) and is served `max-age=31536000, immutable` — a returning
visitor never re-fetches 292 MB. `layout.json` is deliberately *not* immutable, since it
carries that version; it revalidates by ETag, so a repeat visit is a 304. Range requests
stay enabled so an interrupted 292 MB stream resumes rather than restarting, and gzip is
off for the binary (it compresses only 14%) and on for the JSON (3.67 MB → 973 KB).

**Isolation.** The site runs as its own container on its own port with read-only mounts,
sharing no network, volume or compose project with the pre-existing service on the host.

**Production hygiene.** A final audit removed a dev-server URL shown to phone users,
a HUD debug readout (raw zoom, internal label counts), raw exception text in the error
banner, stray mode numbering, and dead CSS.

---

## 11. What is honestly *not* reproduced

Stated plainly, because the strongest version of this project is the one that says so:

- **We have never run ClassyFire.** The paper is fundamentally about a classification
  program — four algorithmic steps, >9000 SMARTS patterns, Markush structures, ~200 IUPAC
  regexes. Not a single molecule was classified here. What was reproduced is the taxonomy
  *artefact* and the *bookkeeping* of the evaluation.
- **R1–R3** (the rule base) — never released. Impossible for anyone outside the group.
- **E2** (249.9 s on their hardware), **S1–S6, S8** (77M compounds, the 424-hour PubChem
  run, 550 ms/compound, 13,100 descriptions) — need their servers.
- **The R5 upper bound.** 596,367 counts pair instances, not distinct rows; the true
  distinct-row count is roughly 1,400 lower.
- **R5's alt-vs-alt clause is never evaluated** by the scan — the sufficient statistics
  do not retain per-row alternative-parent sets. Measured over 200,000 rows, zero rows
  violate that clause without also violating dp-vs-alt, so no published figure depends
  on it; only the completeness claim was wrong, and it has been corrected.

**Scorecard: ~25 of ~41 attemptable claims reproduced, 16 measured and mismatched with
explanations, ~11 structurally impossible.**

---

## 12. Lessons that generalise

1. **Measure before optimising.** 109 CPU-minutes were spent on a scan whose real
   bottleneck (10-million-bit integers) was found only by a later audit. One throughput
   measurement at the start would have exposed it.
2. **Audit your own verification.** The first replication score was inflated by six
   distinct flaws in the checking logic. A verification tool needs verifying.
3. **Look at visual output.** Four visualization designs were wrong in ways no amount of
   code review would have revealed. Rendering a PNG and opening it found each one.
4. **Trust checksums, not tools' self-reports.** `openrsync` claimed success while
   transferring nothing.
5. **A skipped test is not a passing test.** A green CI covering 19 of 31 assertions is
   worse than a red one.
6. **Instrument for the question you cannot answer.** Adding "named / no room /
   off-screen" counters turned "why is this label missing?" from speculation into a
   reading.
7. **Distinguish "the paper is wrong" from "our reproduction is wrong".** Both occurred
   here. E6, the 4825 arithmetic and the 678/674 contradiction are the paper's; the R6
   depth test and the 79,164 subtraction were ours. Only careful measurement separates them.
