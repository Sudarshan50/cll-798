# Chemical Space Atlas — full 73.1M point cloud (Vite + React + WebGL2)

```bash
cd frontend
npm install
npm run dev            # http://localhost:5173
```

## What it renders

**All 73,105,281 compounds.** Not a sample. Each point is one real compound, placed in
the region of the ChemOnt class it was actually assigned to; per-class counts come from
the full 73.1M-row scan and the generator asserts the total reconciles exactly.

Two layers share one camera:

| layer | what | how |
|---|---|---|
| WebGL2 point cloud | 73,105,281 compounds | one `POINTS` draw call, additive blending so overlap reads as density |
| 2D overlay | 4,824 classes, 4,822 taxonomy edges, 14,000 confusability edges | canvas, hover + labels |

## Why it fits in 8 GB

The cloud is 292 MB of vertex data (two int16 per point). Holding that twice would not
fit, so `src/gl.js` pre-allocates the GPU buffer once and then **streams** `points.bin`
through a `ReadableStream`, uploading each chunk with `bufferSubData` and dropping it.
Peak CPU memory is one chunk, not 292 MB. The counter in the corner is the true number
of points resident on the GPU, so if the machine cannot take all of them you see exactly
where it stopped rather than a silent truncation.

## Regenerating the data

```bash
python3 scripts/05_build_graph_data.py     # class graph, from the scan statistics
python3 scripts/06_build_pointcloud.py     # layout + frontend/public/points.bin (292 MB)
```

Neither script rescans the 2 GB dataset; both read the sufficient statistics written by
`scripts/02_scan_dataset.py`.

## Extension view

The third mode, **Extension**, presents the two extension analyses on the same layout.
It reads `public/phase2.json` (Phase 2 confusability overlay) and `public/bipartite.json`
(compound–class bipartite graph and its projections), fetched only the first time the mode
is opened and cached after that; the map and the `points.bin` stream do not wait on them.

- **Map.** The right-hand panel's `edges` selector draws the competing overlay, the
  BiCM-validated class projection, or nothing in place of the confusability bundles; `colour`
  switches between taxonomy, either community partition, and a ramp of ambiguity `H`.
  Amber rings mark the expert-flagged classes. Left click lights a class's extension edges;
  right click adds its Phase 2 and projection scores and top partners to the details card.
- **Panel tabs.** *Phase 2 · Confusability* (hypotheses, LCA depth, H vs size, permutation
  null, null comparison, error propagation, sweep), *Bipartite & projections* (schematic,
  degree distributions, projections, density by level, hub removal, backbones, functional
  groups) and *Method* (plain-language summary and caveats).
- **Missing data.** Every key in either file is optional; anything absent shows as
  "not computed". A file with `meta.fixture: true` is fabricated sample data and the UI
  says so with a *sample data* badge.

`node scripts/make-fixtures.mjs` writes such fixtures for development. It never replaces
a file that lacks `meta.fixture`, so dropping the real outputs into `public/` is enough.

Code: `hooks/useExtensionData.js` (loading), `core/extension.js` + `core/extColors.js`
(indexing, colour schemes), `render/layers/ext*.js` (map layers, used only in this mode),
`components/extension/*` (panel, legend, card section), `components/charts/*` (SVG charts).

## Honest limits

- Position **within** a class is a deterministic Vogel spiral, not a structural embedding
  — we have no 2D coordinates for molecular structure. Position **between** classes is
  meaningful: it is the force-directed layout of the taxonomy plus its confusability
  overlay.
- `points.bin` is a build artefact, not source. Do not commit it.
