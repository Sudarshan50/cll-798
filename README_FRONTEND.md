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

## Honest limits

- Position **within** a class is a deterministic Vogel spiral, not a structural embedding
  — we have no 2D coordinates for molecular structure. Position **between** classes is
  meaningful: it is the force-directed layout of the taxonomy plus its confusability
  overlay.
- `points.bin` is a build artefact, not source. Do not commit it.
