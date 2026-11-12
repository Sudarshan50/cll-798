import { toScreenX, toScreenY } from "../../core/camera.js";
import { ALPHA_BUCKETS, sideFor } from "../../core/extension.js";
import { EDGE_RGB, HILITE_RGB } from "../../core/extColors.js";

// While panning/zooming only the strongest batches are drawn, up to this many edges.
const LOD_MAX = 900;
// Control points are pulled toward the layout centre so edges arc like the bundles.
const BEND = 0.62;

/** Append one curve to the current path; false if it is entirely off screen. */
function curve(ctx, view, a, b) {
  const ax = toScreenX(view, a.px), ay = toScreenY(view, a.py);
  const bx = toScreenX(view, b.px), by = toScreenY(view, b.py);
  const cx = toScreenX(view, (a.px + b.px) * 0.5 * BEND);
  const cy = toScreenY(view, (a.py + b.py) * 0.5 * BEND);
  if (Math.max(ax, bx, cx) < 0 || Math.min(ax, bx, cx) > view.w ||
      Math.max(ay, by, cy) < 0 || Math.min(ay, by, cy) > view.h) return false;
  ctx.moveTo(ax, ay);
  ctx.quadraticCurveTo(cx, cy, bx, by);
  return true;
}

/** The selected extension edge layer: one stroke per (colour, alpha bucket). */
export default {
  enabled: (f) => !!sideFor(f.ext?.idx, f.ext?.layer),
  draw(ctx, { view, hv, lod, ext }) {
    const side = sideFor(ext.idx, ext.layer);
    const dim = hv ? 0.3 : 1;
    let budget = lod ? LOD_MAX : Infinity, drawn = 0;
    // strongest first under LOD, so the budget keeps what matters
    const order = lod ? [...side.batches].reverse() : side.batches;
    for (const g of order) {
      if (budget <= 0) break;
      const u = (g.bucket + 0.5) / ALPHA_BUCKETS;
      const rgb = EDGE_RGB[side.kind][g.c] || EDGE_RGB[side.kind][Object.keys(EDGE_RGB[side.kind])[0]];
      ctx.strokeStyle = `rgba(${rgb},${((0.05 + 0.42 * Math.pow(u, 1.4)) * dim).toFixed(3)})`;
      ctx.lineWidth = 0.55 + 0.6 * u;
      ctx.beginPath();
      const list = g.edges.length > budget ? g.edges.slice(-budget) : g.edges;
      for (const e of list) if (curve(ctx, view, e.a, e.b)) drawn++;
      budget -= list.length;
      ctx.stroke();
    }
    if (hv) {
      const byId = ext.idx.d.byId;
      for (const x of side.adj.get(hv.i) || []) {
        const m = byId.get(x.id); if (!m) continue;
        ctx.strokeStyle = `rgba(${HILITE_RGB},${(0.35 + 0.55 * x.t).toFixed(3)})`;
        ctx.lineWidth = 0.9 + 0.9 * x.t;
        ctx.beginPath();
        if (curve(ctx, view, hv, m)) drawn++;
        ctx.stroke();
      }
    }
    return { extDrawn: drawn };
  },
};
