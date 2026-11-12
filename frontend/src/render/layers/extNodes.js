import { toScreenX, toScreenY } from "../../core/camera.js";
import { nodeColours } from "../../core/extension.js";
import { FLAG_RING, NO_DATA } from "../../core/extColors.js";

/** Class discs recoloured by the chosen extension scheme; expert-flagged classes ringed. */
export default {
  draw(ctx, { d, view, hv, lit, ext }) {
    const { w, h, zoom } = view;
    const cols = ext.color === "tax" ? null : nodeColours(ext.idx, ext.color);
    const p2 = ext.idx.p2?.nodes;
    for (const n of d.nodes) {
      const x = toScreenX(view, n.px), y = toScreenY(view, n.py);
      const r = Math.max(1.0, n.nr * zoom);
      if (x + r + 4 < 0 || y + r + 4 < 0 || x - r - 4 > w || y - r - 4 > h) continue;
      ctx.globalAlpha = hv ? (n.i === hv.i ? 1 : (lit.has(n.i) ? 0.92 : 0.26)) : 1;
      ctx.beginPath(); ctx.arc(x, y, r, 0, 6.2832);
      ctx.fillStyle = cols ? (cols.get(n.i) || NO_DATA) : n.col; ctx.fill();
      if (hv && n.i === hv.i) {
        ctx.lineWidth = 2.2; ctx.strokeStyle = "#FFF3E0"; ctx.stroke();
      }
      if (p2?.get(n.i)?.flag === 1) {
        ctx.beginPath(); ctx.arc(x, y, Math.max(r + 2.2, 4), 0, 6.2832);
        ctx.lineWidth = 1.5; ctx.strokeStyle = FLAG_RING; ctx.stroke();
      }
    }
    ctx.globalAlpha = 1;
  },
};
