import { toScreenX, toScreenY } from "../../core/camera.js";

export default {
  draw(ctx, { d, view, hv, lit }) {
    const { w, h, zoom } = view;
    for (const n of d.nodes) {
      const x = toScreenX(view, n.px), y = toScreenY(view, n.py);
      const r = Math.max(1.0, n.nr * zoom);
      if (x + r < 0 || y + r < 0 || x - r > w || y - r > h) continue;
      ctx.globalAlpha = hv ? (n.i === hv.i ? 1 : (lit.has(n.i) ? 0.92 : 0.30)) : 1;
      ctx.beginPath(); ctx.arc(x, y, r, 0, 6.2832);
      ctx.fillStyle = n.col; ctx.fill();
      if (hv && n.i === hv.i) {
        ctx.lineWidth = 2.2; ctx.strokeStyle = "#FFF3E0"; ctx.stroke();
      }
    }
    ctx.globalAlpha = 1;
  },
};
