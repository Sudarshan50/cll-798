import { toScreenX, toScreenY } from "../../core/camera.js";

export default {
  enabled: (f) => !f.lod,
  draw(ctx, { d, view, hv }) {
    const X = (x) => toScreenX(view, x), Y = (y) => toScreenY(view, y);
    ctx.lineWidth = Math.max(0.4, 0.55);
    for (const b of d.bundles) {
      const c = b.c;
      if (c.length < 3) continue;
      const on = hv && (b.a === hv.i || b.b === hv.i);
      ctx.strokeStyle = on ? "rgba(255,205,140,.75)" : "rgba(116,152,224,.085)";
      ctx.lineWidth = on ? 1.5 : 0.55;
      ctx.beginPath();
      ctx.moveTo(X(c[0][0]), Y(c[0][1]));
      for (let i = 1; i < c.length - 1; i++) {
        const mx = (c[i][0] + c[i + 1][0]) / 2, my = (c[i][1] + c[i + 1][1]) / 2;
        ctx.quadraticCurveTo(X(c[i][0]), Y(c[i][1]), X(mx), Y(my));
      }
      ctx.quadraticCurveTo(X(c[c.length - 1][0]), Y(c[c.length - 1][1]),
                           X(c[c.length - 1][0]), Y(c[c.length - 1][1]));
      ctx.stroke();
    }
  },
};
