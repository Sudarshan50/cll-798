import { toScreenX, toScreenY } from "../../core/camera.js";

/** The radial dendrogram skeleton; opacity scaled by fan size. */
export default {
  draw(ctx, { d, view, lit }) {
    const X = (x) => toScreenX(view, x), Y = (y) => toScreenY(view, y);
    for (const t of d.treelinks) {
      const p = t.p;
      if (lit && (lit.has(t.pa) || lit.has(t.ch))) {
        ctx.strokeStyle = "rgba(255,214,160,.55)";
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.moveTo(X(p[0][0]), Y(p[0][1]));
        for (let i = 1; i < p.length; i++) ctx.lineTo(X(p[i][0]), Y(p[i][1]));
        ctx.stroke();
        continue;
      }
      const fan = 1 / (1 + t.s / 40);
      // Floor of 0.055: the previous 0.012 floor made 48% of the hierarchy invisible
      // (a parent with 259 children dropped to alpha 0.023).
      ctx.strokeStyle = `rgba(168,180,206,${Math.min(0.40,
        (t.d <= 3 ? 0.34 : 0.24) * fan + 0.055).toFixed(3)})`;
      ctx.lineWidth = t.d <= 3 ? 0.9 : 0.6;
      ctx.beginPath();
      ctx.moveTo(X(p[0][0]), Y(p[0][1]));
      for (let i = 1; i < p.length; i++) ctx.lineTo(X(p[i][0]), Y(p[i][1]));
      ctx.stroke();
    }
  },
};
