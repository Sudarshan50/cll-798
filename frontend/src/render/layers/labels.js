import { layoutLabels } from "../../core/labels.js";

const font = (weight) => `${weight} 11px "IBM Plex Sans",system-ui,sans-serif`;

/** Rotated labels on the ring, de-conflicted by angle. */
export default {
  enabled: (f) => f.showLabels && !f.lod,
  draw(ctx, { d, view, hv, link }) {
    ctx.textBaseline = "middle";
    // measure with the weight that is actually drawn, or the collision boxes lie
    const measure = (text, weight) => { ctx.font = font(weight); return ctx.measureText(text).width; };
    // hovering: name the lineage, the children and the confusability partners only
    const { placements, named, noRoom, offScreen } = layoutLabels({
      view,
      list: hv ? [hv, ...link.anc, ...link.rest] : d.byCount,
      hardCap: hv ? Infinity : 46,
      selected: hv,
      measure,
    });

    for (const { n, lx, ly, a, big, weight } of placements) {
      const fill = hv
        ? (n.i === hv.i ? "#FFF3E0"
           : (link.anc.includes(n) ? "rgba(255,214,160,.95)" : "rgba(232,237,247,.92)"))
        : "rgba(232,237,247,.9)";
      ctx.save();
      ctx.translate(lx, ly);
      if (!big) {
        const deg = (a * 180) / Math.PI;
        const flip = deg > 90 || deg < -90;
        ctx.rotate(flip ? a + Math.PI : a);
        ctx.textAlign = flip ? "right" : "left";
      } else {
        ctx.textAlign = "center";
      }
      ctx.font = font(weight);
      ctx.lineWidth = 3.4; ctx.strokeStyle = "rgba(6,9,16,.95)";
      ctx.strokeText(n.n, 0, 0);
      ctx.fillStyle = fill;
      ctx.fillText(n.n, 0, 0);
      ctx.restore();
    }
    return { named, noRoom, offScreen };
  },
};
