import hierarchy from "./hierarchy.js";

/** The taxonomy skeleton, dimmed so the extension edges read on top of it. Only the
 *  selected class's own tree edges are lit; its extension partners are not tree kin. */
export default {
  draw(ctx, frame) {
    ctx.globalAlpha = 0.5;
    hierarchy.draw(ctx, { ...frame, lit: frame.hv ? new Set([frame.hv.i]) : null });
    ctx.globalAlpha = 1;
  },
};
