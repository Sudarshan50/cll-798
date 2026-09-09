import bundles from "./layers/bundles.js";
import hierarchy from "./layers/hierarchy.js";
import nodes from "./layers/nodes.js";
import labels from "./layers/labels.js";

const LAYERS = [bundles, hierarchy, nodes, labels];

/** Draw every enabled layer, back to front, and sum the counters they report. */
export function renderScene(canvas, frame) {
  const { dpr, w, h } = frame.view;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  let stats = { named: 0, noRoom: 0, offScreen: 0 };
  for (const layer of LAYERS) {
    if (layer.enabled && !layer.enabled(frame)) continue;
    const s = layer.draw(ctx, frame);
    if (s) stats = { ...stats, ...s };
  }
  return stats;
}

export function sizeOverlay(ov, w, h, dpr) {
  ov.width = w * dpr; ov.height = h * dpr;
  ov.style.width = w + "px"; ov.style.height = h + "px";
}
