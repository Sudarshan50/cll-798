// Compounds fade in from here. Chosen together with K_HALO in 08_build_network_points.py
// so that on-screen density at this zoom is ~30 points/pixel for EVERY class.
export const POINT_ZOOM = 1500;

/** Upload nothing; just pick the visible class ranges and let the GPU draw them. */
export function renderPoints(cloud, d, view) {
  cloud.clear();
  const { pan, zoom, dpr, w, h } = view;
  if (!cloud.count || zoom < POINT_ZOOM) return 0;
  const fade = Math.min(1, (zoom - POINT_ZOOM) / (POINT_ZOOM * 1.4));
  const ranges = [];
  for (const n of d.nodes) {
    if (!n.cnt) continue;
    const x = n.px * zoom + pan[0], y = n.py * zoom + pan[1];
    // cull against the halo the generator actually used; nr*4.2 was up to 19.6x too
    // small, so 28.5% of classes were dropped while still visibly on screen
    const r = (n.hr ?? n.nr * 1.05) * zoom;
    if (x + r < 0 || y + r < 0 || x - r > w || y - r > h) continue;
    ranges.push([n.off, n.cnt]);
  }
  return cloud.drawRanges(ranges, {
    pan, zoom, dpr, scale: d.meta.scale,
    size: Math.max(1, Math.min(2.5, zoom / 26000)),
    alpha: 0.055 * fade, color: [1.0, 0.76, 0.45],
  });
}
