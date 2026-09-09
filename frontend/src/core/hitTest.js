/** Nearest node whose (padded) disc contains the screen point, or null. */
export function pickNode(nodes, view, mx, my) {
  let best = null, bd = Infinity;
  for (const n of nodes) {
    const dx = n.px * view.zoom + view.pan[0] - mx, dy = n.py * view.zoom + view.pan[1] - my;
    const rr = Math.max(7, n.nr * view.zoom + 5), dd = dx * dx + dy * dy;
    if (dd < rr * rr && dd < bd) { bd = dd; best = n; }
  }
  return best;
}
