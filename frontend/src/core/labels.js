import { toScreenX, toScreenY } from "./camera.js";

/** How many non-ancestor names we can afford to draw at this zoom.
 *  6 @ 0.5k, ~14 @ 2k, ~39 @ 10k, uncapped past 30k. */
export function labelBudget(zoom) {
  if (zoom > 30000) return Infinity;
  return Math.max(6, Math.round(6 * Math.pow(zoom / 500, 0.62)));
}

function fits(boxes, bx, by, bw, bh) {
  for (const r of boxes)
    if (bx < r[0] + r[2] && bx + bw > r[0] && by < r[1] + r[3] && by + bh > r[1])
      return false;
  boxes.push([bx, by, bw, bh]); return true;
}

/**
 * Decide where each name goes and which ones lose. `measure(text, weight)` must use
 * the same font weight the caller will actually draw with, or the boxes are wrong.
 * Returns screen-space placements plus the counters the HUD reports.
 */
export function layoutLabels({ view, list, hardCap, selected, measure }) {
  const { w, h } = view;
  const boxes = [], taken = [], placements = [];
  let placed = 0, named = 0, noRoom = 0, offScreen = 0;

  for (const n of list) {
    if (placed >= hardCap) break;
    const x = toScreenX(view, n.px), y = toScreenY(view, n.py);
    const r = Math.max(1.0, n.nr * view.zoom);
    // Cull by the node's CIRCLE, matching the node pass. Culling by centre with a
    // fixed 40px margin dropped labels for nodes that were plainly visible: at 99k
    // zoom a disc can be 1020px across with its centre far outside the viewport.
    if (x + r < -40 || y + r < -40 || x - r > w + 40 || y - r > h + 40) {
      offScreen++; continue;   // was an uncounted `continue`, so the readout said
    }                          // "46 linked · 0 named" with no explanation
    const a = Math.atan2(n.py, n.px);
    if (!selected) {
      if (taken.some((t) => Math.abs(a - t) < 0.075 ||
                            Math.abs(Math.abs(a - t) - 6.2832) < 0.075)) continue;
      taken.push(a);
    }
    const weight = selected && n.i === selected.i ? 600 : 500;
    const tw = measure(n.n, weight);

    // A big disc gets its name set upright in the middle and clamped into view.
    // Offsetting by r + 9 put the label 1029px away from a zoomed-in node, i.e.
    // off-screen, which is why highlighted nodes appeared to have no name.
    // Past ~50k zoom EVERY node is >30px, so clamping all of them piled the labels
    // of off-screen-centre nodes onto the viewport edge, where they collided and
    // vanished. Clamp only when the disc actually covers the middle of the view.
    const big = r > 30;
    const centreVisible = x >= 0 && x <= w && y >= 0 && y <= h;
    const dominates = Math.hypot(x - w / 2, y - h / 2) < r;
    let lx, ly, bx, by, bw, bh;
    if (big && (centreVisible || dominates)) {
      lx = centreVisible ? Math.min(Math.max(x, tw / 2 + 10), w - tw / 2 - 10)
                         : Math.min(Math.max(w / 2, tw / 2 + 10), w - tw / 2 - 10);
      ly = centreVisible ? Math.min(Math.max(y, 16), h - 16) : h / 2;
      bx = lx - tw / 2; by = ly - 7; bw = tw; bh = 14;
    } else if (big) {
      offScreen++; continue;              // large, but nowhere sensible to put a name
    } else {
      const deg = (a * 180) / Math.PI;
      const flip = deg > 90 || deg < -90;
      lx = x + Math.cos(a) * (r + 9);
      ly = y + Math.sin(a) * (r + 9);
      // the collision box is the axis-aligned hull of the ROTATED rectangle's corners
      const ca = Math.cos(a), sa = Math.sin(a);
      const dx = flip ? -tw : tw, hh = 7;
      const xs = [0, dx * ca + hh * sa, dx * ca - hh * sa, hh * sa, -hh * sa];
      const ys = [0, dx * sa - hh * ca, dx * sa + hh * ca, -hh * ca, hh * ca];
      bx = lx + Math.min(...xs); by = ly + Math.min(...ys);
      bw = Math.max(...xs) - Math.min(...xs);
      bh = Math.max(...ys) - Math.min(...ys);
    }
    if (bx + bw < 0 || by + bh < 0 || bx > w || by > h) { offScreen++; continue; }
    if (!fits(boxes, bx, by, bw, bh)) { noRoom++; continue; }
    placed++; named++;
    placements.push({ n, lx, ly, a, big, weight });
  }
  return { placements, named, noRoom, offScreen };
}
