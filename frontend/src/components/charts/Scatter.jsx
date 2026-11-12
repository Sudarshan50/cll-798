import React from "react";
import { num } from "../../core/format.js";
import { Axes, Empty, INK, PAD, W, finite, linear, niceTicks } from "./axes.jsx";

/**
 * Linear x/y scatter. `points` = [{ x, y, color?, r?, title? }]. Points are drawn in
 * order, so put the ones that must stay visible last. `line` joins them by x.
 */
export function Scatter({ points = [], h = 190, xLabel, yLabel, xFmt, yFmt, line = false,
                          color = INK.ink2 }) {
  const ps = points.filter((p) => p && finite(p.x) && finite(p.y));
  if (!ps.length) return <Empty />;
  let x0 = Math.min(...ps.map((p) => p.x)), x1 = Math.max(...ps.map((p) => p.x));
  let y0 = Math.min(...ps.map((p) => p.y)), y1 = Math.max(...ps.map((p) => p.y));
  if (line) y0 = Math.min(0, y0);
  if (x1 === x0) { x0 -= 1; x1 += 1; }
  if (y1 === y0) { y0 -= 1; y1 += 1; }
  const xt = niceTicks(x0, x1, 5), yt = niceTicks(y0, y1, 4);
  if (yt.length > 1 && yt[yt.length - 1] < y1) yt.push(+(yt[yt.length - 1] + yt[1] - yt[0]).toPrecision(10));
  x0 = Math.min(x0, xt[0]); y0 = Math.min(y0, yt[0]); y1 = Math.max(y1, yt[yt.length - 1]);
  const x = linear(x0, x1, PAD.l + 4, W - PAD.r - 4);
  const y = linear(y0, y1, h - PAD.b - 4, PAD.t + 2);
  const path = line
    ? [...ps].sort((a, b) => a.x - b.x).map((p, i) => `${i ? "L" : "M"}${x(p.x).toFixed(1)} ${y(p.y).toFixed(1)}`).join("")
    : null;
  return (
    <figure className="x-chart">
      <svg viewBox={`0 0 ${W} ${h}`} role="img">
        <Axes w={W} h={h} x={x} y={y} xTicks={xt} yTicks={yt} xFmt={xFmt} yFmt={yFmt}
              xLabel={xLabel} yLabel={yLabel} />
        {path && <path d={path} fill="none" stroke={color} strokeWidth="1.5" />}
        {ps.map((p, i) => (
          <circle key={i} cx={x(p.x)} cy={y(p.y)} r={p.r ?? (line ? 2.6 : 1.6)}
                  fill={p.color || color} fillOpacity={p.opacity ?? 0.8}>
            <title>{p.title ?? `${num(p.x)}, ${num(p.y)}`}</title>
          </circle>))}
      </svg>
    </figure>
  );
}
