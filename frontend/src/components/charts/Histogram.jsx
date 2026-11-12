import React from "react";
import { fmt, num } from "../../core/format.js";
import { Axes, Empty, INK, PAD, W, finite, linear, niceMax, niceTicks } from "./axes.jsx";

/**
 * Contiguous histogram with an optional observed-value marker. `bins` may be the
 * counts.length + 1 edges or the counts.length centres; both are accepted.
 */
export function Histogram({ bins = [], counts = [], observed, h = 170, xLabel, yLabel = "draws" }) {
  const c = Array.isArray(counts) ? counts.map((v) => (finite(v) ? v : 0)) : [];
  const b = Array.isArray(bins) ? bins.filter(finite) : [];
  if (!c.length || b.length < 2) return <Empty />;
  let edges;
  if (b.length === c.length + 1) edges = b;
  else if (b.length === c.length) {
    const step = b[1] - b[0];
    edges = [...b.map((v) => v - step / 2), b[b.length - 1] + step / 2];
  } else return <Empty what="histogram bins do not match counts" />;
  let lo = edges[0], hi = edges[edges.length - 1];
  if (finite(observed)) { lo = Math.min(lo, observed); hi = Math.max(hi, observed); }
  const pad = (hi - lo) * 0.03;
  const x = linear(lo - pad, hi + pad, PAD.l + 2, W - PAD.r - 2);
  const top = niceMax(Math.max(...c));
  const y = linear(0, top, h - PAD.b, PAD.t + 10);
  return (
    <figure className="x-chart">
      <svg viewBox={`0 0 ${W} ${h}`} role="img">
        <Axes w={W} h={h} x={x} y={y} xTicks={niceTicks(lo, hi, 5)} yTicks={niceTicks(0, top, 3)}
              xFmt={(v) => num(v, 2)} xLabel={xLabel} yLabel={yLabel} />
        {c.map((v, i) => (
          <rect key={i} x={x(edges[i]) + 0.4} width={Math.max(0.5, x(edges[i + 1]) - x(edges[i]) - 0.8)}
                y={y(v)} height={y(0) - y(v)} fill={INK.blue} fillOpacity="0.7">
            <title>{`${num(edges[i])} – ${num(edges[i + 1])}: ${fmt(v)}`}</title>
          </rect>))}
        {finite(observed) && (
          <g>
            <line x1={x(observed)} x2={x(observed)} y1={PAD.t + 4} y2={h - PAD.b}
                  stroke={INK.org} strokeWidth="1.8" />
            <text x={x(observed)} y={PAD.t + 4} dy="-0.1em" textAnchor="middle"
                  className="mark">observed {num(observed)}</text>
          </g>)}
      </svg>
    </figure>
  );
}
