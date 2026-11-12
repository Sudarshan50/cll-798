import React from "react";
import { fmt } from "../../core/format.js";
import { Axes, Empty, Legend, PAD, W, decades, finite, log10 } from "./axes.jsx";

const sup = (v) => {
  const e = Math.round(Math.log10(v));
  return e >= -1 && e <= 3 ? fmt(v) : `1e${e}`;
};

/** Log-log scatter for distributions. `series` = [{ name, color, points: [[x, y], ...] }]. */
export function LogLog({ series = [], h = 190, xLabel, yLabel }) {
  const ss = series.map((s) => ({
    ...s, points: (Array.isArray(s.points) ? s.points : [])
      .filter((p) => Array.isArray(p) && finite(p[0]) && finite(p[1]) && p[0] > 0 && p[1] > 0),
  })).filter((s) => s.points.length);
  if (!ss.length) return <Empty />;
  const all = ss.flatMap((s) => s.points);
  const xs = decades(Math.min(...all.map((p) => p[0])), Math.max(...all.map((p) => p[0])));
  const ys = decades(Math.min(...all.map((p) => p[1])), Math.max(...all.map((p) => p[1])));
  const x = log10(xs[0], xs[xs.length - 1] > xs[0] ? xs[xs.length - 1] : xs[0] * 10, PAD.l + 4, W - PAD.r - 4);
  const y = log10(ys[0], ys[ys.length - 1] > ys[0] ? ys[ys.length - 1] : ys[0] * 10, h - PAD.b - 3, PAD.t + 2);
  const thin = (t) => (t.length > 6 ? t.filter((_, i) => i % 2 === 0) : t);
  return (
    <figure className="x-chart">
      {ss.length > 1 && <Legend items={ss.map((s) => ({ label: s.name, color: s.color }))} />}
      <svg viewBox={`0 0 ${W} ${h}`} role="img">
        <Axes w={W} h={h} x={x} y={y} xTicks={thin(xs)} yTicks={thin(ys)} xFmt={sup} yFmt={sup}
              xLabel={xLabel} yLabel={yLabel} />
        {ss.map((s) => s.points.map((p, i) => (
          <circle key={`${s.name}${i}`} cx={x(p[0])} cy={y(p[1])} r="2.1" fill={s.color} fillOpacity="0.85">
            <title>{`${s.name}: ${fmt(p[0])} → ${fmt(p[1])}`}</title>
          </circle>)))}
      </svg>
    </figure>
  );
}
