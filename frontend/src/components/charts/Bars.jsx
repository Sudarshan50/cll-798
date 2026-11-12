import React from "react";
import { num } from "../../core/format.js";
import { Axes, Empty, Legend, PAD, W, finite, linear, niceMax, niceTicks } from "./axes.jsx";

/**
 * Grouped vertical bars. `series` = [{ name, color, values }] aligned with `categories`.
 * `tip(series, category, value)` writes the hover text; missing values draw nothing.
 */
export function Bars({ categories = [], series = [], h = 170, yLabel, xLabel, yFmt,
                       tip = (s, c, v) => `${s.name} · ${c}: ${num(v)}` }) {
  const ss = series.filter((s) => Array.isArray(s.values));
  const vals = ss.flatMap((s) => s.values.filter(finite));
  if (!categories.length || !vals.length) return <Empty />;
  const top = niceMax(Math.max(...vals, 0));
  const y = linear(0, top, h - PAD.b, PAD.t);
  const band = (W - PAD.l - PAD.r) / categories.length;
  const x = (i) => PAD.l + band * (i + 0.5);
  const bw = Math.min(26, (band * 0.78) / ss.length);
  return (
    <figure className="x-chart">
      {ss.length > 1 && <Legend items={ss.map((s) => ({ label: s.name, color: s.color }))} />}
      <svg viewBox={`0 0 ${W} ${h}`} role="img">
        <Axes w={W} h={h} x={x} y={y} yTicks={niceTicks(0, top, 4)} yFmt={yFmt}
              xTicks={categories.map((c, i) => ({ v: i, label: String(c) }))}
              xLabel={xLabel} yLabel={yLabel} />
        {categories.map((c, i) => ss.map((s, j) => {
          const v = s.values[i];
          if (!finite(v)) return null;
          const bx = x(i) - (bw * ss.length) / 2 + j * bw;
          return (
            <rect key={`${i}-${j}`} x={bx + 0.5} width={bw - 1} y={y(Math.max(0, v))}
                  height={Math.max(0, y(0) - y(Math.max(0, v)))} fill={s.color} rx="1.5">
              <title>{tip(s, c, v)}</title>
            </rect>);
        }))}
      </svg>
    </figure>
  );
}
