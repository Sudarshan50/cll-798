import React from "react";
import { compact } from "../../core/format.js";

export const W = 360;
export const PAD = { l: 42, r: 10, t: 10, b: 32 };
export const INK = { org: "#F2A65A", inorg: "#4FC3B0", blue: "#7498E0", violet: "#B28CE8",
                     mute: "#6B7891", ink2: "#A8B3C9" };

export const finite = (v) => typeof v === "number" && Number.isFinite(v);

export const linear = (d0, d1, r0, r1) => {
  const k = d1 === d0 ? 0 : (r1 - r0) / (d1 - d0);
  return (v) => r0 + (v - d0) * k;
};
export const log10 = (d0, d1, r0, r1) => {
  const s = linear(Math.log10(d0), Math.log10(d1), r0, r1);
  return (v) => s(Math.log10(v));
};

/** Round tick positions covering [lo, hi]. */
export function niceTicks(lo, hi, n = 4) {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / n, p = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw) || raw;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(10));
  return out;
}
export function niceMax(v) {
  if (!(v > 0)) return 1;
  const t = niceTicks(0, v, 4), last = t[t.length - 1];
  return last >= v ? last : +(last + (t.length > 1 ? t[1] - t[0] : v)).toPrecision(10);
}
export function decades(lo, hi) {
  const out = [];
  for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) out.push(Math.pow(10, e));
  return out;
}

export function Empty({ what = "not computed" }) {
  return <div className="x-empty">{what}</div>;
}

export function Legend({ items }) {
  return (
    <div className="x-legend">
      {items.map((it) => (
        <span key={it.label}><i style={{ background: it.color }} />{it.label}</span>))}
    </div>
  );
}

/**
 * Grid, tick labels and axis titles. `xTicks` may carry `{v, label}` objects when the
 * x axis is categorical; otherwise ticks are numbers placed through `x`.
 */
export function Axes({ w, h, x, y, xTicks = [], yTicks = [], xFmt = compact, yFmt = compact,
                       xLabel, yLabel }) {
  const x0 = PAD.l, x1 = w - PAD.r, y0 = PAD.t, y1 = h - PAD.b;
  return (
    <g className="axes">
      {yTicks.map((t) => (
        <g key={`y${t}`}>
          <line x1={x0} x2={x1} y1={y(t)} y2={y(t)} className="grid" />
          <text x={x0 - 5} y={y(t)} dy="0.32em" textAnchor="end">{yFmt(t)}</text>
        </g>))}
      <line x1={x0} x2={x1} y1={y1} y2={y1} className="axis" />
      <line x1={x0} x2={x0} y1={y0} y2={y1} className="axis" />
      {xTicks.map((t, i) => {
        const v = typeof t === "object" ? t.v : t;
        return (
          <text key={`x${i}`} x={x(v)} y={y1 + 12} textAnchor="middle">
            {typeof t === "object" ? t.label : xFmt(t)}
          </text>);
      })}
      {xLabel && <text className="title" x={(x0 + x1) / 2} y={h - 4} textAnchor="middle">{xLabel}</text>}
      {yLabel && (
        <text className="title" transform={`translate(10 ${(y0 + y1) / 2}) rotate(-90)`}
              textAnchor="middle">{yLabel}</text>)}
    </g>
  );
}
