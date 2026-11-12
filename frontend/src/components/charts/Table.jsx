import React from "react";
import { num, text } from "../../core/format.js";
import { Empty } from "./axes.jsx";

/**
 * Compact numeric table. `columns` = [{ key, label, fmt?, align?, title? }]; values
 * that are missing render as an em dash. `maxHeight` makes the body scroll.
 */
export function Table({ columns = [], rows = [], maxHeight }) {
  const rs = Array.isArray(rows) ? rows.filter((r) => r && typeof r === "object") : [];
  if (!rs.length) return <Empty />;
  return (
    <div className="x-table" style={maxHeight ? { maxHeight } : undefined}>
      <table>
        <thead>
          <tr>{columns.map((c) => (
            <th key={c.key} title={c.title} className={c.align === "left" ? "l" : ""}>{c.label}</th>))}
          </tr>
        </thead>
        <tbody>
          {rs.map((r, i) => (
            <tr key={i}>{columns.map((c) => {
              const v = r[c.key];
              const s = c.fmt ? c.fmt(v, r) : typeof v === "number" ? num(v) : text(v);
              return <td key={c.key} className={c.align === "left" ? "l" : ""}>{s}</td>;
            })}</tr>))}
        </tbody>
      </table>
    </div>
  );
}
