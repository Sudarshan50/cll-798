import React from "react";
import { arr } from "../../core/extension.js";
import { num, text } from "../../core/format.js";

export function SampleBadge() {
  return <span className="x-badge" title="Fabricated fixture values for UI development — not results">sample data</span>;
}

export function Section({ title, note, children }) {
  return (
    <section className="x-sec">
      <h3>{title}</h3>
      {note && <p className="note">{note}</p>}
      {children}
    </section>
  );
}

const STATUS_TEXT = {
  idle: "not loaded yet",
  loading: "loading…",
  missing: "not computed — the file is not in public/ yet",
  error: "could not be read (see console)",
};

/** Calls `children(data)` only once the file is ready; otherwise says why not. */
export function FileGate({ file, name, children }) {
  if (file?.status === "ready" && file.data) return children(file.data);
  return <div className="x-empty">{name}: {STATUS_TEXT[file?.status] || STATUS_TEXT.idle}</div>;
}

export function Headline({ items }) {
  const xs = arr(items).filter((h) => h && (h.label || h.value != null));
  if (!xs.length) return null;
  return (
    <div className="x-cards">
      {xs.map((h, i) => (
        <div className="x-card" key={i}>
          <div className="lab">{text(h.label)}</div>
          <div className="v">{typeof h.value === "number" ? num(h.value) : text(h.value)}</div>
          {h.note != null && <div className="nt">{text(h.note)}</div>}
        </div>))}
    </div>
  );
}

export function Meta({ meta, fields }) {
  const m = meta && typeof meta === "object" ? meta : {};
  const parts = fields.map(([k, label, f]) => (m[k] == null ? null : `${label} ${f ? f(m[k]) : text(m[k])}`))
    .filter(Boolean);
  return parts.length ? <div className="x-meta">{parts.join(" · ")}</div> : null;
}
