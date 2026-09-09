import React, { useEffect, useRef, useState } from "react";
import { fmt } from "../core/format.js";
import { lineage } from "../core/graph.js";

export function InfoPanel({ n, x, y, d, onClose }) {
  const ref = useRef(null);
  const [pos, setPos] = useState({ left: x + 14, top: y - 60 });
  useEffect(() => {
    const el = ref.current; if (!el) return;
    setPos({ left: Math.min(Math.max(12, x + 14), innerWidth - el.offsetWidth - 14),
             top: Math.min(Math.max(12, y - el.offsetHeight / 2),
                           innerHeight - el.offsetHeight - 14) });
  }, [x, y, n]);
  useEffect(() => {
    const k = (e) => { if (e.key === "Escape") onClose(); };
    addEventListener("keydown", k);
    return () => removeEventListener("keydown", k);
  }, [onClose]);
  const lin = lineage(d, n).map((c) => c.n).reverse();
  const conf = (d.cnb.get(n.i) || []).slice(0, 3);
  return (
    <div className="card pinned" ref={ref} style={pos}>
      <button className="x" onClick={onClose} aria-label="Close">×</button>
      <div className="nm">{n.n}</div>
      <div className="id">{n.c} · level {n.d}</div>
      <div className="path">{lin.slice(0, -1).join(" › ") || "kingdom"}</div>
      <div className="grid">
        <div><div className="lab">compounds here</div><div className="val">{fmt(n.dp)}</div></div>
        <div><div className="lab">with descendants</div><div className="val">{fmt(n.st)}</div></div>
        <div><div className="lab">alt-parent entropy</div><div className="val">{n.e.toFixed(2)}<em> bits</em></div></div>
        <div><div className="lab">superclass</div><div className="val sm">{n.sc || "—"}</div></div>
      </div>
      {conf.length > 0 && (
        <div className="sec">{conf.map(([id, w]) => (
          <div className="cf" key={id}><span>{d.byId.get(id)?.n}</span><b>{fmt(w)}</b></div>))}
        </div>)}
      {n.x?.length > 0 && (
        <div className="sec">{n.x.slice(0, 2).map((e, i) => (
          <div key={i}>
            <div className="ex">{e.k}{e.c ? <em> CID {e.c}</em> : null}</div>
            <div className="ex"><em>{e.s}</em></div>
          </div>))}
        </div>)}
    </div>
  );
}
