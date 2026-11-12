import React from "react";
import { num } from "../../core/format.js";
import { sideFor } from "../../core/extension.js";

const NEIGHBOURS = 5;

function community(side, id) {
  if (id == null) return "—";
  const c = side.communities.get(id);
  return c?.label ? `#${id} · ${c.label}` : `#${id}`;
}

/** Extension rows for the right-click card: Phase 2 scores, projection scores, top partners. */
export function ExtensionDetails({ n, ix, layer }) {
  const p2 = ix.p2, bp = ix.bp;
  const p = p2?.nodes.get(n.i), b = bp?.nodes.get(n.i);
  const nbSide = sideFor(ix, layer) || p2 || bp;
  const nb = (nbSide?.adj.get(n.i) || []).slice(0, NEIGHBOURS);
  return (
    <>
      {p2 && (
        <div className="sec">
          <div className="lab">phase 2 · confusability</div>
          {p ? (
            <div className="grid">
              <div><div className="lab">H filtered</div><div className="val">{num(p.H)}<em> bits</em></div></div>
              <div><div className="lab">H unfiltered</div><div className="val">{num(p.Hu)}<em> bits</em></div></div>
              <div><div className="lab">strength · degree</div><div className="val">{num(p.str)}<em> · {num(p.deg)}</em></div></div>
              <div><div className="lab">expert-flagged</div><div className="val">{p.flag === 1 ? "yes" : "no"}</div></div>
              <div className="wide"><div className="lab">community</div><div className="val sm">{community(p2, p.comm)}</div></div>
            </div>
          ) : <div className="val sm">not in the overlay</div>}
        </div>)}
      {bp && (
        <div className="sec">
          <div className="lab">bipartite projection</div>
          {b ? (
            <div className="grid">
              <div><div className="lab">backbone degree k</div><div className="val">{num(b.k)}</div></div>
              <div><div className="lab">fitness</div><div className="val">{num(b.fit)}</div></div>
              <div className="wide"><div className="lab">community</div><div className="val sm">{community(bp, b.comm)}</div></div>
            </div>
          ) : <div className="val sm">not in the projection</div>}
        </div>)}
      {nb.length > 0 && (
        <div className="sec">
          <div className="lab">{nbSide.kind === "p2" ? "top overlay partners · PMI" : "top backbone partners · score"}</div>
          {nb.map((x) => (
            <div className="cf" key={x.id}><span>{ix.d.byId.get(x.id)?.n ?? x.id}</span><b>{num(x.w)}</b></div>))}
        </div>)}
    </>
  );
}
