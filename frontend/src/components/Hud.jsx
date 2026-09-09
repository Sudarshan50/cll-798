import React from "react";
import { fmt } from "../core/format.js";

export function Hud({ mode, onScreen, loaded, meta, bundles, status }) {
  const pct = meta ? (loaded / meta.points) * 100 : 0;
  return (
    <div className="hud">
      {mode === "compounds" ? (
        <>
          <div className="count">{fmt(onScreen)}</div>
          <div className="of">
            {onScreen ? "compounds drawn at this zoom" : "zoom in to reveal compounds"}
          </div>
          <div className="bar"><i style={{ width: `${pct}%` }} /></div>
          <div className="stat">
            {loaded < (meta?.points ?? 1)
              ? `${fmt(loaded)} of ${fmt(meta?.points)} streamed to GPU`
              : `all ${fmt(meta?.points)} compounds resident`}
          </div>
        </>
      ) : (
        <>
          <div className="count">{fmt(meta?.points)}</div>
          <div className="of">compounds, carried by node size</div>
          <div className="stat">
            {fmt(meta?.classes)} ChemOnt classes · {fmt(bundles)} bundled
            confusability links
          </div>
        </>
      )}
      {status !== "ready" && <div className="stat dim">{status}…</div>}
    </div>
  );
}
