import React, { useState } from "react";
import { COLOURINGS, EDGE_LAYERS } from "../../core/extension.js";
import { SampleBadge } from "./shared.jsx";
import { Phase2Tab } from "./Phase2Tab.jsx";
import { BipartiteTab } from "./BipartiteTab.jsx";
import { MethodTab } from "./MethodTab.jsx";

/** Horizontal room the open panel takes from the right edge, for floating cards to avoid. */
export const PANEL_SPACE = 404 + 24 + 10;

const TABS = [
  { id: "p2", label: "Phase 2 · Confusability" },
  { id: "bp", label: "Bipartite & projections" },
  { id: "method", label: "Method" },
];

/**
 * Right-hand panel for the Extension mode. The map selectors stay visible when it is
 * collapsed; the tabs and charts only render while it is open.
 */
export function ExtensionPanel({ open, onToggle, files, layer, onLayer, colour, onColour }) {
  const [tab, setTab] = useState("p2");
  const fixture = files.phase2.data?.meta?.fixture || files.bipartite.data?.meta?.fixture;
  return (
    <aside className={`xpanel${open ? "" : " shut"}`} aria-label="Extension analyses">
      <header>
        <h2>Extension</h2>
        {fixture && <SampleBadge />}
        <span className="grow" />
        <button className="tog" onClick={onToggle} aria-expanded={open}>
          {open ? "collapse" : "expand"}
        </button>
      </header>
      <div className="x-map">
        <label htmlFor="x-layer">edges</label>
        <select id="x-layer" autoComplete="off" value={layer} onChange={(e) => onLayer(e.target.value)}>
          {EDGE_LAYERS.map((o) => <option key={o.id} value={o.id}>{o.label}</option>)}
        </select>
        <label htmlFor="x-colour">colour</label>
        <select id="x-colour" autoComplete="off" value={colour} onChange={(e) => onColour(e.target.value)}>
          {COLOURINGS.map((o) => <option key={o.id} value={o.id}>{o.label}</option>)}
        </select>
      </div>
      {open && (
        <>
          <div className="x-tabs" role="tablist">
            {TABS.map((t) => (
              <button key={t.id} role="tab" aria-selected={tab === t.id}
                className={tab === t.id ? "on" : ""} onClick={() => setTab(t.id)}>{t.label}</button>))}
          </div>
          <div className="x-body" role="tabpanel" key={tab}>
            {tab === "p2" && <Phase2Tab file={files.phase2} />}
            {tab === "bp" && <BipartiteTab file={files.bipartite} />}
            {tab === "method" && <MethodTab />}
          </div>
        </>
      )}
    </aside>
  );
}
