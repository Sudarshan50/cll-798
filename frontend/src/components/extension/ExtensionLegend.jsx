import React from "react";
import { fmt, num, text } from "../../core/format.js";
import { COLOURINGS, EDGE_LAYERS, sideFor } from "../../core/extension.js";
import { EDGE_LABELS, EDGE_RGB, FLAG_RING, H_GRADIENT, communityColour } from "../../core/extColors.js";

const SHOW_COMMUNITIES = 6;

function Swatch({ color, ring }) {
  return <i className={ring ? "ring" : ""} style={ring ? { borderColor: color } : { background: color }} />;
}

function ColourKey({ ix, colour }) {
  if (colour === "tax")
    return (
      <>
        <div className="lg"><Swatch color="#F2A65A" />organic superclasses (warm hues)</div>
        <div className="lg"><Swatch color="#4FC3B0" />inorganic superclasses (cool hues)</div>
      </>);
  if (colour === "H") {
    if (!ix?.H) return <div className="lg dim">ambiguity H not computed</div>;
    return (
      <div className="lg ramp">
        <span>{num(ix.H.lo, 2)}</span><b style={{ background: H_GRADIENT }} /><span>{num(ix.H.hi, 2)} bits</span>
      </div>);
  }
  const side = colour === "p2c" ? ix?.p2 : ix?.bp;
  if (!side || !side.communities.size) return <div className="lg dim">communities not computed</div>;
  const top = [...side.communities.values()].slice(0, SHOW_COMMUNITIES);
  return (
    <>
      {top.map((c, i) => (
        <div className="lg" key={c.id}><Swatch color={communityColour(i)} />
          <span className="nm">{c.label ? text(c.label) : `community ${text(c.id)}`}</span><em>{num(c.size)}</em></div>))}
      {side.communities.size > SHOW_COMMUNITIES &&
        <div className="lg dim">+ {side.communities.size - SHOW_COMMUNITIES} smaller communities</div>}
    </>);
}

/** Hud contents in Extension mode: edge count, edge key, colour key, flag ring. */
export function ExtensionLegend({ ix, layer, colour, files }) {
  const side = sideFor(ix, layer);
  const layerName = EDGE_LAYERS.find((l) => l.id === layer)?.label;
  const file = layer === "p2" ? files.phase2 : layer === "bp" ? files.bipartite : null;
  const fixture = files.phase2.data?.meta?.fixture || files.bipartite.data?.meta?.fixture;
  return (
    <>
      <div className="count">{side ? fmt(side.edges.length) : "—"}</div>
      <div className="of">
        {layer === "none" ? "no extension edges drawn"
          : side ? `${layerName} edges` : `${layerName}: ${file?.status === "loading" ? "loading…" : "not computed"}`}
      </div>
      <div className="legend">
        {side && Object.entries(EDGE_LABELS[side.kind]).filter(([c]) => side.bySource[c]).map(([c, label]) => (
          <div className="lg" key={c}><Swatch color={`rgb(${EDGE_RGB[side.kind][c]})`} />{label}
            <em>{fmt(side.bySource[c])}</em></div>))}
        {side && <div className="lg dim">opacity ∝ weight rank · left click lights a class's edges</div>}
        <div className="lg head">{COLOURINGS.find((c) => c.id === colour)?.label}</div>
        <ColourKey ix={ix} colour={colour} />
        {ix?.p2?.flagged > 0 && (
          <div className="lg"><Swatch color={FLAG_RING} ring />expert-flagged class<em>{fmt(ix.p2.flagged)}</em></div>)}
        {fixture && <div className="lg sample">sample data — fixture values, not results</div>}
      </div>
    </>
  );
}
