import React from "react";
import { fmt, num, pval, text } from "../../core/format.js";
import { arr, obj } from "../../core/extension.js";
import { Bars } from "../charts/Bars.jsx";
import { Scatter } from "../charts/Scatter.jsx";
import { LogLog } from "../charts/LogLog.jsx";
import { Table } from "../charts/Table.jsx";
import { Empty, INK, finite } from "../charts/axes.jsx";
import { FileGate, Headline, Meta, Section } from "./shared.jsx";

// compound -> classes it is labelled with, for the schematic
const LABELS = [[0, 1], [0], [0, 1, 2], [1, 2], [2]];
const CY = [26, 50, 74, 98, 122], KY = [38, 74, 110];
const KP = [[200, 46], [260, 20], [320, 46]];
const PENT = [0, 1, 2, 3, 4].map((i) => {
  const a = (-90 + i * 72) * Math.PI / 180;
  return [260 + 26 * Math.cos(a), 110 + 22 * Math.sin(a)];
});

function Schematic() {
  const shared = (a, b) => LABELS.filter((l) => l.includes(a) && l.includes(b)).length;
  const cc = [];
  for (let i = 0; i < 5; i++) for (let j = i + 1; j < 5; j++)
    if (LABELS[i].some((k) => LABELS[j].includes(k))) cc.push([i, j]);
  return (
    <figure className="x-chart x-schematic">
      <svg viewBox="0 0 360 150" role="img"
           aria-label="Compounds linked to the classes they carry; projecting onto one side links two nodes that share a neighbour">
        <text x="28" y="11" textAnchor="middle" className="title">compounds</text>
        <text x="100" y="11" textAnchor="middle" className="title">classes</text>
        {LABELS.map((ks, c) => ks.map((k) => (
          <line key={`${c}${k}`} x1="28" y1={CY[c]} x2="100" y2={KY[k]} stroke={INK.mute} strokeWidth="1" />)))}
        {CY.map((y, c) => <circle key={c} cx="28" cy={y} r="5" fill={INK.blue}><title>{`compound ${c + 1}`}</title></circle>)}
        {KY.map((y, k) => (
          <g key={k}><rect x="94" y={y - 6} width="12" height="12" rx="2" fill={INK.org} />
            <text x="114" y={y} dy="0.32em">{"ABC"[k]}</text></g>))}
        <path d="M128 74 C 148 74, 150 46, 172 44 M128 74 C 148 74, 150 106, 172 108"
              fill="none" stroke={INK.ink2} strokeWidth="1" markerEnd="url(#x-arrow)" />
        <defs>
          <marker id="x-arrow" viewBox="0 0 6 6" refX="5" refY="3" markerWidth="6" markerHeight="6" orient="auto">
            <path d="M0 0 L6 3 L0 6 z" fill={INK.ink2} />
          </marker>
        </defs>
        {[[0, 1], [1, 2], [0, 2]].map(([a, b]) => (
          <line key={`k${a}${b}`} x1={KP[a][0]} y1={KP[a][1]} x2={KP[b][0]} y2={KP[b][1]}
                stroke={INK.org} strokeOpacity="0.7" strokeWidth={shared(a, b) * 1.4}>
            <title>{`${"ABC"[a]}–${"ABC"[b]}: ${shared(a, b)} shared compound(s)`}</title>
          </line>))}
        {KP.map(([x, y], k) => (
          <g key={k}><rect x={x - 6} y={y - 6} width="12" height="12" rx="2" fill={INK.org} />
            <text x={x} y={y} dy="0.32em" textAnchor="middle" style={{ fill: "#0D1220" }}>{"ABC"[k]}</text></g>))}
        <text x="260" y="68" textAnchor="middle" className="title">class projection (weighted)</text>
        {cc.map(([a, b]) => (
          <line key={`c${a}${b}`} x1={PENT[a][0]} y1={PENT[a][1]} x2={PENT[b][0]} y2={PENT[b][1]}
                stroke={INK.blue} strokeOpacity="0.6" strokeWidth="1" />))}
        {PENT.map(([x, y], c) => <circle key={c} cx={x} cy={y} r="4.5" fill={INK.blue} />)}
        <text x="260" y="146" textAnchor="middle" className="title">compound projection (cliques)</text>
      </svg>
    </figure>
  );
}

function Degrees({ modes }) {
  const ms = arr(modes).filter((m) => m && typeof m === "object");
  if (!ms.length) return <Empty />;
  const colours = [INK.blue, INK.org, INK.inorg, INK.violet];
  return ms.map((m, i) => (
    <div key={m.name ?? i}>
      <div className="x-cap">
        <b>{m.name == null ? `mode ${i + 1}` : text(m.name)}</b> · n {num(m.n)} · mean degree {num(m.mean_degree)}
      </div>
      <LogLog h={160} series={[{ name: m.name == null ? "degree" : text(m.name), color: colours[i % colours.length], points: m.degree_hist }]}
              xLabel="degree k" yLabel="count" />
    </div>));
}

function FgCard({ fg }) {
  const f = obj(fg);
  const z = finite(f.nodf) && finite(f.nodf_null_mean) && f.nodf_null_sd > 0
    ? (f.nodf - f.nodf_null_mean) / f.nodf_null_sd : null;
  const classes = arr(f.top_classes).filter(Boolean), groups = arr(f.top_groups).filter(Boolean);
  if (!Object.keys(f).length) return <Empty />;
  return (
    <div className="x-fg">
      <div className="x-cards">
        <div className="x-card"><div className="lab">NODF nestedness</div>
          <div className="v">{num(f.nodf)}</div>
          <div className="nt">null {num(f.nodf_null_mean)} ± {num(f.nodf_null_sd, 2)}{z != null ? ` · z ${num(z, 3)}` : ""}</div></div>
        <div className="x-card"><div className="lab">Mantel r</div>
          <div className="v">{num(f.mantel_r)}</div>
          <div className="nt">p {pval(f.mantel_p)} · class vs group similarity</div></div>
      </div>
      <div className="x-two">
        <div>
          <div className="lab">top classes · fitness</div>
          {classes.length ? classes.slice(0, 8).map((c, i) => (
            <div className="cf" key={i}><span>{text(c.name ?? c.id)}</span><b>{num(c.fit)}</b></div>))
            : <div className="x-empty">not computed</div>}
        </div>
        <div>
          <div className="lab">top groups · complexity</div>
          {groups.length ? groups.slice(0, 8).map((g, i) => (
            <div className="cf" key={i} title={`ubiquity ${num(g.ubiquity)}`}>
              <span>{text(g.name ?? g.id)}</span><b>{num(g.complexity)}</b></div>))
            : <div className="x-empty">not computed</div>}
        </div>
      </div>
    </div>
  );
}

const pct = (v) => (finite(v) ? `${Math.round(v * 100)}%` : "—");

export function BipartiteTab({ file }) {
  return (
    <>
      <Section title="From a bipartite graph to its projections"
               note="Compounds carry several class labels (direct parent ∪ alternative parents). Projecting onto one side links two nodes that share a neighbour: class edges are weighted by shared compounds, and every class turns its compounds into a clique.">
        <Schematic />
      </Section>
      <FileGate file={file} name="bipartite.json">
        {(b) => {
          const ch = obj(b.charts), dbl = arr(ch.density_by_level).filter(Boolean);
          return (
            <>
              <Meta meta={b.meta} fields={[["rows", "rows", fmt], ["classes", "classes", fmt],
                ["substituents", "substituents", fmt], ["incidences", "incidences", fmt], ["generated", "built"]]} />
              <Headline items={b.headline} />
              <Section title="Degree distributions" note="Both sides of the compound–class graph, log-log.">
                <Degrees modes={b.modes} />
              </Section>
              <Section title="One-mode projections">
                <Table rows={b.projections} columns={[
                  { key: "name", label: "projection", align: "left", fmt: (v, r) => v ?? r.mode ?? "—" },
                  { key: "nodes", label: "nodes" }, { key: "edges", label: "edges" },
                  { key: "density", label: "density" }, { key: "clustering", label: "C" },
                  { key: "components", label: "comp." },
                ]} />
                {arr(b.projections).filter((p) => p?.note).map((p, i) => (
                  <div className="x-cap dim" key={i}>{text(p.name ?? p.mode)}: {text(p.note)}</div>))}
              </Section>
              <Section title="Projection density by taxonomy level" note="Classes truncated to each level, then projected.">
                <Bars categories={dbl.map((r) => r.level)} yLabel="density" yFmt={(v) => num(v, 2)}
                  series={[{ name: "density", color: INK.inorg, values: dbl.map((r) => r.density) }]}
                  tip={(s, c, v) => `${c}: density ${num(v)} · ${fmt(dbl.find((r) => r.level === c)?.edges)} edges`} />
              </Section>
              <Section title="Hub removal" note="Density after deleting the k highest-degree classes: how much of the projection is a few hubs.">
                <Scatter line color={INK.org} xLabel="hubs removed" yLabel="density" yFmt={(v) => num(v, 2)}
                  points={arr(ch.hub_removal).filter(Boolean).map((r) => ({ x: r.removed, y: r.density,
                    title: `${fmt(r.removed)} removed: density ${num(r.density)}` }))} />
              </Section>
              <Section title="Backbones" note="Raw projection vs three filters. C is observed / null clustering; NMI and ARI against the SuperClass partition.">
                <Table rows={b.backbones} columns={[
                  { key: "method", label: "method", align: "left" },
                  { key: "edges", label: "edges" }, { key: "density", label: "density" },
                  { key: "clustering", label: "C obs/null", fmt: (v, r) => `${num(v, 2)} / ${num(r.clustering_null, 2)}` },
                  { key: "modularity", label: "Q", fmt: (v) => num(v, 2) },
                  { key: "communities", label: "comm." },
                  { key: "nmi_superclass", label: "NMI SC", fmt: (v) => num(v, 2) },
                  { key: "ari_superclass", label: "ARI SC", fmt: (v) => num(v, 2) },
                  { key: "nmi_class", label: "NMI cls", fmt: (v) => num(v, 2) },
                ]} />
              </Section>
              <Section title="Functional groups" note="Class × substituent matrix: nestedness against a null, and economic-complexity style fitness / complexity.">
                <FgCard fg={b.fg} />
              </Section>
              <Section title="Label-signature sizes" note="How many compounds share each exact set of class labels.">
                <LogLog h={160} series={[{ name: "signatures", color: INK.violet, points: ch.signature_sizes }]}
                        xLabel="compounds per signature" yLabel="signatures" />
              </Section>
              <Section title="Backbone communities">
                <Table rows={arr(b.communities).slice().sort((x, y) => (y?.size ?? 0) - (x?.size ?? 0))} maxHeight={176}
                  columns={[{ key: "id", label: "#" }, { key: "size", label: "size" },
                    { key: "purity", label: "purity", fmt: pct },
                    { key: "label", label: "dominant superclasses", align: "left" }]} />
              </Section>
            </>
          );
        }}
      </FileGate>
    </>
  );
}
