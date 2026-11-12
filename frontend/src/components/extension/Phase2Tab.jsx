import React from "react";
import { fmt, num, pval, text } from "../../core/format.js";
import { arr, obj } from "../../core/extension.js";
import { Bars } from "../charts/Bars.jsx";
import { Scatter } from "../charts/Scatter.jsx";
import { Histogram } from "../charts/Histogram.jsx";
import { Table } from "../charts/Table.jsx";
import { Empty, INK, finite } from "../charts/axes.jsx";
import { FileGate, Headline, Meta, Section } from "./shared.jsx";

const VERDICT = { supported: "supported", "not supported": "not", inconclusive: "inconclusive" };
const pct = (v) => `${Math.round(v * 100)}%`;

function Hypotheses({ list }) {
  const hs = arr(list).filter((h) => h && typeof h === "object");
  if (!hs.length) return <Empty />;
  return hs.map((h, i) => (
    <div className="x-hyp" key={h.id ?? i}>
      <div className="top">
        <b>{text(h.id)}</b>
        <span className="t">{h.title == null ? "" : text(h.title)}</span>
        <span className={`x-verdict ${VERDICT[h.verdict] || "inconclusive"}`}>{h.verdict == null ? "not computed" : text(h.verdict)}</span>
      </div>
      <div className="stat">
        {h.statistic == null ? "statistic" : text(h.statistic)}: <b>{num(h.observed)}</b>
        {finite(h.null) && <> vs null {num(h.null)}</>}
        <span className="p">p {pval(h.p)}</span>
      </div>
      {typeof h.null === "string" && <div className="dt">null: {h.null}</div>}
      {h.detail != null && <div className="dt">{text(h.detail)}</div>}
    </div>));
}

function LcaDepth({ rows }) {
  const rs = arr(rows).filter((r) => r && (finite(r.depth) || typeof r.depth === "string"));
  if (!rs.length) return <Empty />;
  const sp = rs.reduce((s, r) => s + (finite(r.pairs) ? r.pairs : 0), 0) || 1;
  const sm = rs.reduce((s, r) => s + (finite(r.mass) ? r.mass : 0), 0) || 1;
  return (
    <Bars categories={rs.map((r) => r.depth)} xLabel="depth of lowest common ancestor"
      yLabel="share" yFmt={pct}
      series={[
        { name: "pairs", color: INK.blue, values: rs.map((r) => (finite(r.pairs) ? r.pairs / sp : null)), raw: rs.map((r) => r.pairs) },
        { name: "co-assignment mass", color: INK.org, values: rs.map((r) => (finite(r.mass) ? r.mass / sm : null)), raw: rs.map((r) => r.mass) },
      ]}
      tip={(s, c, v) => `${s.name} at LCA depth ${c}: ${pct(v)} (${fmt(s.raw[rs.findIndex((r) => r.depth === c)])})`} />
  );
}

function EntropySize({ rows }) {
  const rs = arr(rows).filter((r) => Array.isArray(r) && finite(r[0]) && finite(r[1]));
  const pts = [
    ...rs.filter((r) => r[2] !== 1).map((r) => ({ x: r[0], y: r[1], color: INK.ink2, opacity: 0.32, r: 1.3,
      title: `10^${num(r[0])} compounds · H ${num(r[1])}` })),
    ...rs.filter((r) => r[2] === 1).map((r) => ({ x: r[0], y: r[1], color: INK.org, opacity: 1, r: 3.2,
      title: `expert-flagged · 10^${num(r[0])} compounds · H ${num(r[1])}` })),
  ];
  return <Scatter points={pts} xLabel="log10 class size (direct-parent compounds)" yLabel="ambiguity H (bits)" />;
}

function ErrorProp({ ep }) {
  const e = obj(ep);
  const levels = arr(e.levels);
  const series = [["overlay", INK.org, "overlay-guided"], ["uniform", INK.mute, "uniform"], ["size", INK.blue, "size-weighted"]]
    .filter(([k]) => Array.isArray(e[k])).map(([k, color, name]) => ({ name, color, values: e[k] }));
  return <Bars categories={levels} series={series} xLabel="level at which the error lands" yLabel="share" yFmt={pct}
               tip={(s, c, v) => `${s.name} · ${c}: ${pct(v)}`} />;
}

export function Phase2Tab({ file }) {
  return (
    <FileGate file={file} name="phase2.json">
      {(p) => {
        const ch = obj(p.charts), pn = obj(ch.perm_null);
        return (
          <>
            <Meta meta={p.meta} fields={[["primary", "primary:"], ["nodes", "nodes", fmt],
              ["edges_competing", "competing edges", fmt], ["edges_all", "all pairs", fmt], ["generated", "built"]]} />
            <Headline items={p.headline} />
            <Section title="Hypotheses" note="Pre-registered in docs/phase2_plan.md §1; p from the stated null.">
              <Hypotheses list={p.hypotheses} />
            </Section>
            <Section title="Where co-assignments sit in the tree"
                     note="Depth of the lowest common ancestor of each competing pair: count of pairs vs compound mass.">
              <LcaDepth rows={ch.lca_depth} />
            </Section>
            <Section title="Ambiguity is not just size" note="Per-class H against class size; expert-flagged classes in amber.">
              <EntropySize rows={ch.entropy_vs_size} />
            </Section>
            <Section title="H1a permutation null" note="Size-stratified draws of 20 classes; amber line is the flagged set.">
              <Histogram bins={pn.bins} counts={pn.counts} observed={pn.observed} xLabel="mean H of a 20-class draw" />
            </Section>
            <Section title="Observed vs degree-preserving null">
              <Table rows={ch.null_compare} columns={[
                { key: "measure", label: "measure", align: "left" },
                { key: "observed", label: "observed" },
                { key: "null_mean", label: "null", fmt: (v, r) => (finite(v) ? `${num(v)} ± ${num(r.null_sd, 2)}` : "—") },
                { key: "z", label: "z", fmt: (v) => num(v, 3) },
              ]} />
            </Section>
            <Section title="Error propagation" note="Where a misassignment lands when guided by the overlay vs uniform / size-weighted baselines.">
              <ErrorProp ep={p.errorProp} />
            </Section>
            <Section title="Robustness sweep" note="Filter thresholds vs edge count, Jaccard overlap with the primary graph, and H1a p.">
              <Table rows={p.sweep} maxHeight={176} columns={[
                { key: "theta", label: "θ" }, { key: "tau", label: "τ" }, { key: "nmin", label: "n min" },
                { key: "edges", label: "edges" }, { key: "jaccard", label: "Jaccard" },
                { key: "p_h1a", label: "p H1a", fmt: pval },
              ]} />
            </Section>
            <Section title="Overlay communities" note="Louvain on the competing overlay; purity = share of the dominant SuperClass.">
              <Table rows={arr(p.communities).slice().sort((a, b) => (b?.size ?? 0) - (a?.size ?? 0))} maxHeight={176}
                columns={[{ key: "id", label: "#" }, { key: "size", label: "size" },
                  { key: "purity", label: "purity", fmt: (v) => (finite(v) ? pct(v) : "—") },
                  { key: "label", label: "dominant superclasses", align: "left" }]} />
            </Section>
          </>
        );
      }}
    </FileGate>
  );
}
