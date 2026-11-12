// Writes SAMPLE public/phase2.json and public/bipartite.json for UI development.
// Every value is fabricated (but shaped like the real outputs, with real class ids so
// the edges land on the map) and the files carry meta.fixture = true. A file that
// exists without that flag is real output and is never touched.
//
//   node scripts/make-fixtures.mjs
import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const pub = join(dirname(fileURLToPath(import.meta.url)), "..", "public");
const layout = JSON.parse(readFileSync(join(pub, "layout.json"), "utf8"));

let seed = 798;
const rnd = () => {                         // mulberry32, so reruns are identical
  seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
  let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
};
const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(6.2832 * rnd());
const r = (v, d = 3) => Number(v.toFixed(d));
const today = new Date().toISOString().slice(0, 10);

function writeIfFixture(name, obj) {
  const p = join(pub, name);
  if (existsSync(p)) {
    let cur = null;
    try { cur = JSON.parse(readFileSync(p, "utf8")); } catch { /* unreadable: leave it */ }
    if (!cur?.meta?.fixture) { console.log(`skip ${name}: real output present`); return; }
  }
  writeFileSync(p, JSON.stringify(obj));
  console.log(`wrote ${name} (${(JSON.stringify(obj).length / 1024).toFixed(0)} KB)`);
}

const nodes = layout.nodes;
const byId = new Map(nodes.map((n) => [n.i, n]));
const dpNodes = nodes.filter((n) => n.dp > 0 && n.e > 0);
const scNames = [...new Set(dpNodes.map((n) => n.sc).filter(Boolean))];

function communities(members, count) {
  const comm = new Map();
  for (const n of members) {
    let c = Math.abs(scNames.indexOf(n.sc) * 7 + (n.d % 2)) % count;
    if (rnd() < 0.27) c = Math.floor(rnd() * count);
    comm.set(n.i, c);
  }
  const out = [];
  for (let c = 0; c < count; c++) {
    const ms = members.filter((n) => comm.get(n.i) === c);
    if (!ms.length) continue;
    const tally = new Map();
    ms.forEach((n) => tally.set(n.sc || "kingdom-level", (tally.get(n.sc || "kingdom-level") || 0) + 1));
    const top = [...tally.entries()].sort((a, b) => b[1] - a[1]);
    out.push({ id: c, size: ms.length,
      label: top.slice(0, 2).map(([s]) => s).join(" / "), purity: r(top[0][1] / ms.length) });
  }
  return { comm, list: out.sort((a, b) => b.size - a.size) };
}

// ---------------------------------------------------------------- phase2.json
{
  const { comm, list } = communities(dpNodes, 15);
  const flagged = new Set(dpNodes.filter((n) => n.e > 4 && n.dp > 200)
    .sort(() => rnd() - 0.5).slice(0, 20).map((n) => n.i));
  const out = {};
  const deg = new Map(), str = new Map();
  const edges = [];
  for (const l of layout.links) {
    if (edges.length >= 2600) break;
    if (!byId.get(l.s)?.dp || !byId.get(l.t)?.dp || rnd() < 0.35) continue;
    const pmi = r(0.4 + Math.log10(l.w) * 0.55 + gauss() * 0.4);
    const roll = rnd();
    edges.push([l.s, l.t, pmi, l.w, roll < 0.55 ? 1 : roll < 0.96 ? 2 : 3]);
    for (const id of [l.s, l.t]) {
      deg.set(id, (deg.get(id) || 0) + 1); str.set(id, (str.get(id) || 0) + pmi);
    }
  }
  for (const n of dpNodes) {
    const H = Math.max(0, n.e * 0.62 + gauss() * 0.35 + (flagged.has(n.i) ? 0.9 : 0));
    out[n.i] = { H: r(H), Hu: r(H + n.e * 0.3 + rnd() * 0.4), comm: comm.get(n.i),
      str: r(str.get(n.i) || 0, 2), deg: deg.get(n.i) || 0, flag: flagged.has(n.i) ? 1 : 0 };
  }
  const hist = Array(30).fill(0);
  for (let i = 0; i < 10000; i++) {
    const v = 0.5 + gauss() * 0.075;
    hist[Math.max(0, Math.min(29, Math.floor((v - 0.2) / 0.02)))]++;
  }
  const sweep = [];
  for (const theta of [0.5, 1, 2]) for (const tau of [0.01, 0.05]) for (const nmin of [50, 200]) {
    const e = Math.round(42000 / theta * (tau === 0.01 ? 0.6 : 1) * (nmin === 50 ? 1.3 : 1));
    sweep.push({ theta, tau, nmin, edges: e, jaccard: r(0.62 + rnd() * 0.3),
      p_h1a: r(0.004 + rnd() * 0.06, 4) });
  }
  writeIfFixture("phase2.json", {
    meta: { generated: today, fixture: true, primary: "criterion A ∩ B, θ=1, τ=0.05, n≥50",
      nodes: dpNodes.length, edges_all: 1045955, edges_competing: 41873 },
    headline: [
      { label: "competing edges", value: "41,873", note: "of 1,045,955 dp–alt pairs" },
      { label: "H1a p (size-matched)", value: "0.012", note: "n = 20 expert-flagged classes" },
      { label: "modularity Q", value: "0.47", note: "null 0.21 ± 0.01" },
      { label: "NMI vs SuperClass", value: "0.58", note: "< 0.7 target" },
    ],
    hypotheses: [
      { id: "H1a", title: "Flagged classes are more ambiguous", statistic: "mean H, flagged − matched",
        observed: 0.71, null: 0.50, p: 0.012, verdict: "supported",
        detail: "20 expert-flagged classes vs 10,000 size-stratified draws." },
      { id: "H1b", title: "Structure beyond the degree sequence", statistic: "weighted clustering (Onnela)",
        observed: 0.31, null: 0.12, p: 0.001, verdict: "supported",
        detail: "z = 19 against 500 degree-preserving rewirings." },
      { id: "H1c", title: "Communities recut the taxonomy", statistic: "NMI(Louvain, SuperClass)",
        observed: 0.58, null: 0.7, p: null, verdict: "inconclusive",
        detail: "Below the 0.7 target but sensitive to the resolution parameter." },
    ],
    edges, nodes: out, communities: list, sweep,
    errorProp: { levels: ["kingdom", "superclass", "class", "subclass", "none"],
      overlay: [0.04, 0.21, 0.38, 0.27, 0.10], uniform: [0.02, 0.09, 0.22, 0.31, 0.36],
      size: [0.03, 0.14, 0.29, 0.30, 0.24] },
    charts: {
      lca_depth: [1, 2, 3, 4, 5, 6, 7, 8].map((depth) => ({ depth,
        pairs: Math.round(4000 * Math.exp(-((depth - 4) ** 2) / 3)),
        mass: Math.round(9e6 * Math.exp(-((depth - 3.3) ** 2) / 2.2)) })),
      entropy_vs_size: dpNodes.filter((n, i) => flagged.has(n.i) || i % 2 === 0)
        .map((n) => [r(Math.log10(n.dp), 3), out[n.i].H, out[n.i].flag]),
      perm_null: { bins: Array.from({ length: 31 }, (_, i) => r(0.2 + i * 0.02)), counts: hist, observed: 0.71 },
      null_compare: [
        { measure: "weighted clustering", observed: 0.31, null_mean: 0.12, null_sd: 0.01, z: 19.0 },
        { measure: "modularity Q", observed: 0.47, null_mean: 0.21, null_sd: 0.012, z: 21.7 },
        { measure: "assortativity (depth)", observed: 0.38, null_mean: 0.02, null_sd: 0.03, z: 12.0 },
        { measure: "assortativity (size)", observed: -0.11, null_mean: -0.09, null_sd: 0.02, z: -1.0 },
      ],
    },
  });
}

// ------------------------------------------------------------- bipartite.json
{
  const members = nodes.filter((n) => n.ip > 0 && n.d >= 2);
  const { comm, list } = communities(members, 12);
  const ids = new Set(members.map((n) => n.i));
  const edges = [];
  const kdeg = new Map();
  for (let i = 1500; i < layout.links.length && edges.length < 2200; i++) {
    const l = layout.links[i];
    if (!ids.has(l.s) || !ids.has(l.t) || rnd() < 0.55) continue;
    edges.push([l.s, l.t, l.w, r(3 + Math.abs(gauss()) * 14, 2)]);
    kdeg.set(l.s, (kdeg.get(l.s) || 0) + 1); kdeg.set(l.t, (kdeg.get(l.t) || 0) + 1);
  }
  const outNodes = {};
  for (const n of members)
    outNodes[n.i] = { k: kdeg.get(n.i) || 0, comm: comm.get(n.i),
      fit: r(Math.log10(n.ip + 1) / 7 + rnd() * 0.15) };
  const logBin = (vals) => {
    const m = new Map();
    for (const v of vals) { if (v <= 0) continue; const b = Math.round(Math.pow(10, Math.round(Math.log10(v) * 5) / 5)); m.set(b, (m.get(b) || 0) + 1); }
    return [...m.entries()].sort((a, b) => a[0] - b[0]);
  };
  const scanM = [0, 2285, 6043, 14323, 17866, 25892, 31721, 38474, 47296, 55253, 64591, 73281, 81210,
    81191, 78143, 70322, 57764, 47288, 35638, 26442, 19685, 14150, 10033, 7247, 4931, 3739, 2843,
    1883, 1405, 936, 734, 444, 281, 204, 127, 102, 63, 43, 35, 31, 12, 12, 8, 5, 6, 2, 5];
  const groups = ["aromatic homomonocyclic compound", "carboxylic acid derivative", "organic oxide",
    "hydrocarbon derivative", "azacycle", "organopnictogen compound", "carbonyl group",
    "amine", "ether", "heteroaromatic compound"];
  const top = [...members].sort((a, b) => outNodes[b.i].fit - outNodes[a.i].fit).slice(0, 8);
  writeIfFixture("bipartite.json", {
    meta: { generated: today, fixture: true, rows: 924007, classes: 3338, substituents: 362, incidences: 11646523 },
    headline: [
      { label: "incidences", value: "11.6M", note: "compound–class labels (DP ∪ alt)" },
      { label: "class projection", value: "451,027", note: "edges before validation" },
      { label: "BiCM backbone", value: "18,940", note: "edges at FDR 0.01" },
      { label: "NODF", value: "0.41", note: "null 0.29 ± 0.01" },
    ],
    modes: [
      { name: "compounds", n: 924007, mean_degree: 12.6,
        degree_hist: scanM.map((c, k) => [k, c]).filter(([k, c]) => k > 0 && c > 0) },
      { name: "classes", n: 3338, mean_degree: 3489.2, degree_hist: logBin(members.map((n) => n.ip)) },
    ],
    projections: [
      { name: "compound × compound", mode: "compound", nodes: 924007, edges: null, density: null,
        clustering: null, components: null, note: "not materialised: ~10¹¹ pairs" },
      { name: "class × class", mode: "class", nodes: 3338, edges: 451027, density: 0.081,
        clustering: 0.74, components: 3, note: "co-labelled on ≥ 1 compound" },
      { name: "substituent × substituent", mode: "substituent", nodes: 362, edges: 14622,
        density: 0.224, clustering: 0.81, components: 1, note: "" },
    ],
    backbones: [
      { method: "raw projection", edges: 451027, density: 0.081, clustering: 0.74, clustering_null: 0.71,
        modularity: 0.18, communities: 9, nmi_superclass: 0.31, ari_superclass: 0.12, nmi_class: 0.27 },
      { method: "hypergeometric", edges: 38211, density: 0.0069, clustering: 0.42, clustering_null: 0.19,
        modularity: 0.51, communities: 22, nmi_superclass: 0.55, ari_superclass: 0.33, nmi_class: 0.49 },
      { method: "BiCM", edges: 18940, density: 0.0034, clustering: 0.37, clustering_null: 0.08,
        modularity: 0.63, communities: 31, nmi_superclass: 0.61, ari_superclass: 0.41, nmi_class: 0.57 },
      { method: "disparity filter", edges: 24502, density: 0.0044, clustering: 0.22, clustering_null: 0.11,
        modularity: 0.44, communities: 18, nmi_superclass: 0.47, ari_superclass: 0.26, nmi_class: 0.40 },
    ],
    edges, nodes: outNodes, communities: list,
    fg: { nodf: 0.41, nodf_null_mean: 0.29, nodf_null_sd: 0.011, mantel_r: 0.34, mantel_p: 0.001,
      top_classes: top.map((n) => ({ id: n.i, name: n.n, fit: outNodes[n.i].fit })),
      top_groups: groups.map((name, i) => ({ id: i, name, complexity: r(1.8 - i * 0.12), ubiquity: 40 + i * 23 })) },
    charts: {
      density_by_level: [["kingdom", 1, 1], ["superclass", 412, 0.89], ["class", 9211, 0.31],
        ["subclass", 51204, 0.12], ["level 5", 141882, 0.071], ["level 6+", 247317, 0.052]]
        .map(([level, e, density]) => ({ level, edges: e, density })),
      hub_removal: [0, 5, 10, 20, 40, 80, 160, 320].map((removed) => ({ removed,
        density: r(0.081 * Math.exp(-removed / 90) + 0.004, 4) })),
      signature_sizes: [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 778]
        .map((s) => [s, Math.max(1, Math.round(338273 * Math.pow(s, -2.1)))]),
    },
  });
}
