import { NO_DATA, communityColour, rampRgb } from "./extColors.js";

export const arr = (v) => (Array.isArray(v) ? v : []);
export const obj = (v) => (v && typeof v === "object" && !Array.isArray(v) ? v : {});
const finite = (v) => typeof v === "number" && Number.isFinite(v);

export const EDGE_LAYERS = [
  { id: "p2", label: "Competing overlay (Phase 2)" },
  { id: "bp", label: "Validated projection (BiCM, lift ≥ 2)" },
  { id: "none", label: "none" },
];
export const COLOURINGS = [
  { id: "tax", label: "Taxonomy" },
  { id: "p2c", label: "Phase 2 communities" },
  { id: "bpc", label: "Projection communities" },
  { id: "H", label: "Ambiguity H" },
];

/** Strokes are batched per (colour, alpha bucket); this many buckets per colour. */
export const ALPHA_BUCKETS = 6;

/**
 * One side of the extension data (phase2.json or bipartite.json), resolved against the
 * layout once: node records keyed by layout id, edges bound to layout nodes, adjacency
 * sorted by weight, and edges pre-grouped into stroke batches.
 */
function indexSide(d, src, kind) {
  if (!src || typeof src !== "object") return null;
  const nodes = new Map();
  for (const [k, v] of Object.entries(obj(src.nodes))) {
    const id = Number(k);
    if (d.byId.has(id) && v && typeof v === "object") nodes.set(id, v);
  }

  const communities = new Map(), commRank = new Map();
  arr(src.communities).filter((c) => c && c.id != null)
    .sort((a, b) => (b.size ?? 0) - (a.size ?? 0))
    .forEach((c, i) => { communities.set(c.id, c); commRank.set(c.id, i); });

  const raw = [];
  for (const e of arr(src.edges)) {
    if (!Array.isArray(e) || e[0] == null || e[1] == null) continue;   // Number(null) is node 0
    const a = d.byId.get(Number(e[0])), b = d.byId.get(Number(e[1]));
    if (!a || !b || a === b) continue;
    // the overlay is weighted by PMI; the backbone by its significance score
    const w = kind === "p2" ? e[2] : (finite(e[3]) ? e[3] : e[2]);
    if (!finite(w)) continue;
    const src3 = kind === "p2" ? ([1, 2, 3].includes(e[4]) ? e[4] : 1) : 0;
    raw.push({ a, b, w, c: src3, t: 0 });
  }
  // rank, not min-max: PMI can be negative and scores are heavy-tailed
  raw.sort((x, y) => x.w - y.w);
  const n = raw.length;
  raw.forEach((e, i) => { e.t = n > 1 ? i / (n - 1) : 1; });

  const adj = new Map();
  const push = (id, other, e) => {
    let l = adj.get(id); if (!l) adj.set(id, (l = []));
    l.push({ id: other, w: e.w, t: e.t, c: e.c });
  };
  const batches = new Map();
  for (const e of raw) {
    push(e.a.i, e.b.i, e); push(e.b.i, e.a.i, e);
    const bucket = Math.min(ALPHA_BUCKETS - 1, Math.floor(e.t * ALPHA_BUCKETS));
    const key = e.c * ALPHA_BUCKETS + bucket;
    let g = batches.get(key); if (!g) batches.set(key, (g = { c: e.c, bucket, edges: [] }));
    g.edges.push(e);
  }
  adj.forEach((l) => l.sort((x, y) => y.w - x.w));

  const bySource = {};
  for (const e of raw) bySource[e.c] = (bySource[e.c] || 0) + 1;
  let flagged = 0;
  nodes.forEach((v) => { if (v.flag === 1) flagged++; });

  return {
    kind, src, nodes, communities, commRank, adj, edges: raw, bySource, flagged,
    batches: [...batches.values()].sort((x, y) => x.bucket - y.bucket),
    fixture: !!src.meta?.fixture,
  };
}

function hRange(p2) {
  const hs = [];
  p2?.nodes.forEach((v) => { if (finite(v.H)) hs.push(v.H); });
  if (!hs.length) return null;
  hs.sort((a, b) => a - b);
  const lo = hs[Math.floor(hs.length * 0.02)], hi = hs[Math.floor(hs.length * 0.98)];
  return { lo, hi: hi > lo ? hi : lo + 1 };
}

/** Index both extension files against the layout. Either may be null. */
export function indexExtension(d, phase2, bipartite) {
  const p2 = indexSide(d, phase2, "p2");
  const bp = indexSide(d, bipartite, "bp");
  return { d, p2, bp, H: hRange(p2), colours: new Map() };
}

/** Map from node id to fill colour for a colouring scheme; built once per scheme. */
export function nodeColours(ix, scheme) {
  if (ix.colours.has(scheme)) return ix.colours.get(scheme);
  let m = null;
  const side = scheme === "p2c" ? ix.p2 : scheme === "bpc" ? ix.bp : null;
  if (side) {
    m = new Map();
    side.nodes.forEach((v, id) => {
      const rank = side.commRank.get(v.comm);
      m.set(id, v.comm == null ? NO_DATA : communityColour(rank ?? Infinity));
    });
  } else if (scheme === "H" && ix.p2 && ix.H) {
    m = new Map();
    const { lo, hi } = ix.H;
    ix.p2.nodes.forEach((v, id) => {
      if (finite(v.H)) m.set(id, rampRgb((v.H - lo) / (hi - lo)));
    });
  }
  ix.colours.set(scheme, m);
  return m;
}

/** Selection neighbourhood on an extension layer, shaped like graph.linkedTo(). */
export function extLinked(ix, layer, n, budget) {
  const side = layer === "p2" ? ix.p2 : layer === "bp" ? ix.bp : null;
  const nb = side?.adj.get(n.i) || [];
  const all = new Set([n.i]);
  const rest = [];
  for (const x of nb) {
    if (all.has(x.id)) continue;
    all.add(x.id);
    const m = ix.d.byId.get(x.id);
    if (m) rest.push(m);
  }
  return { anc: [], rest: rest.slice(0, budget === Infinity ? rest.length : budget), all };
}

export const sideFor = (ix, layer) => (layer === "p2" ? ix?.p2 : layer === "bp" ? ix?.bp : null);
