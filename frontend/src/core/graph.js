export function buildPalette(nodes) {
  const byId = new Map(nodes.map((n) => [n.i, n]));
  const scs = [...new Set(nodes.map((n) => n.s).filter((v) => v != null))]
    .sort((a, b) => (byId.get(a).k - byId.get(b).k) || a - b);
  const org = scs.filter((s) => byId.get(s).k === 0);
  const ino = scs.filter((s) => byId.get(s).k === 1);
  const hue = new Map();
  org.forEach((s, i) => hue.set(s, (348 + (i / Math.max(1, org.length - 1)) * 82) % 360));
  ino.forEach((s, i) => hue.set(s, 166 + (i / Math.max(1, ino.length - 1)) * 88));
  return (n) => (n.s == null ? (n.k === 0 ? "#F2A65A" : "#4FC3B0")
                             : `hsl(${hue.get(n.s).toFixed(0)} 62% 64%)`);
}

export function buildIndexes(d) {
  d.byId = new Map(d.nodes.map((n) => [n.i, n]));
  d.byCount = [...d.nodes].sort((a, b) => b.dp - a.dp);
  d.cnb = new Map(d.nodes.map((n) => [n.i, []]));
  for (const l of d.links) {
    d.cnb.get(l.s)?.push([l.t, l.w]); d.cnb.get(l.t)?.push([l.s, l.w]);
  }
  d.cnb.forEach((v) => v.sort((a, b) => b[1] - a[1]));
  d.kids = new Map(d.nodes.map((n) => [n.i, []]));
  for (const n of d.nodes) if (n.p != null && n.p >= 0) d.kids.get(n.p)?.push(n);
  return d;
}

export function lineage(d, n) {
  const out = [];
  let c = n, g = 0;
  while (c && g++ < 12) { out.push(c); const p = c.p; if (p == null || p < 0) break; c = d.byId.get(p); }
  return out;
}

/** Ancestors (always named), plus children and confusability partners ranked by size. */
export function linkedTo(d, n, budget) {
  const anc = [];
  let c = d.byId.get(n.p), g = 0;
  while (c && g++ < 12) { anc.push(c); const p = c.p; if (p == null || p < 0) break; c = d.byId.get(p); }
  const kids = d.kids.get(n.i) || [];
  const conf = (d.cnb.get(n.i) || []).map(([id]) => d.byId.get(id)).filter(Boolean);
  const seen = new Set([n.i, ...anc.map((x) => x.i)]);
  const rest = [];
  for (const x of [...kids, ...conf]) {
    if (seen.has(x.i)) continue;
    seen.add(x.i); rest.push(x);
  }
  rest.sort((a, b) => b.dp - a.dp);
  return { anc, rest: rest.slice(0, budget === Infinity ? rest.length : budget), all: seen };
}

export function boundsOf(d, ids) {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const id of ids) {
    const m = d.byId.get(id); if (!m) continue;
    x0 = Math.min(x0, m.px - m.nr); x1 = Math.max(x1, m.px + m.nr);
    y0 = Math.min(y0, m.py - m.nr); y1 = Math.max(y1, m.py + m.nr);
  }
  return isFinite(x0) ? { x0, y0, x1, y1 } : null;
}
