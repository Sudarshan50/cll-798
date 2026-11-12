export const fmt = (v) => (v ?? 0).toLocaleString("en-US");

const finite = (v) => typeof v === "number" && Number.isFinite(v);

/** A measured quantity, or an em dash when it was not computed. Integers keep `fmt`. */
export const num = (v, digits = 3) => {
  if (!finite(v)) return "—";
  if (Number.isInteger(v) && Math.abs(v) >= 1) return fmt(v);
  const a = Math.abs(v);
  if (a >= 1000) return fmt(Math.round(v));
  if (a !== 0 && a < 1e-3) return v.toExponential(1);
  return Number(v.toPrecision(digits)).toString();
};

/** Anything from a results file as display text; objects must never reach React as children. */
export const text = (v) =>
  v == null ? "—" : typeof v === "object" ? JSON.stringify(v) : String(v);

export const pval = (p) => (!finite(p) ? "—" : p < 1e-3 ? "< 0.001" : p.toFixed(3));

/** 11646523 -> "11.6M", for axis ticks that must stay narrow. */
export const compact = (v) => {
  if (!finite(v)) return "—";
  const a = Math.abs(v);
  if (a >= 1e9) return `${+(v / 1e9).toFixed(1)}B`;
  if (a >= 1e6) return `${+(v / 1e6).toFixed(1)}M`;
  if (a >= 1e4) return `${+(v / 1e3).toFixed(0)}k`;
  if (a >= 1e3) return `${+(v / 1e3).toFixed(1)}k`;
  return num(v, 2);
};
