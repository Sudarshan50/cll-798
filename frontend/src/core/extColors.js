/** Colour schemes shared by the extension map layers and their legends. */

export const NO_DATA = "rgba(107,120,145,.38)";

/** Edge stroke colours as rgb triples, keyed by layer then by source code. */
export const EDGE_RGB = {
  p2: { 1: [116, 152, 224], 2: [178, 140, 232], 3: [242, 112, 96] },
  bp: { 0: [79, 195, 176] },
};
export const EDGE_LABELS = {
  p2: { 1: "S1 · direct parent ↔ alt", 2: "S2 · alt ↔ alt", 3: "S3 · rule violation" },
  bp: { 0: "BiCM-validated co-labelling" },
};
export const HILITE_RGB = [255, 205, 140];
export const FLAG_RING = "#F2A65A";

/** Communities are ranked by size; past this many they share one neutral colour. */
export const COMMUNITY_COLOURS = 16;
export const communityColour = (rank) =>
  rank < COMMUNITY_COLOURS
    ? `hsl(${((rank * 137.508 + 24) % 360).toFixed(0)} 60% 63%)`
    : "hsl(220 12% 46%)";

/** Sequential ramp for ambiguity H: low = deep indigo, high = pale amber. */
export const H_STOPS = [[43, 58, 103], [62, 128, 170], [79, 195, 176], [242, 166, 90], [255, 231, 194]];
export function rampRgb(t) {
  const u = Math.max(0, Math.min(1, t)) * (H_STOPS.length - 1);
  const i = Math.min(H_STOPS.length - 2, Math.floor(u)), f = u - i;
  const a = H_STOPS[i], b = H_STOPS[i + 1];
  return `rgb(${(a[0] + (b[0] - a[0]) * f) | 0},${(a[1] + (b[1] - a[1]) * f) | 0},${(a[2] + (b[2] - a[2]) * f) | 0})`;
}
export const H_GRADIENT = `linear-gradient(90deg,${H_STOPS.map((c) => `rgb(${c})`).join(",")})`;
