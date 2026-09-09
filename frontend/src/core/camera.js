export const ZOOM_MIN = 80;
export const ZOOM_MAX = 6e6;
export const clampZoom = (z) => Math.max(ZOOM_MIN, Math.min(ZOOM_MAX, z));

export const toScreenX = (v, x) => x * v.zoom + v.pan[0];
export const toScreenY = (v, y) => y * v.zoom + v.pan[1];

/** { pan, zoom } that lands on `z` while keeping the screen point (sx, sy) over its world point. */
export function zoomAbout(v, sx, sy, z) {
  const s = z / v.zoom;
  return { zoom: z, pan: [sx - (sx - v.pan[0]) * s, sy - (sy - v.pan[1]) * s] };
}

export function boundsToView(v, x0, y0, x1, y1) {
  const zoom = clampZoom(0.82 * Math.min(v.w / Math.max(x1 - x0, 1e-6),
                                         v.h / Math.max(y1 - y0, 1e-6)));
  return { zoom, pan: [v.w / 2 - ((x0 + x1) / 2) * zoom, v.h / 2 - ((y0 + y1) / 2) * zoom] };
}

export class Camera {
  constructor() { this.pan = [0, 0]; this.zoom = 1; this.dpr = 1; this.w = 0; this.h = 0; }

  get state() {
    return { pan: this.pan, zoom: this.zoom, dpr: this.dpr, w: this.w, h: this.h };
  }

  set(v) { Object.assign(this, v); }

  panBy(dx, dy) { this.pan = [this.pan[0] + dx, this.pan[1] + dy]; }

  zoomTo(sx, sy, z) { this.set(zoomAbout(this.state, sx, sy, z)); }

  fitBounds(x0, y0, x1, y1) { this.set(boundsToView(this.state, x0, y0, x1, y1)); }

  resize(w, h, dpr) {
    this.pan = [this.pan[0] + (w - this.w) / 2, this.pan[1] + (h - this.h) / 2];
    this.dpr = dpr; this.w = w; this.h = h;
  }
}
