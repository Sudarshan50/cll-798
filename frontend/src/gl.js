/**
 * Streaming WebGL2 point-cloud renderer.
 *
 * The whole cloud is 73,105,281 points x 4 bytes = 292 MB. On an 8 GB machine we must
 * never hold that twice, so the file is read as a stream and each chunk is uploaded
 * straight into a pre-allocated GPU buffer and then dropped. Peak CPU memory stays at
 * one chunk; the GPU holds the rest.
 */

const VERT = `#version 300 es
in vec2 a_pos;                 // int16 world coords
uniform vec2  u_res;
uniform vec2  u_pan;
uniform float u_zoom;
uniform float u_scale;         // int16 -> world divisor
uniform float u_size;
void main(){
  vec2 world = a_pos / u_scale;
  vec2 px    = world * u_zoom + u_pan;
  vec2 clip  = (px / u_res) * 2.0 - 1.0;
  gl_Position = vec4(clip.x, -clip.y, 0.0, 1.0);
  gl_PointSize = u_size;
}`;

const FRAG = `#version 300 es
precision mediump float;
uniform vec3  u_col;
uniform float u_alpha;
out vec4 outColor;
void main(){
  vec2 d = gl_PointCoord - 0.5;
  float f = max(0.0, 1.0 - dot(d, d) * 4.0);
  outColor = vec4(u_col * f, u_alpha * f);
}`;

function compile(gl, type, src) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS))
    throw new Error(gl.getShaderInfoLog(s) || "shader compile failed");
  return s;
}

export class PointCloud {
  constructor(canvas) {
    const gl = canvas.getContext("webgl2", {
      alpha: false, antialias: false, depth: false,
      powerPreference: "high-performance", preserveDrawingBuffer: false,
    });
    if (!gl) throw new Error("WebGL2 is not available in this browser");
    this.gl = gl;
    this.canvas = canvas;
    this.count = 0;
    this.capacity = 0;

    const p = gl.createProgram();
    gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, VERT));
    gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, FRAG));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS))
      throw new Error(gl.getProgramInfoLog(p) || "link failed");
    this.prog = p;
    this.u = {
      res: gl.getUniformLocation(p, "u_res"),
      pan: gl.getUniformLocation(p, "u_pan"),
      zoom: gl.getUniformLocation(p, "u_zoom"),
      scale: gl.getUniformLocation(p, "u_scale"),
      size: gl.getUniformLocation(p, "u_size"),
      col: gl.getUniformLocation(p, "u_col"),
      alpha: gl.getUniformLocation(p, "u_alpha"),
    };

    this.vao = gl.createVertexArray();
    this.buf = gl.createBuffer();
    gl.bindVertexArray(this.vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buf);
    const loc = gl.getAttribLocation(p, "a_pos");
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.SHORT, false, 4, 0);
    gl.bindVertexArray(null);

    gl.disable(gl.DEPTH_TEST);
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE);   // additive: overlap reads as density
  }

  /** Pre-allocate the full buffer. Returns false if the driver refuses the size. */
  allocate(points) {
    const gl = this.gl;
    const bytes = points * 4;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buf);
    while (gl.getError() !== gl.NO_ERROR) { /* drain */ }
    gl.bufferData(gl.ARRAY_BUFFER, bytes, gl.STATIC_DRAW);
    if (gl.getError() !== gl.NO_ERROR) return false;
    this.capacity = points;
    return true;
  }

  /**
   * Stream the binary in, uploading each chunk and discarding it.
   * onProgress(pointsUploaded) is called per chunk so the UI can count up live.
   */
  async stream(url, total, onProgress, signal) {
    const gl = this.gl;
    const res = await fetch(url, { signal });
    if (!res.ok) throw new Error(`points.bin: HTTP ${res.status}`);
    const reader = res.body.getReader();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buf);

    let offset = 0;          // byte offset into the GPU buffer
    let carry = null;        // bytes left over from a chunk that split a point
    let sinceYield = 0;
    const YIELD_BYTES = 8 << 20;   // yielding per network chunk cost ~18s of clamped
                                   // setTimeout and thousands of React renders
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      let bytes = value;
      if (carry && carry.length) {
        const merged = new Uint8Array(carry.length + bytes.length);
        merged.set(carry, 0);
        merged.set(bytes, carry.length);
        bytes = merged;
        carry = null;
      }
      const usable = bytes.length - (bytes.length % 4);
      if (usable < bytes.length) carry = bytes.slice(usable);
      if (!usable) continue;
      const room = this.capacity * 4 - offset;
      const n = Math.min(usable, room);
      if (n <= 0) break;
      gl.bufferSubData(gl.ARRAY_BUFFER, offset, bytes.subarray(0, n));
      offset += n;
      this.count = offset / 4;
      sinceYield += n;
      if (sinceYield >= YIELD_BYTES) {
        sinceYield = 0;
        onProgress(this.count, total);
        // one yield per 8 MB keeps the frame loop alive without stalling the transfer
        await new Promise((r) => setTimeout(r, 0));
      }
      if (offset >= this.capacity * 4) {
        try { await reader.cancel(); } catch { /* already closed */ }
        break;
      }
    }
    onProgress(this.count, total);
    return this.count;
  }

  resize(w, h, dpr) {
    this.canvas.width = Math.floor(w * dpr);
    this.canvas.height = Math.floor(h * dpr);
    this.gl.viewport(0, 0, this.canvas.width, this.canvas.height);
  }

  clear() {
    const gl = this.gl;
    gl.clearColor(0.031, 0.043, 0.071, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
  }

  /**
   * Draw only the classes whose ranges are given, as [offset, count] pairs.
   * One call per visible class keeps on-screen density at a few points per pixel
   * instead of drawing all 73.1M into whatever happens to be in view.
   */
  drawRanges(ranges, opts) {
    const gl = this.gl;
    if (!this.count || !ranges.length) return 0;
    const { pan, zoom, scale, dpr, size, alpha, color } = opts;
    gl.useProgram(this.prog);
    gl.bindVertexArray(this.vao);
    gl.uniform2f(this.u.res, this.canvas.width, this.canvas.height);
    gl.uniform2f(this.u.pan, pan[0] * dpr, pan[1] * dpr);
    gl.uniform1f(this.u.zoom, zoom * dpr);
    gl.uniform1f(this.u.scale, scale);
    gl.uniform1f(this.u.size, size * dpr);
    gl.uniform1f(this.u.alpha, alpha);
    gl.uniform3f(this.u.col, color[0], color[1], color[2]);
    let drawn = 0;
    for (const [off, n] of ranges) {
      if (off >= this.count) continue;
      const c = Math.min(n, this.count - off);
      if (c <= 0) continue;
      gl.drawArrays(gl.POINTS, off, c);
      drawn += c;
    }
    gl.bindVertexArray(null);
    return drawn;
  }

  draw({ pan, zoom, scale, dpr, size, alpha, color }) {
    const gl = this.gl;
    gl.clearColor(0.031, 0.043, 0.071, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
    if (!this.count) return;
    gl.useProgram(this.prog);
    gl.bindVertexArray(this.vao);
    gl.uniform2f(this.u.res, this.canvas.width, this.canvas.height);
    gl.uniform2f(this.u.pan, pan[0] * dpr, pan[1] * dpr);
    gl.uniform1f(this.u.zoom, zoom * dpr);
    gl.uniform1f(this.u.scale, scale);
    gl.uniform1f(this.u.size, size * dpr);
    gl.uniform1f(this.u.alpha, alpha);
    gl.uniform3f(this.u.col, color[0], color[1], color[2]);
    gl.drawArrays(gl.POINTS, 0, this.count);
    gl.bindVertexArray(null);
  }
}
