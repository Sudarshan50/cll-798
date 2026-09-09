import React, { useCallback, useEffect, useRef, useState } from "react";
import { Camera, clampZoom, zoomAbout } from "./core/camera.js";
import { labelBudget } from "./core/labels.js";
import { boundsOf, linkedTo } from "./core/graph.js";
import { renderScene, sizeOverlay } from "./render/scene.js";
import { POINT_ZOOM, renderPoints } from "./render/points.js";
import { pickNode } from "./core/hitTest.js";
import { useRenderLoop } from "./hooks/useRenderLoop.js";
import { useAtlasSource } from "./hooks/useAtlasSource.js";
import { Hud } from "./components/Hud.jsx";
import { InfoPanel } from "./components/InfoPanel.jsx";
import { ErrorBanner, Hint, LabelToggle, ModeToggle, SelectionBar } from "./components/Chrome.jsx";

export default function App() {
  const glRef = useRef(null), ovRef = useRef(null), wrapRef = useRef(null);
  const dataRef = useRef(null), cloudRef = useRef(null), hoverRef = useRef(null);
  const cameraRef = useRef(null);
  if (!cameraRef.current) cameraRef.current = new Camera();
  const camera = cameraRef.current;
  const press = useRef(null), drag = useRef(null), anim = useRef(0);

  const [mode, setMode] = useState("network");     // "network" | "compounds"
  const [info, setInfo] = useState(null);   // right-click -> pinned info panel
  const [onScreen, setOnScreen] = useState(0);
  const [labels, setLabels] = useState(true);
  // camera is a ref, so zoom changes never re-render. Mirror it into state
  // (from inside the draw, where the live value is known) or the readout freezes and
  // the zoom-scaled label cap looks like it is doing nothing.
  const [ui, setUi] = useState({ zoom: 0, named: 0, budget: 0, noRoom: 0, off: 0 });
  const [sel, setSel] = useState(null);      // the left-clicked class

  const { schedule, interact, busy } = useRenderLoop(() => {
    const d = dataRef.current, cloud = cloudRef.current;
    if (cloud && d) {
      let drawn = 0;
      if (mode === "compounds") drawn = renderPoints(cloud, d, camera.state);
      else cloud.clear();
      setOnScreen((prev) => (Math.abs(prev - drawn) / Math.max(drawn, 1) < 0.01 ? prev : drawn));
    }
    if (!d || !ovRef.current) return;
    const view = camera.state;
    const hv = hoverRef.current;                   // the left-click selection
    const lod = busy.current;                      // cheap pass while interacting
    const budget = labelBudget(view.zoom);
    const link = hv ? linkedTo(d, hv, budget) : null;
    const { named, noRoom, offScreen } = renderScene(ovRef.current, {
      d, view, hv, link, lit: link ? link.all : null, lod, showLabels: labels,
    });
    if (!lod) {
      const nz = Math.round(view.zoom), nb = budget === Infinity ? -1 : budget;
      setUi((u) => (u.zoom === nz && u.named === named && u.budget === nb
        && u.noRoom === noRoom && u.off === offScreen
        ? u : { zoom: nz, named, budget: nb, noRoom, off: offScreen }));
    }
  });

  const { status, loaded, meta, err } =
    useAtlasSource({ wrapRef, glRef, ovRef, camera, dataRef, cloudRef, schedule });

  // Switching to compound mode below the reveal zoom used to look like nothing
  // happened. Fly to the threshold so the difference is visible at once.
  useEffect(() => {
    if (mode !== "compounds" || camera.zoom >= POINT_ZOOM * 1.05) { schedule(); return; }
    const base = { ...camera.state, pan: [...camera.pan] };
    const z1 = POINT_ZOOM * 1.5, cx = camera.w / 2, cy = camera.h / 2;
    const t0 = performance.now(), dur = 900;
    const step = () => {
      const u = Math.min(1, (performance.now() - t0) / dur);
      const e = 1 - Math.pow(1 - u, 3);
      camera.set(zoomAbout(base, cx, cy, base.zoom + (z1 - base.zoom) * e));
      schedule();
      if (u < 1 && anim.current) anim.current = requestAnimationFrame(step);
    };
    anim.current = requestAnimationFrame(step);
    return () => { cancelAnimationFrame(anim.current); anim.current = 0; };
  }, [mode, schedule, camera]);

  useEffect(() => { schedule(); }, [labels, schedule]);

  useEffect(() => {
    const onR = () => {
      const wrap = wrapRef.current; if (!wrap || !cloudRef.current) return;
      const dpr = Math.min(2, window.devicePixelRatio || 1);
      const w = wrap.clientWidth, h = wrap.clientHeight;
      camera.resize(w, h, dpr);
      cloudRef.current.resize(w, h, dpr);
      sizeOverlay(ovRef.current, w, h, dpr);
      schedule();
    };
    addEventListener("resize", onR); return () => removeEventListener("resize", onR);
  }, [schedule, camera]);

  const pick = useCallback((e) => {
    const d = dataRef.current;
    if (!d) return null;
    const r = e.currentTarget.getBoundingClientRect();
    return pickNode(d.nodes, camera.state, e.clientX - r.left, e.clientY - r.top);
  }, [camera]);

  const onMove = (e) => {
    if (drag.current) {
      camera.panBy(e.clientX - drag.current.x, e.clientY - drag.current.y);
      drag.current = { x: e.clientX, y: e.clientY };
      interact(); return;
    }
    // hover no longer selects or opens anything; it only signals that a node is
    // clickable. Use a class, not an inline style, or it overrides :active grabbing.
    e.currentTarget.classList.toggle("over", !!pick(e));
  };

  const onDown = (e) => {
    cancelAnimationFrame(anim.current); anim.current = 0;
    press.current = { x: e.clientX, y: e.clientY, b: e.button };
    if (e.button === 0) {
      drag.current = { x: e.clientX, y: e.clientY };
      e.currentTarget.setPointerCapture(e.pointerId);
    }
  };

  const onUp = (e) => {
    const p = press.current; press.current = null;
    if (e.button === 0) {
      drag.current = null;
      try { e.currentTarget.releasePointerCapture(e.pointerId); } catch {}
    }
    if (!p || p.b !== 0) return;
    const moved = Math.hypot(e.clientX - p.x, e.clientY - p.y);
    if (moved > 4) return;                       // that was a pan, not a click
    const n = pick(e);
    hoverRef.current = n;                        // left click toggles the highlight
    setSel(n);
    if (!n) setInfo(null);
    schedule();
  };

  const onContext = (e) => {
    e.preventDefault();
    const n = pick(e);
    setInfo(n ? { n, x: e.clientX, y: e.clientY } : null);
  };

  const onWheel = (e) => {
    cancelAnimationFrame(anim.current); anim.current = 0;   // user overrides the fly-in
    const r = e.currentTarget.getBoundingClientRect();
    camera.zoomTo(e.clientX - r.left, e.clientY - r.top,
                  clampZoom(camera.zoom * Math.exp(-e.deltaY * 0.0016)));
    interact();
  };

  /** Zoom and pan so the selection and all of its links fit on screen. */
  const fitSelection = () => {
    const d = dataRef.current, n = hoverRef.current;
    if (!d || !n) return;
    const b = boundsOf(d, linkedTo(d, n, Infinity).all);
    if (!b) return;
    cancelAnimationFrame(anim.current); anim.current = 0;
    camera.fitBounds(b.x0, b.y0, b.x1, b.y1);
    schedule();
  };

  return (
    <div className="app">
      <div ref={wrapRef} className="stage"
        onPointerDown={onDown}
        onPointerUp={onUp}
        onPointerMove={onMove}
        onContextMenu={onContext}
        onWheel={onWheel}>
        <canvas ref={glRef} className="layer" />
        <canvas ref={ovRef} className="layer overlay" />
      </div>

      {err && <ErrorBanner message={err} />}

      <ModeToggle mode={mode} onChange={setMode} />

      <Hud mode={mode} onScreen={onScreen} loaded={loaded} meta={meta}
           bundles={dataRef.current?.bundles?.length} status={status} ui={ui} selected={sel} />

      <LabelToggle labels={labels} onChange={setLabels} />

      {info && <InfoPanel n={info.n} x={info.x} y={info.y} d={dataRef.current}
                          onClose={() => setInfo(null)} />}

      {sel && <SelectionBar sel={sel} onFit={fitSelection} onClear={() => {
        hoverRef.current = null; setSel(null); setInfo(null); schedule();
      }} />}

      <Hint />
    </div>
  );
}
