import { useEffect, useState } from "react";
import { PointCloud } from "../gl.js";
import { buildIndexes, buildPalette } from "../core/graph.js";
import { sizeOverlay } from "../render/scene.js";

/**
 * Fetch the layout, frame it, then stream the compound cloud into the GPU.
 * Fills `dataRef` / `cloudRef` and drives the camera's initial framing.
 */
export function useAtlasSource({ enabled, wrapRef, glRef, ovRef, camera, dataRef, cloudRef, schedule }) {
  const [status, setStatus] = useState("loading");
  const [loaded, setLoaded] = useState(0);
  const [meta, setMeta] = useState(null);
  const [err, setErr] = useState(null);

  useEffect(() => {
    // Nothing is fetched on an unsupported device. The gate renders instead of the
    // canvas, so without this the 292 MB stream would still start behind the notice.
    if (!enabled) return undefined;
    let dead = false; const ac = new AbortController();
    (async () => {
      try {
        const d = await (await fetch("/layout.json", { signal: ac.signal })).json();
        if (dead) return;
        const col = buildPalette(d.nodes);
        d.nodes.forEach((n) => { n.col = col(n); });
        buildIndexes(d);
        dataRef.current = d; setMeta(d.meta);

        const dpr = Math.min(2, window.devicePixelRatio || 1);
        const w = wrapRef.current.clientWidth, h = wrapRef.current.clientHeight;
        const ext = Math.max(...d.nodes.map((n) => Math.max(Math.abs(n.px), Math.abs(n.py))));
        camera.set({ pan: [w / 2, h / 2], zoom: (Math.min(w, h) * 0.46) / ext, dpr, w, h });

        cloudRef.current = new PointCloud(glRef.current);
        cloudRef.current.resize(w, h, dpr);
        sizeOverlay(ovRef.current, w, h, dpr);
        setStatus("ready"); schedule();

        if (!cloudRef.current.allocate(d.meta.points))
          throw new Error(`GPU refused a ${(d.meta.points * 4 / 1e6) | 0} MB buffer`);
        setStatus("streaming compounds");
        await cloudRef.current.stream(`/points.bin?v=${d.meta.version ?? "0"}`, d.meta.points, (n) => {
          if (!dead) setLoaded(n);
        }, ac.signal);
        if (dead) return;
        setStatus("ready"); schedule();
      } catch (e) {
        if (!dead && e.name !== "AbortError") {
          // the banner stays generic; the detail belongs in the console, not on screen
          console.error("Atlas failed to load:", e);
          setErr(e.message || String(e));
        }
      }
    })();
    return () => { dead = true; ac.abort(); };
  }, [enabled, wrapRef, glRef, ovRef, camera, dataRef, cloudRef, schedule]);

  return { status, loaded, meta, err };
}
