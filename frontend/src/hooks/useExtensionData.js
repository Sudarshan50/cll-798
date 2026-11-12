import { useEffect, useState } from "react";

const FILES = { phase2: "/phase2.json", bipartite: "/bipartite.json" };

// Survives remounts and mode switches: each file is fetched at most once per page.
const cache = {};

async function load(url, signal) {
  const res = await fetch(url, { signal });
  if (!res.ok) return { status: "missing", data: null };
  // the dev server answers a missing public file with index.html and a 200
  if (!(res.headers.get("content-type") || "").includes("json"))
    return { status: "missing", data: null };
  const data = await res.json();
  if (!data || typeof data !== "object") return { status: "missing", data: null };
  return { status: "ready", data };
}

/**
 * Fetch phase2.json and bipartite.json the first time `enabled` is true, then keep them.
 * Each file reports "idle" | "loading" | "ready" | "missing" | "error" on its own, so one
 * absent analysis never hides the other.
 */
export function useExtensionData(enabled) {
  const [state, setState] = useState(() => ({
    phase2: cache.phase2 || { status: "idle", data: null },
    bipartite: cache.bipartite || { status: "idle", data: null },
  }));

  useEffect(() => {
    if (!enabled) return undefined;
    const ac = new AbortController();
    for (const [key, url] of Object.entries(FILES)) {
      if (cache[key]) continue;
      setState((s) => ({ ...s, [key]: { status: "loading", data: null } }));
      load(url, ac.signal).then((r) => {
        cache[key] = r;
        setState((s) => ({ ...s, [key]: r }));
      }).catch((e) => {
        if (e.name === "AbortError") return;      // retried on the next enable
        console.error(`Extension data ${url} failed:`, e);
        const r = { status: "error", data: null };
        cache[key] = r;
        setState((s) => ({ ...s, [key]: r }));
      });
    }
    return () => ac.abort();
  }, [enabled]);

  return state;
}
