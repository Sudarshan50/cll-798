import { useCallback, useEffect, useRef } from "react";

/**
 * One rAF-coalesced draw. `busy` is true while the user is interacting, so layers
 * can take a cheap pass. `schedule` is stable, so effects may depend on it.
 */
export function useRenderLoop(render) {
  const raf = useRef(0), idle = useRef(0), busy = useRef(false);
  const latest = useRef(render);
  latest.current = render;

  const schedule = useCallback(() => {
    cancelAnimationFrame(raf.current);
    raf.current = requestAnimationFrame(() => latest.current());
  }, []);

  const interact = useCallback(() => {
    busy.current = true; schedule();
    clearTimeout(idle.current);
    idle.current = setTimeout(() => { busy.current = false; schedule(); }, 130);
  }, [schedule]);

  useEffect(() => () => { cancelAnimationFrame(raf.current); clearTimeout(idle.current); }, []);
  return { schedule, interact, busy };
}
