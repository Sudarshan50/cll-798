import { useEffect, useState } from "react";

export const MIN_W = 1024;
export const MIN_H = 600;

const ok = () =>
  window.innerWidth >= MIN_W &&
  window.innerHeight >= MIN_H &&
  // a phone or tablet has no way to right-click the info panel, and cannot
  // allocate the 292 MB vertex buffer the compound mode needs
  !(window.matchMedia("(pointer: coarse)").matches &&
    !window.matchMedia("(any-pointer: fine)").matches);

/** Whether this device can run the atlas at all. */
export default function useViewportGate() {
  const [supported, setSupported] = useState(ok);
  useEffect(() => {
    const on = () => setSupported(ok());
    addEventListener("resize", on);
    addEventListener("orientationchange", on);
    return () => {
      removeEventListener("resize", on);
      removeEventListener("orientationchange", on);
    };
  }, []);
  return supported;
}
