import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  server: { port: 5173, open: true },
  // points.bin is ~292 MB: never inline it, and don't let esbuild touch it
  build: { assetsInlineLimit: 0, chunkSizeWarningLimit: 4000 },
  assetsInclude: ["**/*.bin"],
});
