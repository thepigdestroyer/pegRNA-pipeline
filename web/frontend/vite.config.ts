import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";

// Build straight into web/static/, which FastAPI already serves (app.py mounts
// /static and reads static/index.html for "/"). Keeping that contract means the
// deploy box needs no Node: the built assets are committed.
export default defineConfig({
  plugins: [react(), tailwindcss()],
  // FastAPI mounts the app under /static, so that is the default. A standalone
  // static host (Vercel) serves the same files from the root instead, and sets
  // VITE_BASE=/ — with the wrong base the bundle 404s and nothing renders.
  base: process.env.VITE_BASE ?? "/static/",
  build: {
    outDir: "../static",
    emptyOutDir: false, // index.html is replaced; nothing else lives there
    rollupOptions: {
      output: {
        entryFileNames: "assets/[name]-[hash].js",
        chunkFileNames: "assets/[name]-[hash].js",
        assetFileNames: "assets/[name]-[hash][extname]",
      },
    },
  },
  server: {
    // `npm run dev` proxies the API to a locally running uvicorn.
    proxy: {
      "/api": "http://127.0.0.1:8092",
    },
  },
});
