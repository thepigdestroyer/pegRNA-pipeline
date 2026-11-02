import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";

// Build straight into web/static/, which FastAPI already serves (app.py mounts
// /static and reads static/index.html for "/"). Keeping that contract means the
// deploy box needs no Node: the built assets are committed.
export default defineConfig({
  plugins: [react(), tailwindcss()],
  base: "/static/",
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
