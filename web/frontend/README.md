# pegRNA designer — frontend

React + TypeScript + Tailwind, built with Vite. The build output goes to
`../static/`, which `web/app.py` already serves (`/` reads `static/index.html`,
`/static` is mounted). **The built assets are committed**, so a deploy box needs
no Node — only this directory does.

## Working on it

```bash
npm install

# Terminal 1 — the API (from the repo root)
PEGRNA_REPO="$PWD" uvicorn app:app --app-dir web --port 8092

# Terminal 2 — the UI, with /api proxied to that server
npm run dev
```

## Shipping a change

```bash
npm run build      # cleans ../static/assets, type-checks, then builds
```

Commit the regenerated `web/static/index.html` and `web/static/assets/*`
alongside your source change, or the deployment will serve the old bundle.

## Layout

| Path | What |
|---|---|
| `src/App.tsx` | Page layout, form state, validation, submit |
| `src/api.ts` | Typed client for `/api/*`, including the NDJSON progress stream |
| `src/types.ts` | Shapes returned by `web/app.py` |
| `src/pipeline.ts` | The gate cascade, described in plain language |
| `src/components/SequenceMap.tsx` | The pegRNA drawn to scale, by part |
| `src/components/PipelinePanel.tsx` | Live per-gate progress and measured values |
| `src/components/ResultPanel.tsx` | Result, downloads, stand-in disclosure |
| `src/components/CargoBuilder.tsx` | Motif-library cargo assembly |
| `src/components/Primitives.tsx` | Field, Badge, Callout, Disclosure, InfoHint |

## Things worth knowing

- **Gate ids are a contract.** The `id` values in `pipeline.ts` match what the
  CLI emits on its `@@P` progress lines and writes into `gates.json`. Renaming
  one silently stops that step lighting up.
- **Stand-in gates must stay visible.** `gates.stand_in_gates` lists the checks
  that returned a placeholder instead of running a real model (api mode with no
  NVIDIA key, or no genome FASTA). The pipeline badges them and the result
  panel leads with a callout. Do not bury this: people order DNA off these
  reports.
- **An empty form is incomplete, not invalid.** Errors appear only once a field
  is touched or the design is submitted.
