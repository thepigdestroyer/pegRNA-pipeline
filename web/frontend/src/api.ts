import type {
  DesignConfig,
  DesignResult,
  Health,
  Mode,
  MotifGroups,
  ProgressEvent,
} from "./types";

/**
 * Where the pipeline API lives.
 *
 * Empty (the default) means same-origin -- the app served by the FastAPI box
 * itself. Set VITE_API_BASE at build time to point a separately hosted frontend
 * (e.g. on Vercel) at the machine running the models. It has to be a direct
 * browser->backend call: a real design runs for minutes, well past the
 * execution cap on a serverless proxy.
 */
const API_BASE = (import.meta.env.VITE_API_BASE ?? "").replace(/\/$/, "");

/** Resolve an API path against the configured backend. */
export function apiURL(path: string): string {
  return `${API_BASE}${path}`;
}

const TOKEN_KEY = "pegrna.accessToken";

/** The access token for `api` mode, kept in localStorage between visits. */
export function getToken(): string {
  try {
    return localStorage.getItem(TOKEN_KEY) ?? "";
  } catch {
    return "";
  }
}

export function setToken(value: string): void {
  try {
    if (value) localStorage.setItem(TOKEN_KEY, value);
    else localStorage.removeItem(TOKEN_KEY);
  } catch {
    /* private browsing can block storage; the token still works for this run */
  }
}

/** Human-readable text for the guard responses web/app.py can return. */
export function statusMessage(code: number): string {
  switch (code) {
    case 429:
      return "The server is already running its maximum number of designs. Try again in a moment.";
    case 401:
      return "This deployment needs an access token, and the one supplied was missing or wrong.";
    case 413:
      return "That configuration is too large for the server to accept.";
    case 400:
      return "The server rejected these inputs. Check the target window and site S.";
    default:
      return `The request failed (HTTP ${code}).`;
  }
}

async function getJSON<T>(url: string): Promise<T> {
  const r = await fetch(apiURL(url));
  if (!r.ok) throw new Error(statusMessage(r.status));
  return (await r.json()) as T;
}

export const getHealth = () => getJSON<Health>("/api/health");
export const getModes = () => getJSON<{ modes: Mode[]; default: Mode }>("/api/modes");
export const getExample = () => getJSON<DesignConfig>("/api/example");
export const getRealExample = () => getJSON<DesignConfig>("/api/example/real");
export const getMotifs = () => getJSON<{ groups: MotifGroups }>("/api/motifs");

export function runFileURL(runId: string, name: string): string {
  return apiURL(`/api/runs/${encodeURIComponent(runId)}/${encodeURIComponent(name)}`);
}

/**
 * Stream a design, calling `onProgress` for each gate as it is evaluated.
 *
 * The endpoint returns NDJSON: zero or more progress objects, then a single
 * `{type:"result", ...}` line. A non-OK status carries a FastAPI `detail`.
 */
export async function streamDesign(
  config: DesignConfig,
  mode: Mode,
  onProgress: (ev: ProgressEvent) => void,
  signal?: AbortSignal,
): Promise<DesignResult> {
  const r = await fetch(apiURL("/api/design/stream"), {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      ...(getToken() ? { "X-Access-Token": getToken() } : {}),
    },
    body: JSON.stringify({ config, mode }),
    signal,
  });

  if (!r.ok) {
    let detail = "";
    try {
      detail = ((await r.json()) as { detail?: string }).detail ?? "";
    } catch {
      try {
        detail = await r.text();
      } catch {
        /* no body to read */
      }
    }
    throw new Error(`${statusMessage(r.status)}${detail ? ` ${String(detail).slice(0, 300)}` : ""}`);
  }
  if (!r.body) throw new Error("The server returned no response stream.");

  const reader = r.body.getReader();
  const decoder = new TextDecoder();
  let buffer = "";
  let final: DesignResult | null = null;

  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    buffer += decoder.decode(value, { stream: true });
    let nl: number;
    while ((nl = buffer.indexOf("\n")) >= 0) {
      const line = buffer.slice(0, nl).trim();
      buffer = buffer.slice(nl + 1);
      if (!line) continue;
      let obj: unknown;
      try {
        obj = JSON.parse(line);
      } catch {
        continue; // a partial or non-JSON line is not fatal to the stream
      }
      const kind = (obj as { type?: string }).type;
      if (kind === "result") final = obj as DesignResult;
      // Keepalives exist only to stop an idle proxy dropping a long design.
      else if (kind !== "heartbeat") onProgress(obj as ProgressEvent);
    }
  }

  if (!final) throw new Error("The design stream ended without returning a result.");
  return final;
}
