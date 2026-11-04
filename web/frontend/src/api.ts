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

// Nothing asks for an access token any more. Earlier versions kept one in localStorage; clear
// any copy a returning visitor's browser still holds rather than leave a secret lying around.
try {
  localStorage.removeItem("pegrna.accessToken");
} catch {
  /* storage can be blocked (private browsing); there is nothing to clear then */
}

/** Human-readable text for the guard responses web/app.py can return. */
export function statusMessage(code: number): string {
  switch (code) {
    case 429:
      return "The server is already running its maximum number of designs. Try again in a moment.";
    case 413:
      return "That configuration is too large for the server to accept.";
    case 400:
      return "The server rejected these inputs. Check the target window and site S.";
    default:
      return `The request failed (HTTP ${code}).`;
  }
}

/**
 * The design backend is a lab machine behind a tunnel, so "unreachable" is a
 * normal state, not an exceptional one. `fetch` reports it as the famously
 * unhelpful "Failed to fetch", so translate it into something a user can act on.
 */
export const BACKEND_UNREACHABLE =
  "Cannot reach the design server. It may be offline, or its address may have changed.";

function isNetworkError(err: unknown): boolean {
  return err instanceof TypeError;
}

async function getJSON<T>(url: string): Promise<T> {
  let r: Response;
  try {
    r = await fetch(apiURL(url));
  } catch (err) {
    if (isNetworkError(err)) throw new Error(BACKEND_UNREACHABLE);
    throw err;
  }
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
  let r: Response;
  try {
    r = await fetch(apiURL("/api/design/stream"), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ config, mode }),
      signal,
    });
  } catch (err) {
    if (isNetworkError(err)) throw new Error(BACKEND_UNREACHABLE);
    throw err;
  }

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
