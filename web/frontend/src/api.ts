import type {
  DesignConfig,
  DesignResult,
  Health,
  Mode,
  MotifGroups,
  ProgressEvent,
} from "./types";

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
  const r = await fetch(url);
  if (!r.ok) throw new Error(statusMessage(r.status));
  return (await r.json()) as T;
}

export const getHealth = () => getJSON<Health>("/api/health");
export const getModes = () => getJSON<{ modes: Mode[]; default: Mode }>("/api/modes");
export const getExample = () => getJSON<DesignConfig>("/api/example");
export const getRealExample = () => getJSON<DesignConfig>("/api/example/real");
export const getMotifs = () => getJSON<{ groups: MotifGroups }>("/api/motifs");

export function runFileURL(runId: string, name: string): string {
  return `/api/runs/${encodeURIComponent(runId)}/${encodeURIComponent(name)}`;
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
  const r = await fetch("/api/design/stream", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
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
      if ((obj as DesignResult).type === "result") final = obj as DesignResult;
      else onProgress(obj as ProgressEvent);
    }
  }

  if (!final) throw new Error("The design stream ended without returning a result.");
  return final;
}
