/**
 * Reading a design's progress stream, and giving up on one that has gone silent.
 *
 * The server answers a design with NDJSON: progress lines, a keepalive line every 15 s for as long as the pipeline is
 * quiet (PEGRNA_HEARTBEAT_SECONDS), then one `{"type":"result"}` line. A stream that is alive is therefore never silent
 * for long, and a connection that silently stops delivering (it happened: the server finished a design and wrote its
 * result while the page, which had stopped receiving data a minute in, went on waiting) used to look like a frozen
 * page, for ever. Six missed keepalives in a row is a connection that has gone, not one that is slow.
 *
 * This file imports only types, so a plain Node script can run it (scripts/qa-frontend-logic.mjs).
 */
import type { DesignResult, ProgressEvent } from "./types";

/**
 * How long a stream may stay completely silent before the page stops waiting. It has to stay several times the
 * server's keepalive interval (15 s), or a healthy design would be given up on.
 */
export const STREAM_STALL_MS = 90_000;

/** What the person is told. The design may well have finished on the server, so it says so and what to do. */
export function stalledMessage(ms: number): string {
  return (
    `Lost contact with the design server: nothing arrived for ${Math.round(ms / 1000)} seconds, although it sends ` +
    "a signal every 15. The design may have finished there. Run it again: the same inputs give the same result."
  );
}

/** The part of a stream reader this needs, so that a test can stand in for the network. */
export interface ChunkReader {
  read(): Promise<{ done: boolean; value?: Uint8Array }>;
  cancel(): Promise<void>;
}

/**
 * Read the stream to its end and return the final result, calling `onProgress` for each progress line. Rejects with
 * `stalledMessage` if no data at all (a keepalive counts) arrives for `stallMs`, and with whatever the reader rejects
 * with (an aborted fetch, a dropped connection) unchanged.
 */
export async function readDesignStream(
  reader: ChunkReader,
  onProgress: (ev: ProgressEvent) => void,
  stallMs: number = STREAM_STALL_MS,
): Promise<DesignResult> {
  const decoder = new TextDecoder();
  let buffer = "";
  let final: DesignResult | null = null;

  for (;;) {
    // A fresh deadline for every wait: it is the silence that is limited, not the length of the design.
    let timer: ReturnType<typeof setTimeout> | undefined;
    const silence = new Promise<"silent">((resolve) => {
      timer = setTimeout(() => resolve("silent"), stallMs);
    });
    let chunk: { done: boolean; value?: Uint8Array } | "silent";
    try {
      chunk = await Promise.race([reader.read(), silence]);
    } finally {
      clearTimeout(timer);
    }
    if (chunk === "silent") {
      reader.cancel().catch(() => undefined);
      throw new Error(stalledMessage(stallMs));
    }
    if (chunk.done) break;

    buffer += decoder.decode(chunk.value, { stream: true });
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
