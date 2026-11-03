/**
 * Stable front door for the pegRNA backend.
 *
 * The backend sits behind a free tunnel whose hostname is random and rotates
 * (localhost.run re-rolls it roughly every 20 minutes; a cloudflared quick
 * tunnel re-rolls on every restart). The frontend bakes its backend URL in at
 * build time, so every rotation silently broke the public site: it loaded fine
 * and could not do anything.
 *
 * This Worker gives one permanent hostname. Only the Worker's ORIGIN constant
 * changes when the tunnel rotates, and the backend host updates that itself via
 * the Cloudflare API (web/deploy/update-worker-origin.sh), so nothing needs a
 * redeploy and the site never goes dark.
 *
 * A plain DNS CNAME cannot do this: these tunnel providers route by Host
 * header, so they must receive their OWN hostname, not ours. The Worker
 * rewrites it.
 *
 * The response body is passed straight through, so the NDJSON design stream
 * (and its keepalives) still arrive incrementally.
 */

const ORIGIN = "__ORIGIN__"; // rewritten on deploy

export default {
  async fetch(request) {
    if (!ORIGIN || ORIGIN.includes("__ORIGIN")) {
      return json({ error: "Backend origin is not configured on the proxy." }, 503);
    }

    const incoming = new URL(request.url);
    const origin = new URL(ORIGIN);
    const target = new URL(incoming.pathname + incoming.search, origin);

    // The tunnel providers route on Host, so send theirs, not ours.
    const headers = new Headers(request.headers);
    headers.set("Host", origin.host);
    // Preserve the real client's origin for the backend's CORS check.
    headers.delete("cf-connecting-ip");

    let upstream;
    try {
      upstream = await fetch(target, {
        method: request.method,
        headers,
        body:
          request.method === "GET" || request.method === "HEAD" ? undefined : request.body,
        redirect: "manual",
      });
    } catch (err) {
      return json(
        {
          error:
            "Cannot reach the design server. It may be offline, or its tunnel may be reconnecting.",
          detail: String(err),
        },
        502,
      );
    }

    // Stream the body through untouched so per-gate progress still arrives live.
    const out = new Headers(upstream.headers);
    out.delete("content-encoding");
    out.delete("content-length");
    return new Response(upstream.body, { status: upstream.status, headers: out });
  },
};

function json(body, status) {
  return new Response(JSON.stringify(body), {
    status,
    headers: { "content-type": "application/json", "cache-control": "no-store" },
  });
}
