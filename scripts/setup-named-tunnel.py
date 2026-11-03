#!/usr/bin/env python3
"""Replace the quick tunnel with a permanent, named Cloudflare tunnel.

    scripts/setup-named-tunnel.py --dry-run          # show the plan; needs no token, touches nothing
    scripts/setup-named-tunnel.py                    # do it
    scripts/setup-named-tunnel.py --deploy           # ...and redeploy the site against the new hostname

Why: a quick tunnel gets a random hostname on every restart and Cloudflare rate-limits creating
them, so the site breaks whenever the tunnel restarts. A named tunnel keeps one hostname for good.
Like the quick tunnel it is OUTBOUND-ONLY: nothing listens on the shared DSI box.

What it does (every step is safe to repeat; a second run changes nothing):
  1. finds or creates the tunnel (remotely managed, so its routing lives in Cloudflare)
  2. points its ingress at the local web app, with a catch-all 404 for any other hostname
  3. finds or creates the proxied DNS record  <hostname> -> <tunnel id>.cfargotunnel.com
  4. fetches the tunnel's connector token and installs it on DSI (mode 600, over ssh stdin, so it
     is never on a command line or in a log), sets PEGRNA_PUBLIC_URL, restarts the tunnel service
  5. waits until https://<hostname>/api/health answers

The token it needs
------------------
A Cloudflare API token with:  Account > Cloudflare Tunnel > Edit   and   Zone > DNS > Edit
(for the zone that contains the hostname). Create it at dash.cloudflare.com > My Profile > API
Tokens > Create Token > Create Custom Token, then store it WITHOUT pasting it anywhere:

    read -rs T; security add-generic-password -a "$USER" -s cloudflare-tunnel-token -w "$T" -U; unset T

or export CLOUDFLARE_TUNNEL_TOKEN for a single run. The token is never printed.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

API = "https://api.cloudflare.com/client/v4"
KEYCHAIN_SERVICE = "cloudflare-tunnel-token"
HOSTNAME_RE = re.compile(r"^(?=.{4,253}$)([a-z0-9]([a-z0-9-]*[a-z0-9])?\.)+[a-z]{2,}$")
SITE_ORIGIN = "https://pegrna-designer.vercel.app"

TOKEN_HELP = (
    "No Cloudflare token found. Create one with these permissions:\n"
    "    Account > Cloudflare Tunnel > Edit\n"
    "    Zone    > DNS               > Edit   (for the zone containing the hostname)\n"
    "at dash.cloudflare.com > My Profile > API Tokens > Create Token > Create Custom Token, then\n"
    "store it without pasting it anywhere:\n"
    '    read -rs T; security add-generic-password -a "$USER" -s '
    f'{KEYCHAIN_SERVICE} -w "$T" -U; unset T\n'
    "(or export CLOUDFLARE_TUNNEL_TOKEN for one run)."
)


class SetupError(RuntimeError):
    """Something the operator has to fix. Never carries a secret."""


@dataclass
class Config:
    hostname: str
    tunnel_name: str = "pegrna-dsi"
    service: str = "http://127.0.0.1:8093"
    ssh_host: str = "dsi"
    account_id: Optional[str] = None
    zone_id: Optional[str] = None


@dataclass
class Provisioned:
    account_id: str
    zone_id: str
    tunnel_id: str
    connector_token: str
    actions: list


# ---- Cloudflare API ---------------------------------------------------------------------


class Cloudflare:
    def __init__(self, token: str, opener: Callable = urllib.request.urlopen) -> None:
        self._token = token
        self._open = opener

    def call(self, method: str, path: str, body=None):
        data = json.dumps(body).encode() if body is not None else None
        req = urllib.request.Request(
            API + path,
            data=data,
            method=method,
            headers={
                "Authorization": f"Bearer {self._token}",
                "Content-Type": "application/json",
            },
        )
        try:
            with self._open(req, timeout=60) as resp:
                payload = json.loads(resp.read().decode())
        except urllib.error.HTTPError as exc:
            try:
                payload = json.loads(exc.read().decode())
            except (ValueError, OSError):
                payload = {"success": False, "errors": [{"message": f"HTTP {exc.code}"}]}
            payload.setdefault("_status", exc.code)
        if not payload.get("success", False):
            msgs = "; ".join(e.get("message", "?") for e in payload.get("errors", [])) or "failed"
            raise SetupError(_explain(method, path, msgs, payload.get("_status")))
        return payload.get("result")


def _explain(method: str, path: str, msgs: str, status) -> str:
    base = f"Cloudflare rejected {method} {path}: {msgs}"
    low = msgs.lower()
    if (
        status in (401, 403)
        or "authentication" in low
        or "permission" in low
        or "not authorized" in low
    ):
        need = (
            "Zone > DNS > Edit"
            if "/dns_records" in path or path.startswith("/zones")
            else ("Account > Cloudflare Tunnel > Edit")
        )
        return f"{base}\n  This token lacks a permission it needs here: {need}."
    return base


def read_token(
    env: dict = os.environ, run: Callable = subprocess.run, service: str = KEYCHAIN_SERVICE
) -> str:
    """The Cloudflare token: CLOUDFLARE_TUNNEL_TOKEN, else the macOS keychain item."""
    if env.get("CLOUDFLARE_TUNNEL_TOKEN", "").strip():
        return env["CLOUDFLARE_TUNNEL_TOKEN"].strip()
    try:
        out = run(
            ["security", "find-generic-password", "-s", service, "-w"],
            capture_output=True,
            text=True,
            check=False,
        )
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip()
    except FileNotFoundError:
        pass  # no `security` binary: not macOS
    raise SetupError(TOKEN_HELP)


# ---- provisioning -----------------------------------------------------------------------


def _discover_account(api: Cloudflare) -> str:
    try:
        accounts = api.call("GET", "/accounts?per_page=50") or []
    except SetupError as exc:
        raise SetupError(
            f"{exc}\n  Pass --account-id instead (it is shown on any zone's Overview page in the "
            "Cloudflare dashboard)."
        ) from exc
    if len(accounts) == 1:
        return accounts[0]["id"]
    names = ", ".join(f"{a.get('name')} ({a['id']})" for a in accounts) or "none visible"
    raise SetupError(f"Cannot choose an account automatically ({names}); pass --account-id.")


def _discover_zone(api: Cloudflare, hostname: str) -> str:
    labels = hostname.split(".")
    for i in range(1, len(labels) - 1):
        name = ".".join(labels[i:])
        zones = api.call("GET", f"/zones?name={name}") or []
        if zones:
            return zones[0]["id"]
    raise SetupError(
        f"No Cloudflare zone found for {hostname}. The domain must be on this Cloudflare "
        "account, and the token must be allowed to read it."
    )


def provision(api: Cloudflare, cfg: Config) -> Provisioned:
    """Make Cloudflare's side match the config. Idempotent: reads first, writes only a gap."""
    actions = []
    account = cfg.account_id or _discover_account(api)
    zone = cfg.zone_id or _discover_zone(api, cfg.hostname)

    found = api.call(
        "GET", f"/accounts/{account}/cfd_tunnel?name={cfg.tunnel_name}&is_deleted=false"
    )
    if found:
        tunnel = found[0]
        if tunnel.get("config_src") not in (None, "cloudflare"):
            raise SetupError(
                f"A tunnel named {cfg.tunnel_name!r} already exists but is managed by a local "
                "config file, so its routing cannot be set from here. Delete it in the Cloudflare "
                "dashboard or pass --tunnel-name."
            )
        actions.append(f"tunnel {cfg.tunnel_name!r} already exists ({tunnel['id']})")
    else:
        tunnel = api.call(
            "POST",
            f"/accounts/{account}/cfd_tunnel",
            {"name": cfg.tunnel_name, "config_src": "cloudflare"},
        )
        actions.append(f"created tunnel {cfg.tunnel_name!r} ({tunnel['id']})")
    tid = tunnel["id"]

    api.call(
        "PUT",
        f"/accounts/{account}/cfd_tunnel/{tid}/configurations",
        {
            "config": {
                "ingress": [
                    {"hostname": cfg.hostname, "service": cfg.service},
                    {"service": "http_status:404"},
                ]
            }
        },
    )
    actions.append(f"ingress: {cfg.hostname} -> {cfg.service} (everything else 404)")

    want = f"{tid}.cfargotunnel.com"
    records = api.call("GET", f"/zones/{zone}/dns_records?name={cfg.hostname}") or []
    if records:
        rec = records[0]
        if rec.get("type") != "CNAME":
            raise SetupError(
                f"{cfg.hostname} already has a {rec.get('type')} record. It is not removed "
                "automatically; delete it in the Cloudflare dashboard or pick another hostname."
            )
        if rec.get("content") == want and rec.get("proxied"):
            actions.append(f"DNS {cfg.hostname} already points at the tunnel")
        else:
            api.call(
                "PATCH",
                f"/zones/{zone}/dns_records/{rec['id']}",
                {"type": "CNAME", "name": cfg.hostname, "content": want, "proxied": True},
            )
            actions.append(f"DNS {cfg.hostname} repointed at the tunnel")
    else:
        api.call(
            "POST",
            f"/zones/{zone}/dns_records",
            {
                "type": "CNAME",
                "name": cfg.hostname,
                "content": want,
                "proxied": True,
                "comment": "pegRNA backend (scripts/setup-named-tunnel.py)",
            },
        )
        actions.append(f"DNS created: {cfg.hostname} -> {want} (proxied)")

    token = api.call("GET", f"/accounts/{account}/cfd_tunnel/{tid}/token")
    if not isinstance(token, str) or not token:
        raise SetupError("Cloudflare returned no connector token for the tunnel.")
    return Provisioned(account, zone, tid, token, actions)


def plan(cfg: Config) -> list:
    """The steps, for --dry-run. Reads as the script would run; contains no secret."""
    return [
        f"find or create tunnel {cfg.tunnel_name!r} (remotely managed)  [needs Account > Cloudflare Tunnel > Edit]",
        f"set its ingress: {cfg.hostname} -> {cfg.service}; every other hostname -> 404",
        f"find or create proxied CNAME {cfg.hostname} -> <tunnel id>.cfargotunnel.com  [needs Zone > DNS > Edit]",
        "fetch the tunnel's connector token (held in memory only)",
        f"on {cfg.ssh_host}: write it to ~/.config/pegrna/tunnel.token (mode 600), over ssh stdin",
        f"on {cfg.ssh_host}: set PEGRNA_PUBLIC_URL=https://{cfg.hostname} in ~/.config/pegrna/pegrna.env",
        f"on {cfg.ssh_host}: restart pegrna-tunnel; the unit runs web/deploy/tunnel-run.sh, which uses the token",
        f"wait for https://{cfg.hostname}/api/health",
    ]


# ---- installing on DSI ------------------------------------------------------------------


def _ssh(cfg: Config, command: str, *, stdin: Optional[str] = None, run: Callable = subprocess.run):
    out = run(
        ["ssh", "-o", "ConnectTimeout=20", cfg.ssh_host, command],
        input=stdin,
        capture_output=True,
        text=True,
        check=False,
    )
    if out.returncode != 0:
        raise SetupError(
            f"ssh {cfg.ssh_host} failed ({out.returncode}): {(out.stderr or '').strip()[-300:]}\n"
            "Is Zscaler (ZPA) running?"
        )
    return out.stdout


def install_on_dsi(cfg: Config, connector_token: str, run: Callable = subprocess.run) -> None:
    url = f"https://{cfg.hostname}"
    if not HOSTNAME_RE.match(cfg.hostname):  # it is interpolated into a remote shell command
        raise SetupError(f"Not a valid hostname: {cfg.hostname!r}")

    # The unit only uses the token if tunnel-run.sh knows about it; refuse to switch over to a
    # script that would ignore it and quietly keep running a quick tunnel.
    _ssh(
        cfg,
        "grep -q PEGRNA_TUNNEL_TOKEN_FILE ~/pegrna/prime-editing/web/deploy/tunnel-run.sh "
        "|| { echo 'the DSI copy of tunnel-run.sh is out of date: rsync the repo first' >&2; exit 3; }",
        run=run,
    )
    # The token travels on stdin: not in argv (visible to `ps`), not in shell history, not in logs.
    _ssh(
        cfg,
        "umask 077; mkdir -p ~/.config/pegrna && cat > ~/.config/pegrna/tunnel.token.new "
        "&& mv ~/.config/pegrna/tunnel.token.new ~/.config/pegrna/tunnel.token",
        stdin=connector_token + "\n",
        run=run,
    )
    _ssh(
        cfg,
        "f=~/.config/pegrna/pegrna.env; "
        f"if grep -q '^PEGRNA_PUBLIC_URL=' \"$f\"; then sed -i 's|^PEGRNA_PUBLIC_URL=.*|PEGRNA_PUBLIC_URL={url}|' \"$f\"; "
        f"else printf '\\nPEGRNA_PUBLIC_URL={url}\\n' >> \"$f\"; fi; "
        "systemctl --user restart pegrna-tunnel",
        run=run,
    )


def wait_healthy(
    url: str,
    *,
    attempts: int = 24,
    delay: float = 5.0,
    opener: Callable = urllib.request.urlopen,
    sleep: Callable = time.sleep,
) -> bool:
    for i in range(attempts):
        try:
            with opener(url, timeout=15) as resp:
                if resp.status == 200 and json.loads(resp.read().decode()).get("status") == "ok":
                    return True
        except (urllib.error.URLError, ValueError, OSError):
            pass
        if i < attempts - 1:
            sleep(delay)
    return False


# ---- entry point ------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--hostname", default="pegrna.ukdri.org")
    ap.add_argument("--tunnel-name", default="pegrna-dsi")
    ap.add_argument(
        "--service", default="http://127.0.0.1:8093", help="what the tunnel forwards to"
    )
    ap.add_argument("--ssh-host", default="dsi")
    ap.add_argument("--account-id")
    ap.add_argument("--zone-id")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-install", action="store_true", help="configure Cloudflare only")
    ap.add_argument("--deploy", action="store_true", help="then redeploy the site to Vercel")
    args = ap.parse_args(argv)

    cfg = Config(
        hostname=args.hostname.lower(),
        tunnel_name=args.tunnel_name,
        service=args.service,
        ssh_host=args.ssh_host,
        account_id=args.account_id,
        zone_id=args.zone_id,
    )
    if not HOSTNAME_RE.match(cfg.hostname):
        print(f"error: not a valid hostname: {cfg.hostname!r}", file=sys.stderr)
        return 2

    if args.dry_run:
        print(
            f"DRY RUN: nothing is read from the keychain or changed.\n\nPlan for https://{cfg.hostname}:"
        )
        for i, step in enumerate(plan(cfg), 1):
            print(f"  {i}. {step}")
        print(f"\nThen:  PEGRNA_API_BASE=https://{cfg.hostname} bash web/vercel/deploy.sh --prod")
        return 0

    try:
        api = Cloudflare(read_token())
        done = provision(api, cfg)
        for line in done.actions:
            print("  cloudflare:", line)
        if args.no_install:
            print("Cloudflare is configured; skipped the DSI install (--no-install).")
            return 0
        install_on_dsi(cfg, done.connector_token)
        print(
            f"  dsi: token installed, pegrna-tunnel restarted; waiting for https://{cfg.hostname} ..."
        )
        if not wait_healthy(f"https://{cfg.hostname}/api/health"):
            raise SetupError(
                f"https://{cfg.hostname}/api/health did not answer within two minutes. DNS can take "
                "a little longer: re-run this script (it changes nothing that is already right)."
            )
        print(f"  live: https://{cfg.hostname}/api/health is answering")
    except SetupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    deploy = f"PEGRNA_API_BASE=https://{cfg.hostname} bash web/vercel/deploy.sh --prod"
    if args.deploy:
        root = Path(__file__).resolve().parent.parent
        rc = subprocess.run(
            ["bash", "web/vercel/deploy.sh", "--prod"],
            cwd=root,
            env={**os.environ, "PEGRNA_API_BASE": f"https://{cfg.hostname}"},
        ).returncode
        return rc
    print(f"\nNext: {deploy}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
