"""scripts/setup-named-tunnel.py: provisioning logic, idempotence, and secret handling.

The script is the last step from a throwaway tunnel hostname to a permanent one, and it handles
two credentials (the Cloudflare API token and the tunnel's connector token). The things worth
pinning are that a second run changes nothing, that a conflict is reported instead of overwritten,
and that neither secret ever reaches output, a log, or a command line.
"""

import importlib.util
import io
import json
import urllib.error
from pathlib import Path
from types import SimpleNamespace

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "setup_named_tunnel", Path(__file__).resolve().parents[1] / "scripts" / "setup-named-tunnel.py"
)
snt = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(snt)

API_TOKEN = "cf-api-token-SECRET-1234"
CONNECTOR = "connector-token-SECRET-5678"
HOST = "pegrna.example.org"


class FakeAPI:
    """Just enough of Cloudflare's v4 API, with state, recording every call."""

    def __init__(
        self, *, tunnels=None, records=None, accounts=None, zones=None, zone_account="acc1"
    ):
        self.calls = []
        self.accounts = accounts if accounts is not None else [{"id": "acc1", "name": "Org"}]
        self.zones = zones if zones is not None else {"example.org": "zone1"}
        self.zone_account = zone_account  # None: a zone that does not say who owns it
        self.tunnels = list(tunnels or [])
        self.records = list(records or [])

    def call(self, method, path, body=None):
        self.calls.append((method, path, body))
        if path.startswith("/accounts?"):
            return self.accounts
        if path.startswith("/zones?name="):
            z = self.zones.get(path.split("=", 1)[1])
            return [self._zone(z)] if z else []
        if method == "GET" and path.startswith("/zones/") and path.count("/") == 2:
            zid = path.rsplit("/", 1)[1]
            return self._zone(zid) if zid in self.zones.values() else None
        if "/cfd_tunnel?name=" in path:
            return [t for t in self.tunnels if t["name"] == path.split("name=")[1].split("&")[0]]
        if method == "POST" and path.endswith("/cfd_tunnel"):
            t = {"id": "tun1", "name": body["name"], "config_src": body["config_src"]}
            self.tunnels.append(t)
            return t
        if method == "PUT" and path.endswith("/configurations"):
            return {}
        if path.endswith("/token"):
            return CONNECTOR
        if "/dns_records?name=" in path:
            return [r for r in self.records if r["name"] == path.split("=", 1)[1]]
        if method == "POST" and path.endswith("/dns_records"):
            self.records.append({"id": "rec1", **body})
            return {"id": "rec1"}
        if method == "PATCH" and "/dns_records/" in path:
            return {}
        raise AssertionError(f"unexpected call {method} {path}")

    def _zone(self, zid):
        zone = {"id": zid}
        if self.zone_account:
            zone["account"] = {"id": self.zone_account, "name": "Org"}
        return zone

    def writes(self):
        return [c for c in self.calls if c[0] in ("POST", "PATCH", "DELETE")]


def cfg(**kw):
    return snt.Config(hostname=HOST, **kw)


# ---- provisioning -----------------------------------------------------------------------


def test_a_fresh_account_gets_a_tunnel_ingress_dns_and_a_token():
    api = FakeAPI()
    done = snt.provision(api, cfg())
    assert done.tunnel_id == "tun1" and done.connector_token == CONNECTOR
    create = next(c for c in api.calls if c[0] == "POST" and c[1].endswith("/cfd_tunnel"))
    assert create[2] == {
        "name": "pegrna-dsi",
        "config_src": "cloudflare",
    }, "remotely managed, or its ingress could not be set through the API"
    put = next(c for c in api.calls if c[0] == "PUT")
    assert put[2]["config"]["ingress"] == [
        {"hostname": HOST, "service": "http://127.0.0.1:8093"},
        {"service": "http_status:404"},
    ], "anything but our hostname must 404"
    dns = next(c for c in api.calls if c[0] == "POST" and c[1].endswith("/dns_records"))
    assert dns[2]["type"] == "CNAME" and dns[2]["proxied"] is True
    assert dns[2]["content"] == "tun1.cfargotunnel.com"


def test_a_second_run_creates_nothing():
    first = FakeAPI()
    snt.provision(first, cfg())
    again = FakeAPI(tunnels=first.tunnels, records=first.records)
    snt.provision(again, cfg())
    assert again.writes() == [], "nothing already correct may be rewritten"


def test_a_stale_dns_record_is_repointed_not_duplicated():
    stale = {
        "id": "r9",
        "type": "CNAME",
        "name": HOST,
        "content": "old.cfargotunnel.com",
        "proxied": True,
    }
    api = FakeAPI(
        tunnels=[{"id": "tun1", "name": "pegrna-dsi", "config_src": "cloudflare"}], records=[stale]
    )
    snt.provision(api, cfg())
    (patch,) = [c for c in api.calls if c[0] == "PATCH"]
    assert patch[1].endswith("/dns_records/r9") and patch[2]["content"] == "tun1.cfargotunnel.com"
    assert not [c for c in api.calls if c[0] == "POST" and c[1].endswith("/dns_records")]


def test_an_unproxied_record_is_fixed():
    rec = {
        "id": "r1",
        "type": "CNAME",
        "name": HOST,
        "content": "tun1.cfargotunnel.com",
        "proxied": False,
    }
    api = FakeAPI(
        tunnels=[{"id": "tun1", "name": "pegrna-dsi", "config_src": "cloudflare"}], records=[rec]
    )
    snt.provision(api, cfg())
    assert any(c[0] == "PATCH" and c[2]["proxied"] is True for c in api.calls)


def test_a_conflicting_record_type_is_reported_and_left_alone():
    api = FakeAPI(records=[{"id": "r1", "type": "A", "name": HOST, "content": "1.2.3.4"}])
    with pytest.raises(snt.SetupError, match="already has a A record"):
        snt.provision(api, cfg())
    assert not [
        c for c in api.calls if c[0] in ("PATCH", "DELETE")
    ], "never delete someone's record"


def test_a_locally_managed_tunnel_of_the_same_name_is_refused():
    api = FakeAPI(tunnels=[{"id": "t0", "name": "pegrna-dsi", "config_src": "local"}])
    with pytest.raises(snt.SetupError, match="managed by a local config"):
        snt.provision(api, cfg())


def test_the_zone_is_found_from_the_hostname():
    api = FakeAPI(zones={"example.org": "zoneX"})
    assert snt.provision(api, cfg()).zone_id == "zoneX"


def test_an_unknown_domain_is_explained():
    with pytest.raises(snt.SetupError, match="No Cloudflare zone found"):
        snt.provision(FakeAPI(zones={}), cfg())


class NoAccountListing(FakeAPI):
    """A token scoped to one tunnel and one zone: it can read the zone but not list accounts."""

    def call(self, method, path, body=None):
        if path.startswith("/accounts?"):
            raise snt.SetupError("Cloudflare rejected GET /accounts: forbidden")
        return super().call(method, path, body)


def test_a_token_that_cannot_list_accounts_uses_the_zones_account():
    """The first real run: this is exactly what the token could and could not do."""
    done = snt.provision(NoAccountListing(), cfg())
    assert done.account_id == "acc1" and done.tunnel_id == "tun1"


def test_the_zones_account_is_used_without_listing_accounts():
    api = FakeAPI(
        accounts=[{"id": "a", "name": "One"}, {"id": "b", "name": "Two"}], zone_account="b"
    )
    assert snt.provision(api, cfg()).account_id == "b"
    assert not [c for c in api.calls if c[1].startswith("/accounts?")], "no need to list them"


def test_an_explicit_account_id_still_wins():
    assert snt.provision(FakeAPI(), cfg(account_id="other")).account_id == "other"


def test_a_given_zone_id_is_read_to_find_its_account():
    api = FakeAPI(zones={"example.org": "zoneX"}, zone_account="accZ")
    done = snt.provision(api, cfg(zone_id="zoneX"))
    assert done.zone_id == "zoneX" and done.account_id == "accZ"
    assert ("GET", "/zones/zoneX", None) in api.calls


def test_a_zone_id_that_cannot_be_read_falls_back_to_the_account_list():
    class NoZoneRead(FakeAPI):
        def call(self, method, path, body=None):
            if path == "/zones/zoneX":
                raise snt.SetupError("Cloudflare rejected GET /zones/zoneX: forbidden")
            return super().call(method, path, body)

    api = NoZoneRead(zones={"example.org": "zoneX"})
    assert snt.provision(api, cfg(zone_id="zoneX")).account_id == "acc1"  # the one listed account


def test_when_nothing_names_an_account_the_id_is_requested():
    with pytest.raises(snt.SetupError, match="--account-id"):
        snt.provision(NoAccountListing(zone_account=None), cfg())
    done = snt.provision(NoAccountListing(zone_account=None), cfg(account_id="acc1"))
    assert done.tunnel_id == "tun1"


def test_several_accounts_and_an_unknown_zone_owner_require_an_explicit_choice():
    api = FakeAPI(
        accounts=[{"id": "a", "name": "One"}, {"id": "b", "name": "Two"}], zone_account=None
    )
    with pytest.raises(snt.SetupError, match="--account-id"):
        snt.provision(api, cfg())
    assert snt.provision(FakeAPI(zone_account=None), cfg(account_id="acc1")).account_id == "acc1"


# ---- the API client ---------------------------------------------------------------------


def _http_error(code, payload):
    return urllib.error.HTTPError("u", code, "x", {}, io.BytesIO(json.dumps(payload).encode()))


def test_a_missing_permission_names_the_permission_and_not_the_token():
    def opener(req, timeout=None):
        raise _http_error(403, {"success": False, "errors": [{"message": "Authentication error"}]})

    with pytest.raises(snt.SetupError) as exc:
        snt.Cloudflare(API_TOKEN, opener).call("GET", "/accounts/a/cfd_tunnel")
    assert "Cloudflare Tunnel > Edit" in str(exc.value)
    assert API_TOKEN not in str(exc.value)


def test_a_dns_permission_failure_says_dns():
    def opener(req, timeout=None):
        raise _http_error(403, {"success": False, "errors": [{"message": "forbidden"}]})

    with pytest.raises(snt.SetupError) as exc:
        snt.Cloudflare(API_TOKEN, opener).call("POST", "/zones/z/dns_records", {})
    assert "DNS" in str(exc.value)


def test_the_token_is_sent_as_a_bearer_credential():
    seen = {}

    class R:
        def read(self):
            return b'{"success": true, "result": {"ok": 1}}'

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def opener(req, timeout=None):
        seen["auth"] = req.get_header("Authorization")
        return R()

    assert snt.Cloudflare(API_TOKEN, opener).call("GET", "/x") == {"ok": 1}
    assert seen["auth"] == f"Bearer {API_TOKEN}"


# ---- reading the token ------------------------------------------------------------------


def test_the_environment_variable_wins():
    assert snt.read_token({"CLOUDFLARE_TUNNEL_TOKEN": f" {API_TOKEN} "}) == API_TOKEN


def test_the_keychain_is_used_otherwise():
    run = lambda *a, **k: SimpleNamespace(returncode=0, stdout=API_TOKEN + "\n")  # noqa: E731
    assert snt.read_token({}, run) == API_TOKEN


def test_no_token_gives_the_exact_instructions():
    run = lambda *a, **k: SimpleNamespace(returncode=44, stdout="")  # noqa: E731
    with pytest.raises(snt.SetupError) as exc:
        snt.read_token({}, run)
    msg = str(exc.value)
    assert "Cloudflare Tunnel > Edit" in msg and "security add-generic-password" in msg
    assert snt.KEYCHAIN_SERVICE in msg
    assert "read -rs" in msg, "must not tell people to paste the secret into a command line"


def test_a_host_without_the_security_binary_still_gets_instructions():
    def run(*a, **k):
        raise FileNotFoundError("security")

    with pytest.raises(snt.SetupError, match="No Cloudflare token"):
        snt.read_token({}, run)


# ---- installing on DSI ------------------------------------------------------------------


class FakeSSH:
    def __init__(self, fail_on=None):
        self.calls = []
        self.fail_on = fail_on

    def __call__(self, argv, input=None, **kw):  # noqa: A002
        self.calls.append((argv, input))
        rc = 3 if self.fail_on and self.fail_on in argv[-1] else 0
        return SimpleNamespace(returncode=rc, stdout="", stderr="out of date" if rc else "")


def test_the_connector_token_travels_on_stdin_never_in_argv():
    ssh = FakeSSH()
    snt.install_on_dsi(cfg(), CONNECTOR, run=ssh)
    assert [i for _, i in ssh.calls if i and CONNECTOR in i], "token must be sent via stdin"
    for argv, _ in ssh.calls:
        assert CONNECTOR not in " ".join(argv), "a token in argv shows up in `ps`"
    install = next(a[-1] for a, i in ssh.calls if i)
    assert "umask 077" in install and "tunnel.token" in install


def test_the_remote_env_gets_the_public_url_and_the_service_restarts():
    ssh = FakeSSH()
    snt.install_on_dsi(cfg(), CONNECTOR, run=ssh)
    last = ssh.calls[-1][0][-1]
    assert f"PEGRNA_PUBLIC_URL=https://{HOST}" in last
    assert "systemctl --user restart pegrna-tunnel" in last


def test_a_stale_remote_tunnel_script_stops_the_install_before_anything_changes():
    ssh = FakeSSH(fail_on="PEGRNA_TUNNEL_TOKEN_FILE")
    with pytest.raises(snt.SetupError, match="Zscaler"):
        snt.install_on_dsi(cfg(), CONNECTOR, run=ssh)
    assert len(ssh.calls) == 1, "the token must not be written if the script would ignore it"


def test_a_hostname_that_could_inject_into_the_remote_shell_is_refused():
    for bad in ("x.org; rm -rf ~", "a b.org", "$(id).org", "x.org'"):
        with pytest.raises(snt.SetupError, match="valid hostname"):
            snt.install_on_dsi(snt.Config(hostname=bad), CONNECTOR, run=FakeSSH())


# ---- waiting for health -----------------------------------------------------------------


class _Health:
    def __init__(self, body):
        self.status = 200
        self._b = body

    def read(self):
        return self._b

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_health_is_polled_until_it_answers():
    outcomes = [
        urllib.error.URLError("dns"),
        urllib.error.URLError("dns"),
        _Health(b'{"status": "ok"}'),
    ]

    def opener(url, timeout=None):
        o = outcomes.pop(0)
        if isinstance(o, Exception):
            raise o
        return o

    sleeps = []
    assert snt.wait_healthy("https://x/api/health", opener=opener, sleep=sleeps.append) is True
    assert len(sleeps) == 2


def test_health_gives_up_after_the_attempts():
    def opener(url, timeout=None):
        raise urllib.error.URLError("no route")

    assert snt.wait_healthy("https://x", attempts=3, opener=opener, sleep=lambda s: None) is False


def test_the_default_wait_covers_a_first_time_setup():
    """The first real run needed a little over two minutes; the old default gave up at two."""
    assert snt.HEALTH_WAIT_ATTEMPTS * snt.HEALTH_WAIT_DELAY >= 300

    def opener(url, timeout=None):
        raise urllib.error.URLError("not yet")

    sleeps = []
    assert snt.wait_healthy("https://x", opener=opener, sleep=sleeps.append) is False
    assert sum(sleeps) >= 290, "the default must keep trying for about five minutes"


def test_a_non_ok_status_is_not_healthy():
    opener = lambda url, timeout=None: _Health(b'{"status": "degraded"}')  # noqa: E731
    assert snt.wait_healthy("https://x", attempts=2, opener=opener, sleep=lambda s: None) is False


# ---- the command ------------------------------------------------------------------------


def test_dry_run_reads_no_secret_and_makes_no_call(monkeypatch, capsys):
    def boom(*a, **k):
        raise AssertionError("dry run must not touch the keychain or the network")

    monkeypatch.setattr(snt, "read_token", boom)
    monkeypatch.setattr(snt, "Cloudflare", boom)
    assert snt.main(["--hostname", HOST, "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "DRY RUN" in out and HOST in out and "web/vercel/deploy.sh" in out


def test_full_run_never_prints_either_secret(monkeypatch, capsys):
    api = FakeAPI()
    monkeypatch.setattr(snt, "read_token", lambda *a, **k: API_TOKEN)
    monkeypatch.setattr(snt, "Cloudflare", lambda token: api)
    installed = []
    monkeypatch.setattr(snt, "install_on_dsi", lambda c, t, **k: installed.append(t))
    monkeypatch.setattr(snt, "wait_healthy", lambda *a, **k: True)
    assert snt.main(["--hostname", HOST]) == 0
    captured = capsys.readouterr()
    for secret in (API_TOKEN, CONNECTOR):
        assert secret not in captured.out + captured.err
    assert installed == [CONNECTOR]
    assert f"PEGRNA_API_BASE=https://{HOST}" in captured.out


def test_a_failure_is_an_error_message_and_a_nonzero_exit(monkeypatch, capsys):
    monkeypatch.setattr(
        snt, "read_token", lambda *a, **k: (_ for _ in ()).throw(snt.SetupError("no token"))
    )
    assert snt.main(["--hostname", HOST]) == 1
    assert "error: no token" in capsys.readouterr().err


def test_an_invalid_hostname_is_refused_before_anything_happens(capsys):
    assert snt.main(["--hostname", "not a host"]) == 2


def test_no_install_stops_after_cloudflare(monkeypatch):
    monkeypatch.setattr(snt, "read_token", lambda *a, **k: API_TOKEN)
    monkeypatch.setattr(snt, "Cloudflare", lambda token: FakeAPI())
    monkeypatch.setattr(snt, "install_on_dsi", lambda *a, **k: pytest.fail("must not install"))
    assert snt.main(["--hostname", HOST, "--no-install"]) == 0


def test_a_health_timeout_says_what_is_already_done_and_what_to_do_next(monkeypatch, capsys):
    """Nothing is wrong when the URL is slow to answer, and re-running restarts the tunnel again,
    so the message must say so instead of telling the user to re-run."""
    monkeypatch.setattr(snt, "read_token", lambda *a, **k: API_TOKEN)
    monkeypatch.setattr(snt, "Cloudflare", lambda token: FakeAPI())
    monkeypatch.setattr(snt, "install_on_dsi", lambda *a, **k: None)
    monkeypatch.setattr(snt, "wait_healthy", lambda *a, **k: False)
    assert snt.main(["--hostname", HOST, "--deploy"]) == 1
    err = capsys.readouterr().err
    assert "ARE configured" in err and "cloudflared.log" in err
    assert f"PEGRNA_API_BASE=https://{HOST} bash web/vercel/deploy.sh --prod" in err
    assert "Do not re-run" in err
    for secret in (API_TOKEN, CONNECTOR):
        assert secret not in err
