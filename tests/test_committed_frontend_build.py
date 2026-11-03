"""The COMMITTED frontend build must be the one FastAPI serves.

web/static/ is what `web/app.py` serves at `/`, from the `/static` mount, calling the API on
the same origin. A build for another host (Vercel) uses a different asset base (`/`) and
bakes in that deployment's backend URL. Such a bundle once ended up committed in its place:
the repo's own copy of the UI pointed at `/assets/...` (which FastAPI does not serve) and at
a throwaway tunnel hostname, so it loaded blank. Nothing noticed, because the Vercel site
built from the same working tree kept working.
"""

import re
from pathlib import Path

STATIC = Path(__file__).resolve().parents[1] / "web" / "static"


def _html():
    return (STATIC / "index.html").read_text()


def test_assets_are_referenced_under_the_static_mount():
    refs = re.findall(r'(?:src|href)="(/[^"]+\.(?:js|css))"', _html())
    assert refs, "index.html references no bundle"
    for ref in refs:
        assert ref.startswith("/static/assets/"), (
            f"{ref}: the committed build must use FastAPI's /static base. A root-based bundle "
            "is the Vercel build; set VITE_OUT_DIR so it cannot overwrite web/static."
        )


def test_every_referenced_asset_exists():
    for ref in re.findall(r'(?:src|href)="/static/(assets/[^"]+)"', _html()):
        assert (STATIC / ref).is_file(), f"index.html references a missing file: {ref}"


def test_no_deployment_specific_backend_url_is_baked_into_the_bundle():
    """The default build talks to its own origin. A hardcoded host means it was built for
    somewhere else."""
    bad = re.compile(r"https://[a-z0-9-]+\.(?:trycloudflare\.com|lhr\.life|vercel\.app)")
    for js in (STATIC / "assets").glob("*.js"):
        hits = sorted(set(bad.findall(js.read_text())))
        assert not hits, f"{js.name} has a deployment URL baked in: {hits}"
