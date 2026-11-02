"""`fake` mode is documented as needing no external service, and the web app
offers it as the mode that "always works". It has to actually hold.

The auto-spacer path widens a short genome window via `resolve_pridict_context`,
which previously fell through to a live Ensembl fetch -- so a fake-mode design
with a window even a few bp short would hang on the network and then fail. (Every
shipped example was 4 bp under `min_flank=100`, so this was the default path, not
an edge case.) In fake mode PRIDICT is itself a fake, so the context is padded
locally instead.
"""

import json
from pathlib import Path

import pytest

from pegrna_pipeline.hpc.recoded_driver import design_recoded


@pytest.fixture(autouse=True)
def _ban_network(monkeypatch):
    """Make any outbound HTTP attempt an immediate, obvious failure."""
    import socket

    def _blocked(*args, **kwargs):
        raise AssertionError("fake mode attempted a network connection")

    monkeypatch.setattr(socket.socket, "connect", _blocked)
    monkeypatch.setattr(socket, "create_connection", _blocked)


def _cfg(name: str) -> dict:
    return json.loads(Path(f"examples/ctcf_cterm/{name}").read_text())


def test_autospacer_example_designs_offline():
    cfg = _cfg("config_autospacer.json")
    out = design_recoded(
        rtt_aa=cfg["rtt_aa"],
        cargo=cfg["cargo"],
        genome=cfg["genome"],
        site_s=cfg["site_s"],
        cds_mask=cfg["cds_mask"],
        pbs_lengths=cfg["pbs_lengths"],
        ngrna_distance_min=cfg["ngrna_distance_min"],
        ngrna_distance_max=cfg["ngrna_distance_max"],
        pegrna_strand=cfg["pegrna_strand"],
        arm_side=cfg["arm_side"],
        max_run=cfg["max_run"],
        transcript=cfg.get("transcript"),
        mode="fake",
    )
    assert out.result.spacer
    assert len(out.result.spacer) == 20


def test_short_window_without_transcript_still_designs_offline():
    """A user-supplied window under min_flank must not reach for Ensembl."""
    cfg = _cfg("config_autospacer.json")
    genome = cfg["genome"]
    site_s = cfg["site_s"]
    # Trim the right flank well under min_flank (100) and drop the transcript.
    keep = site_s + 40
    out = design_recoded(
        rtt_aa=cfg["rtt_aa"],
        cargo=cfg["cargo"],
        genome=genome[:keep],
        site_s=site_s,
        cds_mask=cfg["cds_mask"][:keep],
        pbs_lengths=cfg["pbs_lengths"],
        ngrna_distance_min=cfg["ngrna_distance_min"],
        ngrna_distance_max=cfg["ngrna_distance_max"],
        pegrna_strand=cfg["pegrna_strand"],
        arm_side=cfg["arm_side"],
        max_run=cfg["max_run"],
        transcript=None,
        mode="fake",
    )
    assert len(out.result.spacer) == 20


def test_shipped_examples_have_enough_flank_to_skip_any_fetch():
    """The examples should clear min_flank outright, so even api mode runs
    without a network round-trip for context."""
    from pegrna_pipeline.hpc.recoded_driver import design_recoded as _d  # noqa: F401

    min_flank = 100
    for name in ("config.json", "config_autospacer.json", "config_autospacer_flag.json"):
        cfg = _cfg(name)
        left = cfg["site_s"]
        right = len(cfg["genome"]) - cfg["site_s"]
        assert left >= min_flank, f"{name}: left flank {left} < {min_flank}"
        assert right >= min_flank, f"{name}: right flank {right} < {min_flank}"
        assert len(cfg["cds_mask"]) == len(cfg["genome"]), f"{name}: cds_mask length mismatch"
