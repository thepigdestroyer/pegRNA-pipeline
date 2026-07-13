"""Auto-spacer wiring for ``design_recoded`` (Stage 2a of the 2026-07-10 PRIDICT2
wiring spec).

Previously the pipeline required a hand-supplied ``spacer``; it only auto-enumerated
the ngRNA nick / PBS / RTT / recode-arm around it. These tests pin the additive,
backward-compatible wiring:

* AUTO path — with no ``spacer`` and ``auto_spacer=True`` (the default), the pipeline
  builds a pegRNA#1 insertion editseq, calls the injected PRIDICT runner, and uses the
  top pegRNA's spacer. The runner is *injected* (from the scorer bundle in fake mode, or
  an explicit ``pridict_runner=`` kwarg) so CI never launches real PRIDICT2.
* OVERRIDE path — an explicit ``spacer`` uses the historical fixed-spacer path unchanged
  and never touches the runner.
* FAIL-LOUD — auto mode with no runner, or a runner that returns no pegRNA, raises a
  clear error (never silently degrades).

Everything downstream of the spacer (nick enumeration, PBS/RTT, recode, gate cascade) is
unchanged; the design runs end-to-end in ``mode="fake"`` on the deterministic fakes.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc import recoded_driver, scorer_factory
from pegrna_pipeline.hpc.data_types import PegRNADesign

# ---- shared inputs / helpers -------------------------------------------------

#: The fake bundle's PRIDICT stand-in (:class:`PridictRunnerFake`) emits this
#: pegRNA#1 spacer for every row, so the auto path in fake mode resolves to it.
_FAKE_BUNDLE_SPACER = "GAATTCGTGTTCCTAAGCGT"


def _auto_inputs() -> dict:
    """Design inputs WITHOUT a spacer (mirrors the known-good fake-mode fixture in
    ``test_recoded_driver.py``, spacer removed)."""
    return dict(
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome="ACGATGCATGGCATGACTGA" * 16,  # 320 bp -> flanks 150 / 170 at site_s=150
        site_s=150,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
    )


def _peg(spacer: str, *, pe2_k562: float = 0.5, pe2_hek293t: float = 0.5) -> PegRNADesign:
    """A minimal, gate-neutral :class:`PegRNADesign` carrying ``spacer``."""
    return PegRNADesign(
        sequence_name="peg1_p1",
        spacer=spacer,
        pbs="TCGACGCTTAGGAACAC",
        rt="ACGTACGTACGTACGTACGTACGTACGTAC",
        pbs_length=13,
        rt_length=30,
        target_strand="Fw",
        editing_position=50,
        pe2_hek293t=pe2_hek293t,
        pe2_k562=pe2_k562,
        pegrna_full="",
    )


class _OneRunner:
    """Injected runner returning a single pegRNA with a caller-chosen spacer."""

    def __init__(self, spacer: str) -> None:
        self._spacer = spacer

    def run_batch(self, rows):
        name = rows[0]["sequence_name"]
        return {name: {"pegrnas": [_peg(self._spacer)], "nicking": []}}


class _TwoRunner:
    """Returns two pegRNAs so the top-PE2 selection policy can be exercised."""

    def run_batch(self, rows):
        name = rows[0]["sequence_name"]
        # k562 (the default score column) prefers ``high``; hek293t would prefer
        # ``low`` -- so a correct selection proves the k562 policy is used.
        low = _peg("AAAACCCCGGGGTTTTAAAA", pe2_k562=0.10, pe2_hek293t=0.90)
        high = _peg("CAGTCAGTCAGTCAGTCAGT", pe2_k562=0.90, pe2_hek293t=0.10)
        return {name: {"pegrnas": [low, high], "nicking": []}}


class _EmptyRunner:
    """Returns no pegRNA for the pegRNA#1 request (fail-loud trigger)."""

    def run_batch(self, rows):
        name = rows[0]["sequence_name"]
        return {name: {"pegrnas": [], "nicking": []}}


class _BoomRunner:
    """Must never be called: proves the override path skips the runner entirely."""

    def run_batch(self, rows):  # pragma: no cover - asserts it is not invoked
        raise AssertionError("PRIDICT runner must not be called when a spacer is supplied")


# ============================================================== AUTO path


def test_auto_spacer_uses_bundle_runner_when_no_spacer_given():
    """No spacer + fake mode: the fake scorer bundle's injected PRIDICT runner
    designs the spacer, and the finished design carries it."""
    out = recoded_driver.design_recoded(**_auto_inputs(), mode="fake")
    assert out.result.pegRNA_sequence  # a full design was produced
    assert out.result.spacer == _FAKE_BUNDLE_SPACER


def test_auto_spacer_uses_explicit_injected_runner_spacer():
    """An explicit ``pridict_runner=`` overrides the bundle runner; the pipeline
    uses that runner's pegRNA#1 spacer."""
    out = recoded_driver.design_recoded(
        **_auto_inputs(), mode="fake", pridict_runner=_OneRunner("TGCATGCATGCATGCATGCA")
    )
    assert out.result.spacer == "TGCATGCATGCATGCATGCA"


def test_auto_spacer_selects_top_pe2_pegrna():
    """When PRIDICT returns several pegRNAs, the highest default-cell-line (k562)
    PE2 design is chosen."""
    out = recoded_driver.design_recoded(**_auto_inputs(), mode="fake", pridict_runner=_TwoRunner())
    assert out.result.spacer == "CAGTCAGTCAGTCAGTCAGT"


# ============================================================== OVERRIDE path


def test_explicit_spacer_uses_fixed_path_and_skips_runner():
    """A supplied spacer preserves the historical fixed-spacer behaviour and never
    invokes the PRIDICT runner (even when one is injected)."""
    out = recoded_driver.design_recoded(
        **_auto_inputs(),
        spacer="ACGTACGTACGTACGTACGT",
        mode="fake",
        pridict_runner=_BoomRunner(),
    )
    assert out.result.spacer == "ACGTACGTACGTACGTACGT"


# ============================================================== FAIL-LOUD


def test_auto_spacer_fails_loud_when_runner_returns_no_pegrna():
    with pytest.raises(RuntimeError, match="no pegRNA"):
        recoded_driver.design_recoded(**_auto_inputs(), mode="fake", pridict_runner=_EmptyRunner())


def test_auto_spacer_fails_loud_when_no_runner_available(monkeypatch):
    """If auto mode is requested but neither the bundle nor an explicit kwarg
    supplies a runner, the design fails loud rather than guessing a spacer."""

    def _bundle_without_runner(*, mode, **kwargs):  # noqa: ARG001
        bundle = scorer_factory._fake_bundle()
        bundle.pridict_runner = None
        return bundle

    monkeypatch.setattr(recoded_driver, "make_scorers", _bundle_without_runner)
    with pytest.raises(RuntimeError, match="PRIDICT runner"):
        recoded_driver.design_recoded(**_auto_inputs(), mode="fake")


def test_no_spacer_with_auto_disabled_fails_loud():
    with pytest.raises(ValueError, match="auto_spacer"):
        recoded_driver.design_recoded(**_auto_inputs(), mode="fake", auto_spacer=False)


# ============================================================== insert truncation


class _RecordingRunner:
    """Captures the editseq PRIDICT is handed, then returns one gate-neutral pegRNA."""

    def __init__(self, spacer: str = "ACGTACGTACGTACGTACGT") -> None:
        self._spacer = spacer
        self.editseq = None

    def run_batch(self, rows):
        self.editseq = rows[0]["editseq"]
        name = rows[0]["sequence_name"]
        return {name: {"pegrnas": [_peg(self._spacer)], "nicking": []}}


def _bracket_insert(editseq: str) -> str:
    """Extract the ``(+INSERT)`` insertion payload from a PRIDICT editseq."""
    return editseq.split("(+", 1)[1].split(")", 1)[0]


def test_pridict_editseq_truncates_cargo_to_40bp():
    """PRIDICT2 rejects insertions > 40 bp, and here it is only a spacer/nick designer,
    so the design editseq carries at most the cargo's leading 40 bp (the full cargo still
    flows into the tj-PE assembly downstream). Regression for the 503 bp CTCF cargo that
    silently produced no PRIDICT output on the DSI box."""
    rec = _RecordingRunner()
    cargo = "ATGC" * 20  # 80 bp, well over PRIDICT's 40 bp ceiling
    genome = "ACGATGCATGGCATGACTGA" * 16  # 320 bp -> flanks 150/170 at site_s=150
    spacer, _ = recoded_driver._design_peg1_via_pridict(
        genome=genome,
        site_s=150,
        cargo=cargo,
        auto_spacer=True,
        min_flank=100,
        runner=rec,
        transcript=None,
        context_fetcher=None,
    )
    assert spacer == "ACGTACGTACGTACGTACGT"
    ins = _bracket_insert(rec.editseq)
    assert len(ins) == 40
    assert ins == cargo[:40].upper()


def test_pridict_editseq_keeps_short_cargo_whole():
    """A cargo <= 40 bp is passed to PRIDICT unchanged."""
    rec = _RecordingRunner()
    cargo = "ATGAAACCCGGG"  # 12 bp
    genome = "ACGATGCATGGCATGACTGA" * 16
    recoded_driver._design_peg1_via_pridict(
        genome=genome,
        site_s=150,
        cargo=cargo,
        auto_spacer=True,
        min_flank=100,
        runner=rec,
        transcript=None,
        context_fetcher=None,
    )
    assert _bracket_insert(rec.editseq) == cargo.upper()
