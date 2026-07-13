"""api-mode splice gate wires the N-padding SpliceAI adapter (audit fix).

The GPU-less web host (``mode="api"``) usually only has the short local design
window (hundreds of nt), but the real ``SpliceAIScorer`` fail-louds unless handed
>= 10 kb of context (>= 5 kb each side). ``_PaddedSpliceAIScorer`` N-pads short
windows up to that requirement so the SpliceAI advisory cross-check actually runs
on the web host. It was defined but never wired into ``_api_bundle``, so on the
deployment the cross-check was silently skipped (``run_recoded_tjpe`` catches the
short-context ``ValueError`` and records "advisory skip"). These tests pin the
padding behaviour and the api-bundle wiring.

Real mode is unchanged and intentionally so: it supplies the true GRCh38 flanks,
so ``_real_bundle`` keeps the raw ``SpliceAIScorer`` (asserted here as a
regression guard, since broadening the shared workflow's exception handling would
break real mode's fail-loud contract).
"""

from __future__ import annotations

import pegrna_pipeline.hpc.splice_screen as splice_screen
from pegrna_pipeline.hpc import scorer_factory
from pegrna_pipeline.hpc.scorer_factory import _PaddedSpliceAIScorer


class _RecordingInner:
    """Stand-in for ``SpliceAIScorer`` that records the (ref, alt) it is handed
    and never touches TensorFlow, so the padding math is testable in CI."""

    _context = 5000

    def __init__(self) -> None:
        self.calls: list[tuple[str, str]] = []

    def delta(self, ref_seq: str, alt_seq: str) -> float:
        self.calls.append((ref_seq, alt_seq))
        return 0.0


def _padded_with_recording_inner(monkeypatch):
    # _PaddedSpliceAIScorer.__init__ does `from ...splice_screen import SpliceAIScorer`,
    # so patching the attribute on that module injects the recording inner.
    monkeypatch.setattr(splice_screen, "SpliceAIScorer", _RecordingInner)
    scorer = _PaddedSpliceAIScorer()
    return scorer, scorer._inner


def test_padded_scorer_pads_short_window_to_full_context(monkeypatch):
    scorer, inner = _padded_with_recording_inner(monkeypatch)
    scorer.delta("ACGTACGT", "ACGTACGTGGG")  # far below 10 kb; an insertion edit
    assert len(inner.calls) == 1
    ref, alt = inner.calls[0]
    # Padded to >= 2*context on both sides so SpliceAIScorer will not fail-loud.
    assert len(ref) >= 2 * inner._context
    assert len(alt) >= 2 * inner._context
    # N is SpliceAI's zero-signal pad; the edit region is preserved intact between
    # the flanks, so the delta over the real bases is unchanged.
    flank = "N" * inner._context
    assert ref == flank + "ACGTACGT" + flank
    assert alt == flank + "ACGTACGTGGG" + flank


def test_padded_scorer_passes_through_adequate_context(monkeypatch):
    scorer, inner = _padded_with_recording_inner(monkeypatch)
    ref = "A" * (2 * inner._context)
    alt = "A" * (2 * inner._context)
    scorer.delta(ref, alt)
    # Already long enough -> no padding, handed through verbatim.
    assert inner.calls[0] == (ref, alt)


def test_api_bundle_wires_padded_spliceai(monkeypatch):
    # Stub the heavy api scorers (torch/TF/ViennaRNA/biotite absent in CI) so the
    # bundle can be built; leave the light splice scorer real to assert its type.
    # No NVIDIA key -> DockFake path, so complex_dock_nim is never touched.
    from pegrna_pipeline.hpc import cryptic_splice, structure_importance, vienna_structure

    monkeypatch.setattr(vienna_structure, "ViennaRNAStructureScorer", lambda **_: object())
    monkeypatch.setattr(structure_importance, "AFDBImportanceScorer", lambda *_, **__: object())
    monkeypatch.setattr(cryptic_splice, "MaxEntScanScorer", lambda **_: object())

    bundle = scorer_factory.make_scorers(mode="api")
    assert isinstance(bundle.splice_scorer, _PaddedSpliceAIScorer)


def test_real_bundle_keeps_raw_spliceai():
    # Regression guard: real mode supplies true GRCh38 flanks, so it must keep the
    # raw SpliceAIScorer (NOT the padded web-host adapter). make_scorers(mode="real")
    # constructs lazily in CI (see test_real_mode_raises_without_heavy_deps).
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    bundle = scorer_factory.make_scorers(mode="real")
    assert type(bundle.splice_scorer) is SpliceAIScorer
    assert not isinstance(bundle.splice_scorer, _PaddedSpliceAIScorer)
