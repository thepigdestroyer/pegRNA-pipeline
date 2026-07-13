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

import numpy as np
import pytest

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


# ---------------------------------------------------------- SpliceAI indel delta
# The real SpliceAIScorer.delta indel branch is CX3-only (needs TF weights), but its
# coordinate logic is testable by injecting synthetic score tracks. These pin the fix
# for the broadcast bug that made every tjPE insertion advisory-skip: the score tracks
# cover the input minus the trimmed context, so prefix/suffix must align in TRACK
# coordinates, not raw-input coordinates.


def _spliceai_with_fake_tracks(monkeypatch, context, donor_spike_len=None, spike_pos=0):
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    sc = SpliceAIScorer(context=context)

    def _fake_tracks(seq):
        n = len(seq) - 2 * context  # the model trims `context` from each end
        acceptor = np.zeros(max(n, 0))
        donor = np.zeros(max(n, 0))
        if donor_spike_len is not None and n == donor_spike_len:
            donor[spike_pos] = 0.9
        return acceptor, donor

    monkeypatch.setattr(sc, "_tracks", _fake_tracks)
    return sc


def test_spliceai_delta_indel_aligns_in_track_coordinates(monkeypatch):
    # ref core "AC" (track len 2), alt core "GGGAC" (track len 5, a 3-nt insertion),
    # each N-padded by context=3. On raw-input coords the common prefix is the 3 N's,
    # which sliced the length-2 and length-5 tracks into a (2,) vs (3,) broadcast error.
    sc = _spliceai_with_fake_tracks(monkeypatch, context=3)
    ref = "NNN" + "AC" + "NNN"
    alt = "NNN" + "GGGAC" + "NNN"
    assert sc.delta(ref, alt) == 0.0  # aligns cleanly; zero tracks -> zero delta


def test_spliceai_delta_indel_reports_gain_in_inserted_region(monkeypatch):
    # A donor spike inside the inserted "GGG" region is reported as a splice GAIN.
    sc = _spliceai_with_fake_tracks(monkeypatch, context=3, donor_spike_len=5, spike_pos=1)
    ref = "NNN" + "AC" + "NNN"
    alt = "NNN" + "GGGAC" + "NNN"
    assert sc.delta(ref, alt) == pytest.approx(0.9)
