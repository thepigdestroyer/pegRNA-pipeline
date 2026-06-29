"""Task 8: splice-disruption gate (Protocol + fake; real SpliceAI/Pangolin deferred to CX3)."""

import pytest

from pegrna_pipeline.hpc.splice_screen import (
    SpliceFake,
    passes_splice_gate,
)


def test_splice_gate_blocks_high_delta():
    scorer = SpliceFake(delta=0.6)  # strong predicted disruption
    ok, detail = passes_splice_gate("ACGT" * 30, "ACGA" * 30, scorer=scorer, max_delta=0.2)
    assert ok is False
    assert "delta" in detail


def test_splice_gate_passes_low_delta():
    scorer = SpliceFake(delta=0.05)
    ok, _ = passes_splice_gate("ACGT" * 30, "ACGT" * 30, scorer=scorer, max_delta=0.2)
    assert ok is True


def test_splice_gate_boundary_is_blocking():
    scorer = SpliceFake(delta=0.2)  # exactly at threshold => block (>=)
    ok, _ = passes_splice_gate("A" * 50, "A" * 50, scorer=scorer, max_delta=0.2)
    assert ok is False


def test_real_adapter_importerror_without_weights():
    # With adequate, length-matched context the real adapter reaches the weight
    # import and fails loud (ImportError) absent deps — never silently skipped.
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    ref = "A" * 11000
    alt = "A" * 5500 + "T" + "A" * 5499
    with pytest.raises((ImportError, RuntimeError)):
        SpliceAIScorer().delta(ref, alt)


def test_real_adapter_fails_loud_on_short_context():
    # FAIL-LOUD, not fail-open: too-short context must raise, never return Δ=0.
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    with pytest.raises(ValueError, match="context"):
        SpliceAIScorer().delta("ACGT" * 30, "ACGA" * 30)


def test_real_adapter_rejects_length_mismatch():
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    with pytest.raises(ValueError, match="length-matched"):
        SpliceAIScorer().delta("A" * 11000, "A" * 10999)
