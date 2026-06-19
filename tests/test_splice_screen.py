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
    # The real adapter must raise a helpful ImportError (not crash) when the heavy
    # deps / weights are absent — it is never silently skipped.
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    with pytest.raises((ImportError, RuntimeError)):
        SpliceAIScorer().delta("ACGT" * 30, "ACGA" * 30)
