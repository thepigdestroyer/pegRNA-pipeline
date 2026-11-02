"""Task 8: splice-disruption gate (Protocol + fake; real SpliceAI/Pangolin deferred to CX3)."""

import pytest

from pegrna_pipeline.hpc.splice_screen import (
    SpliceFake,
    passes_splice_gate,
)


def _spliceai_installed() -> bool:
    """True when the real SpliceAI stack is importable in this environment."""
    import importlib.util

    try:
        return importlib.util.find_spec("spliceai") is not None and (
            importlib.util.find_spec("tf_keras") is not None
            or importlib.util.find_spec("keras") is not None
        )
    except (ImportError, ValueError):  # pragma: no cover - defensive
        return False


# The two tests below assert the adapter FAILS LOUD when the weights are absent.
# That contract is only observable without the deps, so they are skipped (not
# failed) on a host where the real `[splice]` stack is installed.
needs_no_spliceai = pytest.mark.skipif(
    _spliceai_installed(),
    reason="asserts the missing-dependency path; the real SpliceAI stack is installed here",
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


@needs_no_spliceai
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


@needs_no_spliceai
def test_real_adapter_handles_indel_length_mismatch():
    # Indel edits (recoded-TJPE: the edited allele is longer/shorter than the
    # reference by the cargo insert) are now handled via prefix/suffix alignment,
    # NOT rejected. With adequate context but no weights, delta reaches the model
    # import and fails loud there — it must NOT raise the old "length-matched"
    # ValueError. (The indel-delta math itself is unit-tested with a stubbed
    # _tracks in tests/test_scorer_fixes.py.)
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    with pytest.raises((ImportError, RuntimeError)) as exc:
        SpliceAIScorer().delta("A" * 11000, "A" * 11000 + "GGGG")  # 4-nt insert
    assert "length-matched" not in str(exc.value)
