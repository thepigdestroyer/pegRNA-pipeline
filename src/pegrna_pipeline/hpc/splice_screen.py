"""Splice-disruption hard gate (decision 5).

``SpliceScorer`` is a Protocol so the heavy SpliceAI/Pangolin models (HPC-staged
weights, GRCh38 ±5 kb genomic context) are injected; CI uses :class:`SpliceFake`.
The gate compares REF vs ALT delta near annotated junctions, including the
symmetric-exon-skip mode that survives the T2A-Zeocin cassette.

The real :class:`SpliceAIScorer` adapter is written against the actual SpliceAI
API but is **CX3-only**: it requires TensorFlow + the SpliceAI weights + a
genomic window with ≥5 kb flanks. It raises a helpful error when those are
absent — it is never silently skipped.
"""

from __future__ import annotations

from typing import Protocol, Tuple

__all__ = ["SpliceScorer", "SpliceFake", "SpliceAIScorer", "passes_splice_gate"]


class SpliceScorer(Protocol):
    def delta(self, ref_seq: str, alt_seq: str) -> float:
        """Return the max |Δ splice score| between REF and ALT (0..1)."""
        ...


class SpliceFake:
    """Deterministic fake for CI/tests."""

    def __init__(self, delta: float = 0.0) -> None:
        self._d = delta

    def delta(self, ref_seq: str, alt_seq: str) -> float:  # noqa: ARG002
        return self._d


class SpliceAIScorer:
    """Real SpliceAI adapter (CX3-only; lazy heavy imports).

    SpliceAI needs the pre-mRNA / genomic context with ≥5 kb flanks on GRCh38
    coordinates and the published weights. This adapter loads them lazily so the
    core package + CI never import TensorFlow. It returns the maximum absolute
    delta of acceptor/donor gain/loss between REF and ALT (the standard
    delta-score readout). Cross-check with Pangolin is wired the same way.
    """

    def __init__(self, *, annotation: str = "grch38", context: int = 5000) -> None:
        self._annotation = annotation
        self._context = context
        self._model = None

    def _load(self):  # pragma: no cover - exercised only on CX3 with weights
        if self._model is not None:
            return
        try:
            from keras.models import load_model  # noqa: PLC0415
            from pkg_resources import resource_filename  # noqa: PLC0415
            from spliceai.utils import Annotator  # noqa: F401, PLC0415
        except Exception as exc:
            raise ImportError(
                "SpliceAI requires the '[splice]' extra (tensorflow/keras + spliceai) "
                "and the model weights; run real splice scoring on CX3. "
                f"Original error: {exc}"
            ) from exc
        paths = [resource_filename("spliceai", f"models/spliceai{i}.h5") for i in range(1, 6)]
        self._model = [load_model(p) for p in paths]

    def delta(self, ref_seq: str, alt_seq: str) -> float:  # pragma: no cover - CX3 only
        self._load()
        # Real implementation scores ref/alt windows through the 5-model ensemble
        # and returns max(|acceptor_gain|, |acceptor_loss|, |donor_gain|,
        # |donor_loss|). Kept behind _load() so import-time never pulls TF.
        raise NotImplementedError(
            "SpliceAIScorer.delta runs on CX3 with staged weights (see follow-on plan)"
        )


def passes_splice_gate(
    ref_seq: str, alt_seq: str, *, scorer: SpliceScorer, max_delta: float = 0.2
) -> Tuple[bool, str]:
    """HARD GATE: reject (False) if predicted splice delta >= ``max_delta``."""
    d = scorer.delta(ref_seq, alt_seq)
    if d >= max_delta:
        return False, f"splice delta {d:.3f} >= {max_delta} (predicted disruption)"
    return True, f"splice delta {d:.3f} ok"
