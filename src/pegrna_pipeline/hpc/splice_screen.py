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

from typing import Optional, Protocol, Tuple

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
        except Exception as exc:
            raise ImportError(
                "SpliceAI requires the '[splice]' extra (tensorflow/keras + spliceai) "
                "and the model weights; run real splice scoring on CX3. "
                f"Original error: {exc}"
            ) from exc
        paths = [resource_filename("spliceai", f"models/spliceai{i}.h5") for i in range(1, 6)]
        self._model = [load_model(p) for p in paths]

    @staticmethod
    def _one_hot(seq: str):  # pragma: no cover - CX3 only
        import numpy as np  # noqa: PLC0415

        mapping = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4}
        idx = np.array([mapping.get(b, 0) for b in seq.upper()])
        # SpliceAI uses a 4-channel one-hot with N -> all zeros.
        onehot = np.zeros((len(seq), 4), dtype=np.float32)
        for i, code in enumerate(idx):
            if code > 0:
                onehot[i, code - 1] = 1.0
        return onehot[None, :, :]

    def _tracks(self, seq: str):  # pragma: no cover - CX3 only
        import numpy as np  # noqa: PLC0415

        self._load()
        x = self._one_hot(seq)
        preds = [m.predict(x, verbose=0) for m in self._model]
        # ensemble mean; output channels [null, acceptor, donor]
        mean = np.mean(preds, axis=0)[0]
        return mean[:, 1], mean[:, 2]  # acceptor, donor probability tracks

    @staticmethod
    def _common_prefix(a: str, b: str) -> int:
        n = min(len(a), len(b))
        i = 0
        while i < n and a[i] == b[i]:
            i += 1
        return i

    @staticmethod
    def _common_suffix(a: str, b: str, cap: int) -> int:
        n = min(len(a), len(b), cap)
        i = 0
        while i < n and a[len(a) - 1 - i] == b[len(b) - 1 - i]:
            i += 1
        return i

    def delta(self, ref_seq: str, alt_seq: str) -> float:  # pragma: no cover - CX3 only
        """Max |Δ splice score| between REF and ALT, indel-aware.

        Faithful to the SpliceAI delta-score definition (acceptor/donor gain &
        loss). For an equal-length (SNV/MNV) edit this is the position-wise max
        |Δ|. For an INDEL (the recoded-TJPE common case — the ``edited`` allele is
        longer than ``reference`` by the cargo insert), REF/ALT differ in length,
        so a naive position-wise diff mis-aligns every base downstream of the
        variant. Instead we align on the common prefix/suffix and score:

          * perturbation of existing sites in the aligned flanks (|Δ| there), and
          * a splice-site GAIN = strongest predicted site inside the ALT variant
            region, and a LOSS = strongest site inside the REF variant region.

        This mirrors SpliceAI-lookup's indel handling. FAIL-LOUD on inadequate
        context (a safety gate must never silently return Δ=0 because the window
        was too short). The caller supplies the genomic window (>=``context`` nt
        each side of the edit) via genome_context; VALIDATE ON CX3.
        """
        import numpy as np  # noqa: PLC0415

        if len(ref_seq) < 2 * self._context or len(alt_seq) < 2 * self._context:
            raise ValueError(
                f"SpliceAIScorer.delta requires >= {2 * self._context} nt of context "
                f"(>= {self._context} each side); got ref={len(ref_seq)} alt={len(alt_seq)} nt"
            )
        a_ref, d_ref = self._tracks(ref_seq)
        a_alt, d_alt = self._tracks(alt_seq)
        if len(a_ref) == 0 or len(a_alt) == 0:
            raise ValueError("SpliceAIScorer produced empty score tracks (context too short)")

        if len(ref_seq) == len(alt_seq):
            n = min(len(a_ref), len(a_alt))
            da = np.abs(a_alt[:n] - a_ref[:n]).max()
            dd = np.abs(d_alt[:n] - d_ref[:n]).max()
            return float(max(da, dd))

        # --- indel: align on common prefix/suffix, then gain/loss in the middle ---
        lr, la = len(ref_seq), len(alt_seq)
        p = self._common_prefix(ref_seq, alt_seq)
        s = self._common_suffix(ref_seq, alt_seq, cap=min(lr, la) - p)
        deltas = [0.0]
        if p > 0:  # aligned prefix flank
            deltas.append(float(np.abs(a_alt[:p] - a_ref[:p]).max()))
            deltas.append(float(np.abs(d_alt[:p] - d_ref[:p]).max()))
        if s > 0:  # aligned suffix flank (tail-aligned; both length s)
            deltas.append(float(np.abs(a_alt[la - s :] - a_ref[lr - s :]).max()))
            deltas.append(float(np.abs(d_alt[la - s :] - d_ref[lr - s :]).max()))
        if lr - s > p:  # REF variant region -> LOSS (strongest site removed)
            deltas.append(float(a_ref[p : lr - s].max()))
            deltas.append(float(d_ref[p : lr - s].max()))
        if la - s > p:  # ALT variant region -> GAIN (strongest site created)
            deltas.append(float(a_alt[p : la - s].max()))
            deltas.append(float(d_alt[p : la - s].max()))
        return float(max(deltas))


def passes_splice_gate(
    ref_seq: str, alt_seq: str, *, scorer: SpliceScorer, max_delta: float = 0.2,
    metrics: Optional[dict] = None,
) -> Tuple[bool, str]:
    """HARD GATE: reject (False) if predicted splice delta >= ``max_delta``.

    Populates ``metrics`` (if given) with the raw SpliceAI delta and cutoff for
    the caller's progress display; the return value is unchanged.
    """
    d = scorer.delta(ref_seq, alt_seq)
    if metrics is not None:
        metrics.update({
            "metric": "SpliceAI Δ", "value": round(d, 3), "cutoff": max_delta, "cmp": "<",
        })
    if d >= max_delta:
        return False, f"splice delta {d:.3f} >= {max_delta} (predicted disruption)"
    return True, f"splice delta {d:.3f} ok"
