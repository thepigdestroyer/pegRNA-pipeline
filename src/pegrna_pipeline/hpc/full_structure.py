"""Full-pegRNA 2D structure check (decision 7a).

Primary scorer = **RiNALMo**: its SS head emits an L×L symmetric base-pair
PROBABILITY matrix (RoPE → no fixed length cap, handles ≤~700 nt). We read off
(a) cross-region pairing coupling DISTINCT functional regions (spacer, PBS2,
RTT, PBS1) and (b) per-base accessibility = 1 − Σ pair-probs. Scaffold rows/cols
are EXCLUDED (scaffold is designed to self-structure).

``FullStructureScorer`` is a Protocol so RiNALMo (GPU/Ampere+), RNAformer
(CPU/CI), or a fake are interchangeable. ViennaRNA (``rnafold.py``) stays the
cheap inner-loop pre-screen. The pure scoring logic here is fully testable with
the fake; the real adapters load heavy deps lazily and run on CX3.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Protocol, Tuple

__all__ = [
    "FullStructureScorer",
    "FullStructureFake",
    "RiNALMoScorer",
    "RNAformerScorer",
    "StructureScore",
    "functional_region_pairing",
    "region_accessibility",
]


class FullStructureScorer(Protocol):
    def pair_probs(self, seq: str) -> Dict[Tuple[int, int], float]:
        """Return {(i, j): probability} for i < j (sparse; absent pairs ~ 0)."""
        ...


class FullStructureFake:
    """Deterministic fake for CI/tests."""

    def __init__(self, pairs: Dict[Tuple[int, int], float]) -> None:
        self._pairs = pairs

    def pair_probs(self, seq: str) -> Dict[Tuple[int, int], float]:  # noqa: ARG002
        return self._pairs


@dataclass(frozen=True)
class StructureScore:
    cross_region_pairs: int


def _region_of(i: int, regions: Dict[str, Tuple[int, int]]):
    for name, (lo, hi) in regions.items():
        if lo <= i < hi:
            return name
    return None


def functional_region_pairing(
    seq: str, *, regions: Dict[str, Tuple[int, int]],
    scorer: FullStructureScorer, threshold: float = 0.5,
) -> StructureScore:
    """Count pairings (prob ≥ ``threshold``) that couple two DISTINCT functional
    regions, excluding any pair touching the scaffold."""
    pairs = scorer.pair_probs(seq)
    count = 0
    for (i, j), p in pairs.items():
        if p < threshold:
            continue
        ri, rj = _region_of(i, regions), _region_of(j, regions)
        if ri is None or rj is None:
            continue
        if ri == "scaffold" or rj == "scaffold":
            continue
        if ri != rj:
            count += 1
    return StructureScore(cross_region_pairs=count)


def region_accessibility(
    seq: str, *, region: Tuple[int, int], scorer: FullStructureScorer
) -> float:
    """Mean per-base unpaired probability over ``region`` = 1 − Σ pair-probs."""
    pairs = scorer.pair_probs(seq)
    lo, hi = region
    width = max(1, hi - lo)
    paired = {k: 0.0 for k in range(lo, hi)}
    for (i, j), p in pairs.items():
        if lo <= i < hi:
            paired[i] += p
        if lo <= j < hi:
            paired[j] += p
    return sum(max(0.0, 1.0 - paired[k]) for k in range(lo, hi)) / width


class RiNALMoScorer:
    """Real RiNALMo adapter (CX3-only; Ampere+ GPU + FlashAttention + weights).

    RiNALMo's fine-tuned secondary-structure head outputs an L×L symmetric
    base-pair-probability matrix; we sparsify it to ``{(i, j): p}`` for i < j.
    Heavy deps (torch + the 650M-param weights) load lazily so the core package
    + CI never import them. Raises a helpful error off-GPU — never silently
    skipped.
    """

    def __init__(self, *, weights: str | None = None, threshold: float = 0.5) -> None:
        self.weights = weights
        self.threshold = threshold
        self._model = None

    def _load(self):  # pragma: no cover - exercised only on Ampere+ GPU
        if self._model is not None:
            return
        try:
            import torch  # noqa: F401, PLC0415
            from rinalmo.pretrained import get_pretrained_model  # noqa: F401, PLC0415
        except Exception as exc:
            raise ImportError(
                "RiNALMoScorer requires the '[fullstructure]' extra (torch + rinalmo) "
                "and an Ampere/Ada/Hopper GPU with FlashAttention; run on CX3. "
                f"Original error: {exc}"
            ) from exc
        raise NotImplementedError(
            "RiNALMo SS-head inference + bpp sparsification is wired in the follow-on plan"
        )

    def pair_probs(self, seq: str) -> Dict[Tuple[int, int], float]:  # pragma: no cover
        self._load()
        raise NotImplementedError


class RNAformerScorer:
    """Real RNAformer adapter (conda-forge, CPU-capable → CI-friendly fallback).

    RNAformer predicts a base-pair probability map from sequence alone (no MSA).
    Loaded lazily; used as the cross-check / CI fallback when a GPU + RiNALMo are
    unavailable.
    """

    def __init__(self, *, model_dir: str | None = None) -> None:
        self.model_dir = model_dir
        self._model = None

    def _load(self):  # pragma: no cover - exercised only with rnaformer installed
        if self._model is not None:
            return
        try:
            import RNAformer  # noqa: F401, PLC0415
        except Exception as exc:
            raise ImportError(
                "RNAformerScorer requires the '[fullstructure]' extra (rnaformer). "
                f"Original error: {exc}"
            ) from exc
        raise NotImplementedError(
            "RNAformer inference + bpp extraction is wired in the follow-on plan"
        )

    def pair_probs(self, seq: str) -> Dict[Tuple[int, int], float]:  # pragma: no cover
        self._load()
        raise NotImplementedError
