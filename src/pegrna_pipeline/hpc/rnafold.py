"""ViennaRNA wrapper for pegRNA secondary-structure screening.

Used by the scarless tjPE workflow to (a) pre-screen individual components and
(b) check the fully assembled pegRNA before / alongside PRIDICT ranking.

ViennaRNA is an *optional* dependency: it is imported lazily inside the
functions so importing this module never fails when ViennaRNA is absent.
``is_available()`` lets callers degrade gracefully (skip the structure gate
with a warning) rather than crash.

Thresholds here are deliberately permissive defaults — there is no published
pegRNA ΔG cut-off, so they MUST be calibrated on real efficiency data
(see spec §I / 50-Open-Questions).
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

__all__ = ["is_available", "fold_metrics", "passes_structure"]


def is_available() -> bool:
    """Return True if the ViennaRNA python bindings can be imported."""
    try:
        import RNA  # noqa: F401
    except Exception:
        return False
    return True


def fold_metrics(seq: str, *, regions: Optional[Dict[str, Tuple[int, int]]] = None) -> dict:
    """Return MFE, length-normalised MFE, and mean unpaired-probability per region.

    ``regions`` maps a name → (lo, hi) half-open 0-indexed slice of ``seq`` over
    which to average the per-base unpaired probability (1 = fully accessible).
    Raises ``RuntimeError`` if ViennaRNA is unavailable.
    """
    try:
        import RNA
    except Exception as exc:  # pragma: no cover - exercised only without ViennaRNA
        raise RuntimeError("ViennaRNA not available; install '.[scarless]'") from exc

    rna = seq.upper().replace("T", "U")
    n = len(rna)
    fc = RNA.fold_compound(rna)
    _structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = fc.bpp()  # 1-indexed base-pair probability matrix

    def unpaired_prob(i: int) -> float:
        """Probability that base ``i`` (0-indexed) is unpaired."""
        idx = i + 1
        paired = 0.0
        for j in range(1, n + 1):
            if j == idx:
                continue
            a, b = (idx, j) if idx < j else (j, idx)
            paired += bpp[a][b]
        return max(0.0, 1.0 - paired)

    acc: Dict[str, float] = {}
    for name, (lo, hi) in (regions or {}).items():
        lo = max(0, lo)
        hi = min(n, hi)
        width = max(1, hi - lo)
        acc[name] = sum(unpaired_prob(i) for i in range(lo, hi)) / width

    return {"mfe": float(mfe), "norm_mfe": float(mfe) / max(1, n), "accessibility": acc}


def passes_structure(
    seq: str,
    *,
    regions: Optional[Dict[str, Tuple[int, int]]] = None,
    max_norm_mfe: float = -0.5,
    min_accessibility: float = 0.5,
) -> bool:
    """Permissive structure gate (calibrate thresholds on real data, spec §I).

    Rejects (returns False) if the length-normalised MFE is more negative than
    ``max_norm_mfe`` (too structured) or any region's mean accessibility is below
    ``min_accessibility``.
    """
    m = fold_metrics(seq, regions=regions or {})
    if m["norm_mfe"] < max_norm_mfe:
        return False
    for value in m["accessibility"].values():
        if value < min_accessibility:
            return False
    return True
