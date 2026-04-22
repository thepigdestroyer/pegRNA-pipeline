"""Tethered-junction prime editing (TJ-PE) workflow.

Entry point: :func:`run_tjpe_workflow`. See Plan 2 v1.1 for the pipeline
contract. Tasks here:

* :func:`rank_variants_by_pe2` — T6
* :func:`select_tjpe_ngrna`   — T7
* :func:`derive_rc_pbs2`      — T8
* :func:`assemble_tj_pegrna`  — T9
* :class:`TJPEResult`         — T10
* :func:`run_tjpe_workflow`   — T10
"""

from __future__ import annotations

from typing import List, Optional

from pegrna_pipeline.hpc.data_types import (
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)


OPPOSITE_STRAND = {"Fw": "Rv", "Rv": "Fw"}


# ------------------------------------------------------------------ T6


def rank_variants_by_pe2(variants: List[VariantScore]) -> List[VariantScore]:
    """Return a new list, sorted by ``top_score`` descending, stably.

    Stable sort ensures callers can rely on insertion order for ties.
    """
    return sorted(variants, key=lambda v: v.top_score, reverse=True)


# ------------------------------------------------------------------ T7


def _spacer_passes_sanity(spacer: str) -> bool:
    """Reject unusable spacers on their own sequence content.

    * ``TTTT`` is a RNA-pol III terminator — guide will truncate.
    * A homopolymer run of >= 5 bp kills expression.
    * GC content outside [0.30, 0.75] gives poor DNA binding.
    """
    s = spacer.upper()
    if "TTTT" in s:
        return False
    # Homopolymer >= 5
    run = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            run += 1
            if run >= 5:
                return False
        else:
            run = 1
    if not s:
        return False
    gc = (s.count("G") + s.count("C")) / len(s)
    if gc < 0.30 or gc > 0.75:
        return False
    return True


def select_tjpe_ngrna(
    candidates: List[NgRNADesign],
    *,
    min_distance: int = 40,
    max_distance: int = 100,
    pegrna_strand: Optional[str] = None,
) -> Optional[NgRNADesign]:
    """Pick the best (highest ``deepcas9_score``) surviving ngRNA.

    Parameters
    ----------
    candidates:
        List of :class:`NgRNADesign` rows, typically parsed from a PRIDICT
        nicking-guides CSV.
    min_distance, max_distance:
        Inclusive bounds on ``|distance_to_edit|``.
    pegrna_strand:
        If given (``"Fw"`` or ``"Rv"``), only keep opposite-strand
        candidates: ``{"Fw": "Rv", "Rv": "Fw"}``. Raises :class:`KeyError`
        for any other value — unknown strands are always a bug, not a
        silent filter.

    Returns
    -------
    The candidate with the maximum ``deepcas9_score`` that passes every
    filter, or ``None`` if no candidate survives.
    """
    required_strand: Optional[str] = None
    if pegrna_strand is not None:
        required_strand = OPPOSITE_STRAND[pegrna_strand]  # raises KeyError if unknown

    survivors: List[NgRNADesign] = []
    for c in candidates:
        d = abs(c.distance_to_edit)
        if d < min_distance or d > max_distance:
            continue
        if required_strand is not None and c.target_strand != required_strand:
            continue
        if not _spacer_passes_sanity(c.spacer):
            continue
        survivors.append(c)

    if not survivors:
        return None
    return max(survivors, key=lambda c: c.deepcas9_score)
