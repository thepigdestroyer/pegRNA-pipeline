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
