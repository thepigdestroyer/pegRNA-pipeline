"""Score upstream/downstream pegRNA pairs for Prime Assembly.

Compatibility rules (Plan 2 v1.1):

* Target strands are opposite (one Fw, one Rv).
* ``|edit_pos_up - edit_pos_dn|`` within [``min_distance``, ``max_distance``].
* Distance is also >= :data:`MIN_NICK_SEPARATION_BP` (hard floor).

The joint quality score is the geometric mean of the two PE2 efficiency
predictions, drawn from PRIDICT 2.0's K562 column by default (override
via ``score_cell_line``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List

from pegrna_pipeline.hpc.data_types import (
    DEFAULT_SCORE_CELL_LINE,
    PegRNADesign,
    _score_attr,
)

#: Hard floor on the distance between the two pegRNA nicks. Closer than
#: this and the two edits interfere, so no pair ever survives — even if
#: the caller relaxes ``min_distance`` below this value.
MIN_NICK_SEPARATION_BP = 23


@dataclass(frozen=True)
class PegRNAPair:
    """An upstream/downstream pegRNA combination with a joint score."""

    upstream: PegRNADesign
    downstream: PegRNADesign
    nick_distance: int
    joint_score: float


def filter_compatible_pairs(
    *,
    upstream: List[PegRNADesign],
    downstream: List[PegRNADesign],
    stop_codon_position: int,  # noqa: ARG001 — kept for API clarity / future use
    min_distance: int = 30,
    max_distance: int = 80,
    score_cell_line: str = DEFAULT_SCORE_CELL_LINE,
) -> List[PegRNAPair]:
    """Return all valid upstream/downstream pairs, sorted by joint score desc.

    Parameters
    ----------
    upstream, downstream:
        Candidate pegRNAs for each flap. In practice these come from two
        independent PRIDICT runs (one per flap).
    stop_codon_position:
        Reserved for future use: downstream nick-arithmetic currently
        lives in :func:`run_pa_workflow`, but callers may still want to
        pass the stop codon position for traceability.
    min_distance, max_distance:
        Inclusive bounds on ``|editing_position_up - editing_position_dn|``.

    Returns
    -------
    A list of :class:`PegRNAPair`, sorted by ``joint_score`` descending.
    An empty list is returned when no Cartesian pairing survives the
    compatibility filters.
    """
    effective_min = max(min_distance, MIN_NICK_SEPARATION_BP)
    score_attr = _score_attr(score_cell_line)
    pairs: List[PegRNAPair] = []
    for up in upstream:
        for dn in downstream:
            if up.target_strand == dn.target_strand:
                continue
            distance = abs(up.editing_position - dn.editing_position)
            if distance < effective_min or distance > max_distance:
                continue
            up_score = float(getattr(up, score_attr))
            dn_score = float(getattr(dn, score_attr))
            joint = math.sqrt(max(up_score, 0.0) * max(dn_score, 0.0))
            pairs.append(
                PegRNAPair(
                    upstream=up,
                    downstream=dn,
                    nick_distance=distance,
                    joint_score=joint,
                )
            )
    pairs.sort(key=lambda p: p.joint_score, reverse=True)
    return pairs
