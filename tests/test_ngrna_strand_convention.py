"""Document the {Fw -> Rv, Rv -> Fw} opposite-strand selection assumption.

For PE3/PE3b the nicking guide must lie on the **opposite** strand from the
pegRNA. :func:`select_tjpe_ngrna` encodes this with a tiny mapping dict —
{"Fw": "Rv", "Rv": "Fw"}. If anyone ever changes PRIDICT's strand labels,
this test fails first and the cause is obvious.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import NgRNADesign
from pegrna_pipeline.hpc.tjpe_workflow import (
    OPPOSITE_STRAND,
    select_tjpe_ngrna,
)


def test_opposite_strand_mapping_is_bidirectional():
    assert OPPOSITE_STRAND == {"Fw": "Rv", "Rv": "Fw"}


def _ng(strand: str, distance: int = 50, score: float = 0.8) -> NgRNADesign:
    return NgRNADesign(
        # Realistic spacer: 50% GC, no TTTT, no homopolymer >= 5.
        spacer="AGTGCTACCTGGTACGGATC",
        distance_to_edit=distance,
        target_strand=strand,
        deepcas9_score=score,
        pe3b=False,
        pam_disrupt=False,
    )


def test_select_respects_opposite_strand_for_fw_pegrna():
    """A Fw pegRNA should pair with a Rv ngRNA and vice versa."""
    candidates = [_ng("Fw", 50, 0.9), _ng("Rv", 50, 0.5)]
    picked = select_tjpe_ngrna(candidates, pegrna_strand="Fw")
    assert picked is not None
    assert picked.target_strand == "Rv"


def test_select_respects_opposite_strand_for_rv_pegrna():
    candidates = [_ng("Fw", 50, 0.3), _ng("Rv", 50, 0.9)]
    picked = select_tjpe_ngrna(candidates, pegrna_strand="Rv")
    assert picked is not None
    assert picked.target_strand == "Fw"


def test_select_rejects_unknown_strand():
    """Any strand label outside {Fw,Rv} is caught early so silent bugs can't
    leak through the mapping."""
    with pytest.raises(KeyError):
        select_tjpe_ngrna([_ng("Fw", 50, 0.9)], pegrna_strand="XX")
