"""An injected ngRNA the walk will immediately exclude must not count as a pool.

The auto-spacer path screens PRIDICT's nicking guides and injects the survivors
as the ENTIRE candidate pool; only an EMPTY result falls back to the genome-
enumerated nicks. But the arm-side rule (arm_side=5prime needs the nick 5' of the
insertion site, 3prime needs it 3') was applied later, inside the walk. So when
the only PRIDICT guide that survived off-target screening sat on the wrong side
of S, the injected list was non-empty, no fallback fired, and the walk excluded
that one guide: "no candidate survived all gates after 0 tried" -- while a clean,
usable genome-enumerated nick (nick 101 in the CTCF example) sat unused.

This is exactly what the shipped real-models demo hit once off-target screening
became real: real cas-offinder removed PRIDICT's 5' guides, leaving only a 3' one.
"""

from pegrna_pipeline.hpc.pridict_ngrna import build_pridict_ngrnas

SITE_S = 147
# Distinct 20-mers that pass the spacer sanity filter (no TTTT, no run of >=5, GC 30-75%);
# their sequence is otherwise irrelevant to the arm-side logic.
UPSTREAM = "GCTGCCACAGATGCCGCAAC"
DOWNSTREAM = "TTCTCTGGGCTGTGTTTAAA"
GENOME = "ACGT" * 75  # 300 bp: room for a PBS2 window on either side of S


def _row(spacer: str, distance: int, score: float = 0.6) -> dict:
    return {
        "Nicking-Protospacer": spacer,
        "Nicking-Position-to-edit": str(distance),
        "Nicking-PAMdisrupt": "Nicking_PAM_disrupt",
        "Target_Strand": "Fw",
        "DeepCas9score": str(score),
    }


def _clean(protospacers):
    return {p: (1.0, 0, False) for p in protospacers}


def _build(rows, arm_side):
    return build_pridict_ngrnas(
        rows,
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",  # ngRNA must be on the opposite strand (Fw)
        ngrna_distance_min=40,
        ngrna_distance_max=90,
        offtarget_scorer=_clean,
        arm_side=arm_side,
    )


def test_a_guide_on_the_wrong_side_is_not_injected_for_5prime():
    out = _build([_row(DOWNSTREAM, +56)], arm_side="5prime")
    assert (
        out == []
    ), "3' of S cannot serve arm_side=5prime; an empty list lets the caller fall back"


def test_a_guide_on_the_wrong_side_is_not_injected_for_3prime():
    out = _build([_row(UPSTREAM, -46)], arm_side="3prime")
    assert out == []


def test_the_usable_side_is_kept_and_the_unusable_one_dropped():
    rows = [_row(DOWNSTREAM, +56, score=0.9), _row(UPSTREAM, -46, score=0.4)]
    out = _build(rows, arm_side="5prime")
    assert [n.protospacer for n in out] == [UPSTREAM]
    assert out[0].nick_position == SITE_S - 46


def test_auto_arm_side_keeps_both():
    rows = [_row(DOWNSTREAM, +56, score=0.9), _row(UPSTREAM, -46, score=0.4)]
    out = _build(rows, arm_side="auto")
    assert {n.protospacer for n in out} == {UPSTREAM, DOWNSTREAM}


def test_unusable_guides_are_not_sent_to_the_expensive_offtarget_scan():
    """cas-offinder costs ~100 s per invocation; do not spend it on guides that
    the walk is certain to throw away."""
    seen = []

    def spy(protospacers):
        seen.append(list(protospacers))
        return _clean(protospacers)

    build_pridict_ngrnas(
        [_row(DOWNSTREAM, +56), _row(UPSTREAM, -46)],
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=90,
        offtarget_scorer=spy,
        arm_side="5prime",
    )
    assert seen == [[UPSTREAM]]


def test_default_is_unchanged_for_existing_callers():
    """Omitting arm_side must behave as before (no side filtering)."""
    out = build_pridict_ngrnas(
        [_row(DOWNSTREAM, +56), _row(UPSTREAM, -46)],
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=90,
        offtarget_scorer=_clean,
    )
    assert len(out) == 2
