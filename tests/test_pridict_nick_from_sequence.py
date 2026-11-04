"""A PRIDICT nicking guide is placed by its SEQUENCE, never by PRIDICT's distance (issue #31).

``build_pridict_ngrnas`` used to place the second nick at ``site_s + Nicking-Position-to-edit``. Real
PRIDICT2 output for the CTCF window with a 24 bp insert (scripts/analysis/probe_pridict_nick_convention.py,
39 guides, all consistent) shows that number is not the position of the cut:

    forward-strand guide, upstream of S     distance = true - 1
    forward-strand guide, downstream of S   distance = true + (insert length - 1)
    reverse-strand guide, downstream of S   distance = -true
    reverse-strand guide, upstream of S     distance = -true + insert length

where ``true`` is the cut minus S, in the sense frame, with the cut the first base right of the cut. So the
injected nick was one base off upstream, an insert length off downstream, and mirrored to the wrong side of S
for a reverse-strand guide: an arm one base too long (a duplicated base) in the best case. Nothing downstream
could notice, because the whole design is built around the wrong number. The designer now finds the guide in
the window and computes the cut itself, so PRIDICT's number is not used at all.
"""

import json
from pathlib import Path

import pytest
from pa_helpers import guide_row, make_locus, nicking_row, pridict_distance, rc

from pegrna_pipeline.hpc.nick_search import locate_guide
from pegrna_pipeline.hpc.pa_geometry import find_guides
from pegrna_pipeline.hpc.pridict_ngrna import (
    build_pridict_ngrnas,
    nicking_rows_to_ngrnas,
    place_by_sequence,
)
from pegrna_pipeline.hpc.tjpe_workflow import derive_rc_pbs2

REPO = Path(__file__).resolve().parents[1]
INSERT_LEN = 24


def _clean(spacers):
    return {s: (1.0, 0, False) for s in spacers}


def _build(genome, site_s, rows, *, strand="Rv", arm_side="auto", lo=40, hi=90):
    return build_pridict_ngrnas(
        rows,
        genome=genome,
        site_s=site_s,
        pegrna_strand=strand,
        ngrna_distance_min=lo,
        ngrna_distance_max=hi,
        offtarget_scorer=_clean,
        arm_side=arm_side,
    )


# ---- locating a guide by its sequence ------------------------------------------------------------


def test_a_top_strand_guide_is_located_at_its_true_cut():
    locus = make_locus([("top", 100)], n=300, site_s=150, seed=2)
    _, _, spacer = locus.guides[0]
    assert locate_guide(locus.genome, spacer) == ("Fw", 100)


def test_a_bottom_strand_guide_is_located_at_its_true_cut():
    locus = make_locus([("bottom", 200)], n=300, site_s=150, seed=2)
    _, _, spacer = locus.guides[0]
    assert locate_guide(locus.genome, spacer) == ("Rv", 200)


def test_the_cut_follows_the_same_rule_on_both_strands():
    proto, pad = "GATCGATCGGCTAGCTTACG", "ACACACACACACACACACAC"
    top, bottom = pad + proto + "TGG" + pad, pad + rc(proto + "TGG") + pad
    assert locate_guide(top, proto) == ("Fw", len(pad) + 17)
    assert locate_guide(bottom, proto) == ("Rv", len(pad) + 6)


def test_the_spacer_may_be_given_in_lower_case():
    locus = make_locus([("top", 100)], n=300, site_s=150, seed=2)
    assert locate_guide(locus.genome, locus.guides[0][2].lower()) == ("Fw", 100)


def test_a_spacer_absent_from_the_window_is_not_located():
    locus = make_locus([("top", 100)], n=300, site_s=150, seed=2)
    assert locate_guide(locus.genome, "ACGTACGTACGTACGTACGT") is None


def test_a_spacer_that_occurs_twice_is_not_located():
    locus = make_locus([("top", 100)], n=300, site_s=150, seed=2)
    spacer = locus.guides[0][2]
    twice = locus.genome + "ACACACACAC" + spacer + "TGG" + "ACACAC"
    assert locate_guide(twice, spacer) is None


def test_a_copy_with_no_pam_is_not_a_second_site():
    """Cas9 cannot cut where there is no PAM, so a PAM-less copy does not make the guide ambiguous."""
    locus = make_locus([("top", 100)], n=300, site_s=150, seed=2)
    spacer = locus.guides[0][2]
    with_copy = locus.genome + "ACACACACAC" + spacer + "TAC" + "ACACAC"
    assert locate_guide(with_copy, spacer) == ("Fw", 100)


def test_a_spacer_with_no_pam_after_it_is_not_a_cas9_site():
    spacer = "GATCGATCGGCTAGCTTACG"
    assert locate_guide("ACAC" * 10 + spacer + "TAC" + "ACAC" * 10, spacer) is None
    assert locate_guide("ACAC" * 10 + spacer + "TGG" + "ACAC" * 10, spacer) is not None


@pytest.mark.parametrize("spacer", ["", "ACGT", "ACGTACGTACGTACGTACGN", "ACGTACGTACGTACGTACGTA"])
def test_a_malformed_spacer_is_not_located(spacer):
    assert locate_guide("ACAC" * 30, spacer) is None


def test_a_guide_at_the_very_start_or_end_of_the_window_is_found():
    proto = "GATCGATCGGCTAGCTTACG"
    assert locate_guide(proto + "TGG" + "ACAC" * 5, proto) == ("Fw", 17)
    assert locate_guide("ACAC" * 5 + rc(proto + "TGG"), proto) == ("Rv", 20 + 6)


# ---- the injected nick follows the sequence, whatever PRIDICT says ---------------------------------


def test_the_nick_is_the_true_cut_even_when_pridict_is_one_off_upstream():
    locus = make_locus([("top", 104)], n=300, site_s=150, seed=2)
    rows = [guide_row(locus, 0)]  # PRIDICT says -47; truth is -46
    assert rows[0]["Nicking-Position-to-edit"] == "-47"
    (inj,) = _build(locus.genome, 150, rows, arm_side="5prime")
    assert inj.nick_position == 104
    assert inj.distance_to_edit == -46, "the distance reported is the true one"


def test_the_nick_is_the_true_cut_downstream_where_pridict_adds_the_insert_length():
    locus = make_locus([("top", 205)], n=300, site_s=150, seed=2)
    rows = [guide_row(locus, 0)]  # PRIDICT says +78; truth is +55
    assert rows[0]["Nicking-Position-to-edit"] == "78"
    (inj,) = _build(locus.genome, 150, rows, arm_side="3prime")
    assert inj.nick_position == 205


def test_a_reverse_strand_guide_lands_on_its_own_side_of_s():
    locus = make_locus([("bottom", 205)], n=300, site_s=150, seed=2)
    rows = [guide_row(locus, 0)]  # PRIDICT says -55 (mirrored)
    assert rows[0]["Nicking-Position-to-edit"] == "-55"
    (inj,) = _build(locus.genome, 150, rows, strand="Fw", arm_side="3prime")
    assert inj.nick_position == 205 and inj.target_strand == "Rv"


def test_a_reverse_strand_guide_upstream_is_not_mirrored_to_the_other_side():
    locus = make_locus([("bottom", 95)], n=300, site_s=150, seed=2)
    rows = [guide_row(locus, 0)]  # PRIDICT says +79: it would read as 3' of S
    assert rows[0]["Nicking-Position-to-edit"] == str(55 + INSERT_LEN)
    (inj,) = _build(locus.genome, 150, rows, strand="Fw", arm_side="5prime")
    assert inj.nick_position == 95 and inj.distance_to_edit == -55


def test_the_distance_window_and_the_arm_side_use_the_true_position():
    """Whatever number the row carries, the guide is judged where it really cuts: 30 bp from S."""
    locus = make_locus([("top", 120)], n=300, site_s=150, seed=2)
    rows = [nicking_row(locus.guides[0][2], -45, "Fw")]  # a wrong number that would pass the window
    assert _build(locus.genome, 150, rows, lo=40, hi=70) == []
    assert _build(locus.genome, 150, rows, lo=25, hi=35)[0].nick_position == 120


def test_the_strand_is_the_one_the_guide_is_on_not_the_one_the_row_claims():
    locus = make_locus([("top", 104)], n=300, site_s=150, seed=2)
    spacer = locus.guides[0][2]
    claimed = nicking_rows_to_ngrnas([nicking_row(spacer, 5, "Rv")])  # wrong strand, wrong side
    (placed,) = place_by_sequence(claimed, genome=locus.genome, site_s=150)
    assert placed.target_strand == "Fw" and placed.distance_to_edit == -46


def test_a_guide_pridict_returned_that_is_not_in_the_window_is_dropped():
    locus = make_locus([("top", 104)], n=300, site_s=150, seed=2)
    rows = [nicking_row("ACGTACGTACGTACGTACGT", -46, "Fw")]
    assert _build(locus.genome, 150, rows) == []


def test_a_guide_found_twice_is_dropped_not_guessed():
    locus = make_locus([("top", 104)], n=300, site_s=150, seed=2)
    spacer = locus.guides[0][2]
    repeated = locus.genome + "ACACACAC" + spacer + "TGG" + "ACAC"
    assert _build(repeated, 150, [nicking_row(spacer, -47, "Fw")]) == []


def test_the_rcpbs2_is_taken_at_the_true_cut():
    locus = make_locus([("top", 104)], n=300, site_s=150, seed=2)
    (inj,) = _build(locus.genome, 150, [guide_row(locus, 0)], arm_side="5prime")
    assert inj.rcpbs2 == derive_rc_pbs2(locus.genome, ngrna_nick_position=104, pbs2_length=17)


def test_the_offtarget_scan_sees_only_guides_that_were_placed():
    locus = make_locus([("top", 104)], n=300, site_s=150, seed=2)
    seen = []

    def spy(spacers):
        seen.append(list(spacers))
        return _clean(spacers)

    rows = [guide_row(locus, 0), nicking_row("ACGTACGTACGTACGTACGT", -50, "Fw")]
    build_pridict_ngrnas(
        rows,
        genome=locus.genome,
        site_s=150,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=90,
        offtarget_scorer=spy,
    )
    assert seen == [[locus.guides[0][2]]]


# ---- the real CTCF window ----------------------------------------------------------------------


def _ctcf():
    cfg = json.loads((REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json").read_text())
    return cfg["genome"].upper(), cfg["site_s"]


def test_locate_guide_agrees_with_an_independent_enumeration_on_the_ctcf_window():
    """``find_guides`` (pa_geometry) and ``locate_guide`` (nick_search) are written separately, from the
    SpCas9 rule, and must agree on every guide that is unique in the window."""
    genome, _ = _ctcf()
    checked = 0
    for g in find_guides(genome):
        hits = [x for x in find_guides(genome) if x.spacer == g.spacer]
        if len(hits) != 1:
            continue
        assert locate_guide(genome, g.spacer) == ("Fw" if g.strand == "top" else "Rv", g.cut)
        checked += 1
    assert checked > 40, "the real window must exercise both strands"


def test_every_ctcf_guide_is_placed_at_its_true_cut_under_pridicts_observed_numbers():
    genome, site_s = _ctcf()
    checked, strands = 0, set()
    for g in find_guides(genome):
        label, strand_for_peg = ("Fw", "Rv") if g.strand == "top" else ("Rv", "Fw")
        rel = g.cut - site_s
        if abs(rel) < 25 or abs(rel) > 120:
            continue
        if sum(1 for x in find_guides(genome) if x.spacer == g.spacer) != 1:
            continue
        row = nicking_row(g.spacer, pridict_distance(label, rel, INSERT_LEN), label)
        placed = build_pridict_ngrnas(
            [row],
            genome=genome,
            site_s=site_s,
            pegrna_strand=strand_for_peg,
            ngrna_distance_min=25,
            ngrna_distance_max=120,
            offtarget_scorer=_clean,
        )
        if not placed:  # a spacer that fails the sequence rules, or a PBS2 window off the window
            continue
        assert placed[0].nick_position == g.cut, (g.strand, g.cut, rel)
        strands.add(g.strand)
        checked += 1
    assert checked >= 10 and strands == {"top", "bottom"}, (checked, strands)
