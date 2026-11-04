"""A guide's nick must be reported the same way on both strands (issue #31).

``enumerate_nick_sites`` reported a forward-strand guide's cut as the index of the first base RIGHT of the
cut (``genome[:b] | genome[b:]``) and a reverse-strand guide's as the last base LEFT of it. Everything
downstream reads the number as "first base right of the cut" (``_assemble_edit`` replaces
``genome[nick:S)``, ``_select_rcpbs2`` takes ``genome[nick-L:nick]``), so a forward-strand pegRNA, whose second
nick is a reverse-strand guide, got an arm one base short or one base long. The cuts below are computed
independently of the code under test, from the SpCas9 rule: the PAM strand is nicked 3 bp 5' of the PAM.
"""

import pytest
from pa_helpers import fake_scorers, make_cargo, make_locus, rc

from pegrna_pipeline.hpc.nick_search import enumerate_nick_sites
from pegrna_pipeline.hpc.recoded_tjpe_workflow import run_recoded_tjpe

PROTO = "GATCGATCGGCTAGCTTACG"
PAD = "ACACACACACACACACACAC"


def _true_cut_top(genome: str, protospacer: str) -> int:
    """First base right of the cut for a guide whose protospacer reads 5'->3' on the top strand."""
    return genome.index(protospacer) + 17  # the cut falls between protospacer bases 17 and 18


def _true_cut_bottom(genome: str, protospacer: str) -> int:
    """Same, for a guide on the bottom strand: its protospacer appears reverse-complemented in the sense
    window, the PAM complement ('CCN') sits just left of it, and the cut is 3 bp from that PAM."""
    return genome.index(rc(protospacer)) + 3


def _find(genome: str, strand_of_pegrna: str):
    return {
        c.protospacer: c
        for c in enumerate_nick_sites(
            genome,
            site_s=len(genome) // 2,
            pegrna_strand=strand_of_pegrna,
            distance_min=0,
            distance_max=len(genome),
        )
    }


def test_a_forward_strand_guide_is_reported_at_its_true_cut():
    genome = PAD + PROTO + "TGG" + PAD
    nick = _find(genome, "Rv")[PROTO].nick_position  # the second nick is opposite the pegRNA
    assert nick == _true_cut_top(genome, PROTO)


def test_a_reverse_strand_guide_is_reported_at_its_true_cut():
    genome = PAD + rc(PROTO + "TGG") + PAD
    nick = _find(genome, "Fw")[PROTO].nick_position
    assert nick == _true_cut_bottom(genome, PROTO)


@pytest.mark.parametrize("flank", [0, 1, 7, 33])
def test_the_convention_does_not_depend_on_where_in_the_window_the_guide_sits(flank):
    pad = "AC" * (10 + flank)
    top = pad + PROTO + "TGG" + pad
    bottom = pad + rc(PROTO + "TGG") + pad
    assert _find(top, "Rv")[PROTO].nick_position == _true_cut_top(top, PROTO)
    assert _find(bottom, "Fw")[PROTO].nick_position == _true_cut_bottom(bottom, PROTO)


def test_the_two_strands_agree_on_what_a_cut_is():
    """The same guide placed on either strand at the same cut gives the same number."""
    cut = 40
    top = make_locus([("top", cut)], n=120, site_s=60, seed=3)
    bottom = make_locus([("bottom", cut)], n=120, site_s=60, seed=3)
    (t,) = [c for c in _find(top.genome, "Rv").values()]
    (b,) = [c for c in _find(bottom.genome, "Fw").values()]
    assert t.nick_position == b.nick_position == cut


# ---- a forward-strand pegRNA design, checked against the true cut ---------------------------------


def _design(pegrna_strand: str, second_strand: str, cut: int):
    locus = make_locus([(second_strand, cut)], n=300, site_s=150, seed=5)
    cargo = make_cargo(60)
    result = run_recoded_tjpe(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa=None,
        cargo=cargo,
        genome=locus.genome,
        site_s=locus.site_s,
        cds_mask=locus.cds_mask,
        pbs_lengths=(13,),
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        pegrna_strand=pegrna_strand,
        arm_side="5prime",
        max_run=8,
        max_insert_variants=1,
        **fake_scorers(),
    )
    return locus, cargo, result


@pytest.mark.parametrize(
    "pegrna_strand, second_strand",
    [("Rv", "top"), ("Fw", "bottom")],
    ids=["Rv-pegRNA", "Fw-pegRNA"],
)
def test_the_arm_covers_exactly_the_bases_between_the_two_cuts(pegrna_strand, second_strand):
    cut = 100  # 50 bp 5' of S
    locus, cargo, result = _design(pegrna_strand, second_strand, cut)
    genome, s = locus.genome, locus.site_s
    assert result.nick_position == cut, "the second nick must be where the guide really cuts"
    # the product: untouched left of the second cut, the (recoded) arm, the cargo at S, untouched right of S
    assert result.edited_allele[:cut] == genome[:cut]
    assert result.edited_allele[s + len(cargo) :] == genome[s:]
    assert result.edited_allele[s : s + len(cargo)] == cargo
    assert len(result.edited_allele) == len(genome) + len(cargo), "no base lost or duplicated"
    arm = result.edited_allele[cut:s]
    assert len(arm) == s - cut
    assert result.deletion_size == 0
