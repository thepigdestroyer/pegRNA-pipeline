"""Pure geometry of a Prime Assembly design (docs/STRUCTURE.md, spec 4.1). No models.

Conventions under test: everything 5'->3'; a cut ``b`` is the index of the first base right of the cut,
``genome[:b] | genome[b:]``, for a guide on either strand; guide A is on the top strand left of guide B on the
bottom strand (PAM-in), so both PAMs lie between the nicks.
"""

import pytest
from Bio.Seq import Seq
from pa_helpers import make_cargo, make_locus, rc

from pegrna_pipeline.hpc.microhomology import longest_common_run
from pegrna_pipeline.hpc.pa_geometry import (
    MIN_NICK_SEPARATION_BP,
    Arm,
    Extension,
    Guide,
    GuidePair,
    PamOutcome,
    assemble_donor,
    assemble_pa_pegrna,
    break_pams,
    build_arm,
    choose_extension,
    derive_extension,
    find_guides,
    flap_lengths,
    footprints_clear,
    genomic_run_under_flap,
    pa_regions,
    pam_in_pairs,
    product_allele,
)
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

S = 150
TOPS = (123, 147)  # 5' guides: arms of 27 and 3 bp
BOTTOMS = (164, 190)  # 3' guides: arms of 14 and 40 bp


@pytest.fixture(scope="module")
def locus():
    plants = [("top", c) for c in TOPS] + [("bottom", c) for c in BOTTOMS]
    return make_locus(plants, n=300, site_s=S, seed=7)


# ---- guides and cuts ------------------------------------------------------------------------


def test_the_planted_guides_are_exactly_the_guides_found(locus):
    found = {(g.strand, g.cut, g.spacer) for g in find_guides(locus.genome)}
    assert found == set(locus.guides)


def test_every_guide_is_self_consistent_on_both_strands(locus):
    for g in find_guides(locus.genome):
        assert g.pam[1:] == "GG"
        if g.strand == "top":
            assert locus.genome[slice(*g.protospacer_span)] == g.spacer
            assert locus.genome[slice(*g.pam_span)] == g.pam
            assert g.cut == g.protospacer_span[0] + 17 == g.pam_span[0] - 3
        else:
            assert rc(locus.genome[slice(*g.protospacer_span)]) == g.spacer
            assert rc(locus.genome[slice(*g.pam_span)]) == g.pam
            assert g.cut == g.protospacer_span[0] + 3 == g.pam_span[1] + 3


def test_the_cut_is_the_same_kind_of_number_on_both_strands():
    proto, pad = "GATCGATCGGCTAGCTTACG", "ACACACACACACACACACAC"
    top = pad + proto + "TGG" + pad
    bottom = pad + rc(proto + "TGG") + pad
    (t,) = [g for g in find_guides(top) if g.spacer == proto]
    (b,) = [g for g in find_guides(bottom) if g.spacer == proto]
    assert t.cut == len(pad) + 17 and t.strand == "top"
    assert b.cut == len(pad) + 6 and b.strand == "bottom"


def test_guides_come_back_in_a_stable_order(locus):
    guides = find_guides(locus.genome)
    assert guides == sorted(guides, key=lambda g: (g.cut, g.strand, g.spacer))
    assert guides == find_guides(locus.genome)


def test_a_window_with_an_ambiguous_base_in_a_protospacer_yields_no_guide():
    proto, pad = "GATCGATCGGCTAGCTTACG", "ACACACACACACACACACAC"
    with_n = pad + proto[:5] + "N" + proto[6:] + "TGG" + pad
    assert not [g for g in find_guides(with_n) if g.strand == "top"]


# ---- PAM-in pairs ----------------------------------------------------------------------------


def test_pairs_are_top_left_and_bottom_right_with_s_between(locus):
    pairs = pam_in_pairs(find_guides(locus.genome), site_s=S, min_span=30, max_span=100)
    assert {(p.a.cut, p.b.cut) for p in pairs} == {(123, 164), (123, 190), (147, 190)}
    for p in pairs:
        assert p.a.strand == "top" and p.b.strand == "bottom"
        assert p.a.cut <= S <= p.b.cut
        assert p.span == p.b.cut - p.a.cut and 30 <= p.span <= 100


def test_a_bottom_guide_on_the_left_and_a_top_guide_on_the_right_make_no_pair():
    """PAM-out: the flaps would grow away from each other and never fill the gap."""
    left_bottom = make_locus([("bottom", 120), ("top", 180)], n=300, site_s=150, seed=4)
    assert (
        pam_in_pairs(find_guides(left_bottom.genome), site_s=150, min_span=30, max_span=100) == []
    )


def test_a_guide_on_the_wrong_side_of_s_is_excluded():
    both_left = make_locus([("top", 60), ("bottom", 110)], n=300, site_s=150, seed=4)
    assert pam_in_pairs(find_guides(both_left.genome), site_s=150, min_span=30, max_span=100) == []


def test_a_span_below_the_project_floor_is_refused():
    assert MIN_NICK_SEPARATION_BP == 23
    with pytest.raises(ValueError, match="at least 23"):
        pam_in_pairs([], site_s=S, min_span=22, max_span=100)


def _guide(strand, cut):
    if strand == "top":
        return Guide("top", "A" * 20, "TGG", cut, (cut + 3, cut + 6), (cut - 17, cut + 3))
    return Guide("bottom", "A" * 20, "AGG", cut, (cut - 6, cut - 3), (cut - 3, cut + 17))


def test_footprints_clear_exactly_when_the_two_pams_do_not_collide():
    top = _guide("top", 100)  # footprint [83, 106)
    assert not footprints_clear(top, _guide("bottom", 110))  # [104, 127): the PAMs overlap
    assert not footprints_clear(top, _guide("bottom", 111))  # [105, 128)
    assert footprints_clear(top, _guide("bottom", 112))  # [106, 129): span 12, touching
    assert footprints_clear(top, _guide("bottom", 140))


def test_pairs_are_returned_smallest_span_first(locus):
    pairs = pam_in_pairs(find_guides(locus.genome), site_s=S, min_span=30, max_span=100)
    assert [p.span for p in pairs] == sorted(p.span for p in pairs)
    assert isinstance(pairs[0], GuidePair)


# ---- arms ------------------------------------------------------------------------------------


def test_a_coding_arm_keeps_the_protein_and_breaks_long_identity(locus):
    arm = build_arm(locus.genome, 123, S, locus.cds_mask, 8)
    assert isinstance(arm, Arm) and len(arm.sense) == 27
    assert str(Seq(arm.sense).translate()) == str(Seq(arm.original).translate())
    assert longest_common_run(arm.sense, arm.original) <= 8
    assert arm.n_edits > 0 and arm.recode_span == (123, S) and arm.locked_span is None


def test_a_non_coding_arm_is_byte_exact_and_flagged(locus):
    arm = build_arm(locus.genome, S, 164, locus.cds_mask, 8)
    assert arm.sense == locus.genome[S:164]
    assert arm.locked_span == (S, 164) and arm.recode_span is None
    assert arm.n_edits == 0 and arm.microhomology_flag


def test_an_arm_that_crosses_the_coding_boundary_gets_both_treatments(locus):
    arm = build_arm(locus.genome, 141, 159, locus.cds_mask, 8)
    assert arm.sense[9:] == locus.genome[S:159]  # the non-coding part is untouched
    assert arm.recode_span == (141, S) and arm.locked_span == (S, 159)


def test_an_empty_arm_is_empty(locus):
    arm = build_arm(locus.genome, S, S, locus.cds_mask, 8)
    assert arm.sense == "" and arm.n_edits == 0 and arm.lo == arm.hi == S


def test_building_an_arm_twice_gives_the_same_arm(locus):
    assert build_arm(locus.genome, 123, S, locus.cds_mask, 8) == build_arm(
        locus.genome, 123, S, locus.cds_mask, 8
    )


# ---- PAM breaking, on crafted codons ---------------------------------------------------------
#
# codon:  [27,30)=GCT [30,33)=CTG [33,36)=GCT ...   the top guide's PAM is TGG at 31..33, so its
# critical bases are 32 (third base of the Leu codon CTG, which has synonyms) and 33.

CRAFTED = "GCT" * 10 + "CTG" + "GCT" * 9 + "GCT" * 20  # 120 nt, coding up to S = 60
S_CRAFTED = 60
MASK_CRAFTED = [i < S_CRAFTED for i in range(len(CRAFTED))]


def _crafted_arm(lo, hi, genome=CRAFTED, coding=True, sense=None):
    seq = sense if sense is not None else genome[lo:hi]
    return Arm(
        lo=lo,
        hi=hi,
        sense=seq,
        original=genome[lo:hi],
        recode_span=(lo, hi) if coding else None,
        locked_span=None if coding else (lo, hi),
        locked_expected=None if coding else genome[lo:hi],
        n_edits=sum(a != b for a, b in zip(seq, genome[lo:hi])),
        max_run_after=longest_common_run(seq, genome[lo:hi]) if seq else 0,
        microhomology_flag=not coding,
    )


def _top_guide(pam_start):
    return Guide(
        "top",
        "A" * 20,
        "TGG",
        pam_start - 3,
        (pam_start, pam_start + 3),
        (pam_start - 20, pam_start),
    )


def test_a_pam_in_a_coding_arm_is_broken_by_a_synonymous_change():
    assert CRAFTED[31:34] == "TGG" and CRAFTED[30:33] == "CTG"
    arm5, arm3 = _crafted_arm(28, S_CRAFTED), _crafted_arm(S_CRAFTED, 60, coding=False)
    new5, new3, outcomes = break_pams(
        CRAFTED, arm5, arm3, [("A", _top_guide(31))], MASK_CRAFTED, max_run=8
    )
    assert outcomes == [PamOutcome("A", "broken", 32)]
    changed = [i for i, (a, b) in enumerate(zip(arm5.sense, new5.sense)) if a != b]
    assert changed == [32 - 28], "exactly one base moves, and it is the PAM's first G"
    assert str(Seq(new5.sense[2:5]).translate()) == "L"  # CTG -> CTA/CTC/CTT, still leucine
    assert not (new5.sense[32 - 28] == "G" and new5.sense[33 - 28] == "G")
    assert new3 == arm3 and new5.n_edits == arm5.n_edits + 1


def test_a_pam_the_recode_already_broke_is_left_alone():
    broken = CRAFTED[:32] + "T" + CRAFTED[33:]  # CTG -> CTT
    arm5 = _crafted_arm(28, S_CRAFTED, sense=broken[28:S_CRAFTED])
    arm3 = _crafted_arm(S_CRAFTED, 60, coding=False)
    new5, _, outcomes = break_pams(CRAFTED, arm5, arm3, [("A", _top_guide(31))], MASK_CRAFTED, 8)
    assert outcomes == [PamOutcome("A", "broken-by-recode", None)]
    assert new5 == arm5


def test_a_pam_in_a_non_coding_arm_cannot_be_changed():
    genome = "GCT" * 20 + "AAATGGAAAC" + "GCT" * 10
    mask = [i < 60 for i in range(len(genome))]
    assert genome[63:66] == "TGG"
    arm5, arm3 = _crafted_arm(60, 60, genome), _crafted_arm(60, 90, genome, coding=False)
    guide = _top_guide(63)  # cut 60 = S: its PAM lies in the locked 3' arm
    _, new3, outcomes = break_pams(genome, arm5, arm3, [("A", guide)], mask, 8)
    assert outcomes == [PamOutcome("A", "intact-noncoding", None)]
    assert new3.sense == arm3.sense, "non-coding DNA is byte-exact (Decision 15)"


def test_a_pam_with_no_synonymous_alternative_is_reported():
    """TGG at a codon boundary is tryptophan: nothing synonymous exists."""
    genome = "GCT" * 11 + "TGG" + "GCT" * 20
    mask = [i < 60 for i in range(len(genome))]
    arm5, arm3 = _crafted_arm(30, 60, genome), _crafted_arm(60, 60, genome, coding=False)
    assert genome[33:36] == "TGG"
    new5, _, outcomes = break_pams(genome, arm5, arm3, [("A", _top_guide(33))], mask, 8)
    assert outcomes == [PamOutcome("A", "intact-no-synonymous-change", None)]
    assert new5.sense == arm5.sense


def test_a_bottom_guides_pam_is_broken_at_its_cc():
    """Sense GCC|CTG: the bottom PAM complement is CCT (26..28); the third base of Ala GCC is free."""
    genome = "GCT" * 8 + "GCC" + "CTG" + "GCT" * 25
    mask = [i < 60 for i in range(len(genome))]
    pam_start = 26
    assert genome[pam_start : pam_start + 3] == "CCT"
    guide = Guide(
        "bottom",
        "A" * 20,
        "AGG",
        pam_start + 6,
        (pam_start, pam_start + 3),
        (pam_start + 3, pam_start + 23),
    )
    arm5, arm3 = _crafted_arm(21, 60, genome), _crafted_arm(60, 60, genome, coding=False)
    new5, _, outcomes = break_pams(genome, arm5, arm3, [("B", guide)], mask, 8)
    assert outcomes[0].status == "broken" and outcomes[0].side == "B"
    assert new5.sense[pam_start - 21 : pam_start - 21 + 2] != "CC"
    assert str(Seq(new5.sense).translate()) == str(Seq(arm5.sense).translate())


# a coding stretch holding BOTH guides' PAMs (the CTCF best pair has B's PAM inside arm5):
#   27..29 GCT | 30..32 CTG | 33..35 GCT -> top PAM TGG at 31..33, critical 32, 33
#   63..65 GCC | 66..68 CTG               -> bottom PAM complement CCT at 65..67, critical 65, 66
BOTH = "GCT" * 10 + "CTG" + "GCT" * 10 + "GCC" + "CTG" + "GCT" * 27  # 150 nt, coding up to S = 90
S_BOTH = 90
MASK_BOTH = [i < S_BOTH for i in range(len(BOTH))]


def test_both_guides_pams_in_the_same_arm_are_each_broken():
    assert BOTH[31:34] == "TGG" and BOTH[65:68] == "CCT"
    bottom = Guide("bottom", "A" * 20, "AGG", 71, (65, 68), (68, 88))
    arm5, arm3 = _crafted_arm(27, S_BOTH, BOTH), _crafted_arm(S_BOTH, S_BOTH, BOTH, coding=False)
    new5, _, outcomes = break_pams(
        BOTH, arm5, arm3, [("A", _top_guide(31)), ("B", bottom)], MASK_BOTH, max_run=8
    )
    assert {(o.side, o.status) for o in outcomes} == {("A", "broken"), ("B", "broken")}
    sense = new5.sense  # index i is genome position 27 + i
    assert not (sense[32 - 27] == "G" and sense[33 - 27] == "G"), "guide A's PAM is gone"
    assert not (sense[65 - 27] == "C" and sense[66 - 27] == "C"), "guide B's PAM is gone"
    assert str(Seq(new5.sense).translate()) == str(Seq(arm5.sense).translate())
    assert (
        sum(a != b for a, b in zip(arm5.sense, new5.sense)) == 2
    ), "one base each, nothing else moves"


def test_breaking_a_pam_never_lets_a_recoded_arm_exceed_max_run():
    arm5 = build_arm(CRAFTED, 27, S_CRAFTED, MASK_CRAFTED, 8)
    assert arm5 is not None
    arm3 = _crafted_arm(S_CRAFTED, S_CRAFTED, coding=False)
    new5, _, outcomes = break_pams(CRAFTED, arm5, arm3, [("A", _top_guide(31))], MASK_CRAFTED, 8)
    assert outcomes[0].status in {"broken", "broken-by-recode"}
    assert longest_common_run(new5.sense, new5.original) <= 8


# ---- the donor and the product --------------------------------------------------------------


def test_the_donor_is_the_two_arms_around_the_cargo_and_the_product_restores_the_span(locus):
    cargo = make_cargo(120)
    a5 = build_arm(locus.genome, 123, S, locus.cds_mask, 8)
    a3 = build_arm(locus.genome, S, 164, locus.cds_mask, 8)
    donor = assemble_donor(a5.sense, cargo, a3.sense)
    assert donor == a5.sense + cargo + a3.sense
    product = product_allele(locus.genome, 123, 164, donor)
    assert len(product) == len(locus.genome) + len(cargo), "net deletion is zero"
    assert product[:123] == locus.genome[:123] and product[123 + len(donor) :] == locus.genome[164:]
    assert product[123 + len(a5.sense) : 123 + len(a5.sense) + len(cargo)] == cargo


def test_a_donor_part_that_is_not_dna_is_refused():
    with pytest.raises(ValueError):
        assemble_donor("ACGU", "ACGT", "ACGT")


# ---- flaps, extensions, the pegRNA ------------------------------------------------------------


def test_flap_lengths_scan_outward_from_the_default():
    assert flap_lengths(35, 30, 50)[:11] == [35, 36, 34, 37, 33, 38, 32, 39, 31, 40, 30]
    assert flap_lengths(35, 30, 50)[11:] == list(range(41, 51))
    assert len(set(flap_lengths(35, 30, 50))) == 21
    assert flap_lengths(30, 30, 32) == [30, 31, 32]


def _locus_and_donor(locus):
    cargo = make_cargo(120)
    a5 = build_arm(locus.genome, 123, S, locus.cds_mask, 8)
    a3 = build_arm(locus.genome, S, 164, locus.cds_mask, 8)
    return assemble_donor(a5.sense, cargo, a3.sense)


def test_guide_a_derivation_pairs_with_the_top_strand_and_writes_the_donors_start(locus):
    donor = _locus_and_donor(locus)
    ext = derive_extension("A", locus.genome, 123, donor, flap_len=35, pbs_len=13)
    assert isinstance(ext, Extension)
    assert rc(ext.pbs) == locus.genome[123 - 13 : 123]
    assert rc(ext.rtt) == donor[:35]
    assert (ext.flap_len, ext.pbs_len) == (35, 13) and len(ext.rtt) == 35


def test_guide_b_derivation_pairs_with_the_bottom_strand_and_writes_the_donors_end(locus):
    donor = _locus_and_donor(locus)
    ext = derive_extension("B", locus.genome, 164, donor, flap_len=35, pbs_len=13)
    assert ext.pbs == locus.genome[164 : 164 + 13]
    assert ext.rtt == donor[-35:]


def test_an_extension_that_leaves_the_genome_or_outruns_the_donor_is_none(locus):
    donor = _locus_and_donor(locus)
    assert (
        derive_extension("A", locus.genome, 5, donor, 35, 13) is None
    )  # PBS window off the left edge
    assert derive_extension("B", locus.genome, 295, donor, 35, 13) is None  # off the right edge
    assert (
        derive_extension("A", locus.genome, 123, "ACGT" * 5, 35, 13) is None
    )  # flap longer than donor


def test_the_sides_must_be_a_or_b(locus):
    with pytest.raises(ValueError):
        derive_extension("C", locus.genome, 123, "ACGT" * 30, 35, 13)


def _donor_with(first=None, last=None, filler="ACGTAGCATCAA"):
    core = (filler * 20)[:200]
    if first:
        core = first + core[len(first) :]
    if last:
        core = core[: -len(last)] + last
    return core


def test_the_first_flap_length_with_a_clean_non_c_start_is_chosen(locus):
    donor = _donor_with()
    ext = choose_extension("A", locus.genome, 123, donor, pbs_len=13)
    assert ext.flap_len == 35 and "TTTT" not in ext.rtt and ext.rtt[0] != "C"


def test_a_c_start_moves_to_the_next_flap_length():
    genome = "ACGTAGCATC" * 30
    donor = list(_donor_with())
    donor[34] = "G"  # flap 35's RTT starts with complement(donor[34]) = C
    donor[35] = "A"  # flap 36's starts with T
    ext = choose_extension("A", genome, 123, "".join(donor), pbs_len=13)
    assert ext.flap_len == 36 and ext.rtt[0] == "T"


def test_when_every_length_starts_with_c_the_first_clean_one_is_still_used():
    genome = "ACGTAGCATC" * 30
    donor = (
        "G" * 80
    )  # every RTT is all C: nothing better exists, so the first terminator-free one wins
    ext = choose_extension("A", genome, 123, donor, pbs_len=13)
    assert ext is not None and ext.flap_len == 35 and ext.rtt[0] == "C"


def test_a_terminator_in_the_rtt_rejects_every_flap_length():
    genome = "ACGTAGCATC" * 30
    donor = "AAAA" + "ACGTAGCATC" * 12  # RC(donor[:F]) contains TTTT for every F >= 30
    assert choose_extension("A", genome, 123, donor, pbs_len=13) is None


def test_a_terminator_in_the_pbs_rejects_the_extension():
    genome = (
        "ACGTAGCATC" * 12 + "AAAAAAAAAAAAA" + "ACGTAGCATC" * 12
    )  # RC of the PBS window is TTTT-rich
    donor = _donor_with()
    assert choose_extension("A", genome, 133, donor, pbs_len=13) is None


def test_the_genomic_run_under_each_flap(locus):
    a5 = build_arm(locus.genome, 123, S, locus.cds_mask, 8)
    a3 = build_arm(locus.genome, S, 164, locus.cds_mask, 8)
    cargo = make_cargo(120)
    donor = assemble_donor(a5.sense, cargo, a3.sense)
    run_a = genomic_run_under_flap("A", locus.genome, 123, 164, donor, 35)
    run_b = genomic_run_under_flap("B", locus.genome, 123, 164, donor, 35)
    assert run_a <= 8, "the coding arm was recoded"
    assert run_b >= len(a3.sense) - 1, "the non-coding arm is a copy of the genome under flap B"


def test_the_pegrna_is_spacer_scaffold_rtt_pbs_tail(locus):
    donor = _locus_and_donor(locus)
    ext = derive_extension("A", locus.genome, 123, donor, 35, 13)
    spacer = next(g for g in find_guides(locus.genome) if g.cut == 123).spacer
    peg = assemble_pa_pegrna(spacer=spacer, rtt=ext.rtt, pbs=ext.pbs)
    assert peg == spacer + STANDARD_SGRNA_SCAFFOLD + ext.rtt + ext.pbs + "UUU"
    assert len(peg) == 20 + 76 + 35 + 13 + 3
    regions = pa_regions(spacer, ext.rtt, ext.pbs)
    assert regions == {
        "spacer": (0, 20),
        "scaffold": (20, 96),
        "rtt": (96, 131),
        "pbs1": (131, 144),
    }


def test_the_longest_possible_pegrna_is_inside_the_calibrated_docking_range():
    assert 20 + len(STANDARD_SGRNA_SCAFFOLD) + 50 + 17 + 3 == 166 <= 232


def test_a_pegrna_part_that_is_not_dna_is_refused():
    with pytest.raises(ValueError):
        assemble_pa_pegrna(spacer="A" * 20, rtt="ACGU", pbs="ACGT")
    with pytest.raises(ValueError):
        assemble_pa_pegrna(spacer="A" * 20, rtt="ACGT", pbs="ACGT", tail="AAA")
