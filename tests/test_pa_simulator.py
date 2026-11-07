"""``simulate_pa`` rebuilds a Prime Assembly outcome from what would be ordered and checks it against the intent.

The strand arithmetic of a two-pegRNA design is easy to get silently wrong (ADR 0004 records how it went wrong in
tj-PE), so, like ``simulate_tprt``, an executable check works from the two ORDERED pegRNAs and the donor and
reproduces the molecular mechanism rather than re-running the formulas that made them: the RT copies a pegRNA's
RTT 3'->5' onto the 3' end of the nicked strand, to which the PBS has annealed. If any strand operation (reverse,
complement, reverse complement) is wrong on either pegRNA, it must raise.

This file pins three things: the happy path, a property test over many random loci, and a mutation matrix in which
EVERY single wrong strand operation on either pegRNA's RTT or PBS is caught.
"""

import random

import pytest
from pa_helpers import make_cargo, make_locus, rc

from pegrna_pipeline.hpc.pa_geometry import (
    assemble_donor,
    assemble_pa_pegrna,
    build_arm,
    derive_extension,
    find_guides,
    footprints_clear,
    product_allele,
)
from pegrna_pipeline.hpc.pa_simulator import PAMismatch, simulate_pa, simulate_side
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

S = 150
SCAFFOLD = len(STANDARD_SGRNA_SCAFFOLD)


def _complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))


class Design:
    """One complete, valid design on a planted locus, with every piece the simulator is given."""

    def __init__(self, a_cut=123, b_cut=164, f_a=35, l_a=13, f_b=35, l_b=15, cargo_len=120, seed=7):
        self.locus = make_locus(
            [("top", a_cut), ("bottom", b_cut)], n=300, site_s=S, seed=seed, coding_end=S
        )
        self.genome = self.locus.genome
        guides = find_guides(self.genome)
        self.a = next(g for g in guides if g.strand == "top" and g.cut == a_cut)
        self.b = next(g for g in guides if g.strand == "bottom" and g.cut == b_cut)
        arm5 = build_arm(self.genome, a_cut, S, self.locus.cds_mask, 8)
        arm3 = build_arm(self.genome, S, b_cut, self.locus.cds_mask, 8)
        self.cargo = make_cargo(cargo_len)
        self.donor = assemble_donor(arm5.sense, self.cargo, arm3.sense)
        self.product = product_allele(self.genome, a_cut, b_cut, self.donor)
        self.f_a, self.l_a, self.f_b, self.l_b = f_a, l_a, f_b, l_b
        self.ext_a = derive_extension("A", self.genome, a_cut, self.donor, f_a, l_a)
        self.ext_b = derive_extension("B", self.genome, b_cut, self.donor, f_b, l_b)

    def peg(self, side, rtt=None, pbs=None, spacer=None):
        ext, guide = (self.ext_a, self.a) if side == "A" else (self.ext_b, self.b)
        return assemble_pa_pegrna(
            spacer=spacer or guide.spacer, rtt=rtt or ext.rtt, pbs=pbs or ext.pbs
        )

    def run(self, **over):
        args = dict(
            genome=self.genome,
            a=self.a,
            b=self.b,
            donor=self.donor,
            peg_a=self.peg("A"),
            peg_b=self.peg("B"),
            flap_a=self.f_a,
            pbs_a=self.l_a,
            flap_b=self.f_b,
            pbs_b=self.l_b,
        )
        args.update(over)
        return simulate_pa(**args)


@pytest.fixture(scope="module")
def d():
    return Design()


# ---- the happy path -------------------------------------------------------------------------------


def test_a_correct_design_reproduces_the_intended_product(d):
    assert d.run() == d.product


def test_the_product_is_the_donor_standing_where_the_span_was(d):
    out = d.run()
    assert out == d.genome[: d.a.cut] + d.donor + d.genome[d.b.cut :]
    assert len(out) == len(d.genome) + len(d.cargo), "net deletion is zero: only the cargo is added"


def test_the_expected_product_is_compared_when_given(d):
    assert d.run(expected_product=d.product) == d.product
    with pytest.raises(PAMismatch, match="expected product"):
        d.run(expected_product=d.product[:-1] + ("A" if d.product[-1] != "A" else "C"))


@pytest.mark.parametrize(
    "f_a, l_a, f_b, l_b", [(30, 13, 50, 17), (50, 17, 30, 13), (35, 8, 35, 20)]
)
def test_other_flap_and_pbs_lengths_work(f_a, l_a, f_b, l_b):
    des = Design(f_a=f_a, l_a=l_a, f_b=f_b, l_b=l_b, cargo_len=150)
    assert des.run() == des.product


def test_a_guide_whose_cut_is_at_s_leaves_no_arm_on_that_side():
    des = Design(a_cut=S - 0, b_cut=190, cargo_len=130)  # guide A cuts exactly at S: arm5 is empty
    assert des.run() == des.product


def test_simulate_side_returns_the_flap_in_sense_orientation(d):
    flap_a = simulate_side(
        side="A",
        genome=d.genome,
        guide=d.a,
        donor=d.donor,
        peg=d.peg("A"),
        flap_len=d.f_a,
        pbs_len=d.l_a,
    )
    flap_b = simulate_side(
        side="B",
        genome=d.genome,
        guide=d.b,
        donor=d.donor,
        peg=d.peg("B"),
        flap_len=d.f_b,
        pbs_len=d.l_b,
    )
    assert flap_a == d.donor[: d.f_a] and flap_b == d.donor[-d.f_b :]


# ---- a property test over many random loci --------------------------------------------------------


def _random_case(seed):
    """A random genome, a random PAM-in pair of its own guides, a random donor and random lengths."""
    rng = random.Random(seed)
    n = rng.randint(240, 420)
    genome = "".join(rng.choice("ACGT") for _ in range(n))
    guides = find_guides(genome)
    tops = [g for g in guides if g.strand == "top" and g.cut >= 22]
    bottoms = [g for g in guides if g.strand == "bottom" and g.cut <= n - 40]
    pairs = [
        (a, b)
        for a in tops
        for b in bottoms
        if footprints_clear(a, b) and 23 <= b.cut - a.cut <= 130
    ]
    if not pairs:
        return None
    a, b = rng.choice(pairs)
    s = rng.randint(a.cut, b.cut)
    cargo = "".join(rng.choice("ACGT") for _ in range(rng.randint(100, 220)))
    donor = genome[a.cut : s] + cargo + genome[s : b.cut]
    f_a = rng.randint(20, min(50, len(donor) // 2))
    f_b = rng.randint(20, min(50, len(donor) // 2))
    l_a, l_b = rng.randint(8, 20), rng.randint(8, 20)
    ext_a = derive_extension("A", genome, a.cut, donor, f_a, l_a)
    ext_b = derive_extension("B", genome, b.cut, donor, f_b, l_b)
    if ext_a is None or ext_b is None:
        return None
    return dict(
        genome=genome,
        a=a,
        b=b,
        donor=donor,
        ext_a=ext_a,
        ext_b=ext_b,
        peg_a=assemble_pa_pegrna(spacer=a.spacer, rtt=ext_a.rtt, pbs=ext_a.pbs),
        peg_b=assemble_pa_pegrna(spacer=b.spacer, rtt=ext_b.rtt, pbs=ext_b.pbs),
        f_a=f_a,
        f_b=f_b,
        l_a=l_a,
        l_b=l_b,
    )


def _simulate(case, **over):
    args = dict(
        genome=case["genome"],
        a=case["a"],
        b=case["b"],
        donor=case["donor"],
        peg_a=case["peg_a"],
        peg_b=case["peg_b"],
        flap_a=case["f_a"],
        pbs_a=case["l_a"],
        flap_b=case["f_b"],
        pbs_b=case["l_b"],
    )
    args.update(over)
    return simulate_pa(**args)


SEEDS = range(300)


def test_a_correct_design_always_reproduces_the_product_on_random_loci():
    ran = 0
    for seed in SEEDS:
        case = _random_case(seed)
        if case is None:
            continue
        expected = product_allele(case["genome"], case["a"].cut, case["b"].cut, case["donor"])
        assert _simulate(case) == expected, seed
        ran += 1
    assert ran >= 250, f"only {ran} of {len(SEEDS)} random loci produced a design"


def _mutations(side, case):
    """Every single wrong strand operation on one side's RTT or PBS: (label, part, wrong sequence)."""
    ext = case["ext_a"] if side == "A" else case["ext_b"]
    out = []
    for part in ("rtt", "pbs"):
        right = getattr(ext, part)
        # A string has exactly four strand operations (identity, reverse, complement, reverse complement) and
        # exactly one is right for each part, so applying the other three to the right answer is every way to
        # get it wrong.
        out += [
            (f"{side} {part} reversed", part, right[::-1]),
            (f"{side} {part} complemented", part, _complement(right)),
            (f"{side} {part} reverse-complemented once too often", part, rc(right)),
        ]
    return out


def test_every_wrong_strand_operation_on_either_pegrna_is_caught():
    checked = 0
    for seed in range(60):
        case = _random_case(seed)
        if case is None:
            continue
        for side in ("A", "B"):
            ext = case["ext_a"] if side == "A" else case["ext_b"]
            guide = case["a"] if side == "A" else case["b"]
            for label, part, wrong in _mutations(side, case):
                if wrong == getattr(ext, part):
                    continue  # a palindromic stretch: the "wrong" operation changed nothing
                peg = assemble_pa_pegrna(
                    spacer=guide.spacer,
                    rtt=wrong if part == "rtt" else ext.rtt,
                    pbs=wrong if part == "pbs" else ext.pbs,
                )
                with pytest.raises(PAMismatch):
                    _simulate(case, **{"peg_a" if side == "A" else "peg_b": peg})
                checked += 1
    assert checked >= 500, f"the mutation matrix only ran {checked} cases"


def test_the_mutation_that_matters_most_is_not_caught_by_a_weaker_check():
    """A reversed RTT is a real, plausible slip. It must name the strand, not just fail."""
    case = next(c for c in map(_random_case, SEEDS) if c is not None)
    bad = assemble_pa_pegrna(
        spacer=case["a"].spacer, rtt=case["ext_a"].rtt[::-1], pbs=case["ext_a"].pbs
    )
    with pytest.raises(PAMismatch, match="flap A"):
        _simulate(case, peg_a=bad)


# ---- everything else that can be wrong ------------------------------------------------------------


def test_pegrnas_given_in_the_wrong_order_are_caught(d):
    with pytest.raises(PAMismatch):
        d.run(peg_a=d.peg("B"), peg_b=d.peg("A"))


def test_a_wrong_spacer_is_caught(d):
    other = "GATCGATCGGCTAGCTTACG"
    assert other != d.a.spacer
    with pytest.raises(PAMismatch, match="spacer"):
        d.run(peg_a=d.peg("A", spacer=other))


@pytest.mark.parametrize("side", ["A", "B"])
def test_a_guide_whose_spacer_is_not_in_the_genome_is_caught_even_if_the_pegrna_uses_it(d, side):
    """The pegRNA and the guide agree with each other and disagree with the genome: only reading the genome
    at the protospacer can catch this (corrupted guide data, or a guide from another window)."""
    from dataclasses import replace

    guide = d.a if side == "A" else d.b
    other = "GATCGATCGGCTAGCTTACG"
    assert other != guide.spacer
    odd = replace(guide, spacer=other)
    args = dict(peg_a=d.peg("A"), peg_b=d.peg("B"))
    args["peg_a" if side == "A" else "peg_b"] = d.peg(side, spacer=other)
    args.update(a=odd if side == "A" else d.a, b=odd if side == "B" else d.b)
    with pytest.raises(PAMismatch, match="genome holds at its protospacer"):
        d.run(**args)


def test_a_guide_whose_cut_does_not_match_the_genome_is_caught(d):
    from dataclasses import replace

    with pytest.raises(PAMismatch, match="cut"):
        d.run(a=replace(d.a, cut=d.a.cut + 1))
    with pytest.raises(PAMismatch, match="cut"):
        d.run(b=replace(d.b, cut=d.b.cut - 1))


def test_a_pbs_built_for_a_neighbouring_nick_is_caught(d):
    shifted = derive_extension("A", d.genome, d.a.cut + 1, d.donor, d.f_a, d.l_a)
    with pytest.raises(PAMismatch, match="PBS"):
        d.run(peg_a=d.peg("A", pbs=shifted.pbs))
    shifted_b = derive_extension("B", d.genome, d.b.cut - 1, d.donor, d.f_b, d.l_b)
    with pytest.raises(PAMismatch, match="PBS"):
        d.run(peg_b=d.peg("B", pbs=shifted_b.pbs))


def test_a_guide_on_the_wrong_strand_is_caught(d):
    from dataclasses import replace

    with pytest.raises(PAMismatch, match="top"):
        d.run(a=replace(d.a, strand="bottom"))
    with pytest.raises(PAMismatch, match="bottom"):
        d.run(b=replace(d.b, strand="top"))


def test_a_before_b_is_required():
    """Two guides that are each consistent with the genome, but with the top-strand guide on the RIGHT: the flaps
    would grow away from each other (PAM-out), which never fills the gap."""
    locus = make_locus([("top", 123), ("bottom", 80)], n=300, site_s=S, seed=11, coding_end=S)
    guides = find_guides(locus.genome)
    a = next(g for g in guides if g.strand == "top" and g.cut == 123)
    b = next(g for g in guides if g.strand == "bottom" and g.cut == 80)
    donor = make_cargo(120)
    ext_a = derive_extension("A", locus.genome, a.cut, donor, 35, 13)
    ext_b = derive_extension("B", locus.genome, b.cut, donor, 35, 13)
    with pytest.raises(PAMismatch, match="left of"):
        simulate_pa(
            genome=locus.genome,
            a=a,
            b=b,
            donor=donor,
            peg_a=assemble_pa_pegrna(spacer=a.spacer, rtt=ext_a.rtt, pbs=ext_a.pbs),
            peg_b=assemble_pa_pegrna(spacer=b.spacer, rtt=ext_b.rtt, pbs=ext_b.pbs),
            flap_a=35,
            pbs_a=13,
            flap_b=35,
            pbs_b=13,
        )


def test_a_donor_whose_end_does_not_match_a_flap_is_caught(d):
    first = "A" if d.donor[0] != "A" else "C"
    with pytest.raises(PAMismatch, match="flap A"):
        d.run(donor=first + d.donor[1:])
    last = "A" if d.donor[-1] != "A" else "C"
    with pytest.raises(PAMismatch, match="flap B"):
        d.run(donor=d.donor[:-1] + last)


def test_flaps_that_would_overlap_each_other_are_refused():
    des = Design(cargo_len=120)
    short = des.donor[: des.f_a + des.f_b - 1]  # shorter than the two flaps together
    with pytest.raises(PAMismatch, match="overlap"):
        des.run(donor=short)


def test_a_pegrna_of_the_wrong_length_is_caught(d):
    with pytest.raises(PAMismatch, match="length"):
        d.run(flap_a=d.f_a + 1)
    with pytest.raises(PAMismatch, match="length"):
        d.run(peg_b=d.peg("B")[:-1])


def test_a_wrong_scaffold_is_caught(d):
    peg = d.peg("A")
    broken = peg[:20] + ("A" if peg[20] != "A" else "C") + peg[21:]
    with pytest.raises(PAMismatch, match="scaffold"):
        d.run(peg_a=broken)


def test_a_wrong_tail_is_caught(d):
    peg = d.peg("A")
    with pytest.raises(PAMismatch, match="tail"):
        d.run(peg_a=peg[:-3] + "AAA")
    with pytest.raises(PAMismatch, match="tail"):
        d.run(peg_a=peg[:-3] + "UUU", tail_len=2)


def test_a_custom_scaffold_and_tail_length_are_honoured():
    des = Design()
    custom = STANDARD_SGRNA_SCAFFOLD[:-1] + ("A" if STANDARD_SGRNA_SCAFFOLD[-1] != "A" else "C")
    pa = assemble_pa_pegrna(
        spacer=des.a.spacer, rtt=des.ext_a.rtt, pbs=des.ext_a.pbs, scaffold=custom, tail="UU"
    )
    pb = assemble_pa_pegrna(
        spacer=des.b.spacer, rtt=des.ext_b.rtt, pbs=des.ext_b.pbs, scaffold=custom, tail="UU"
    )
    out = des.run(peg_a=pa, peg_b=pb, scaffold=custom, tail_len=2)
    assert out == des.product
    with pytest.raises(PAMismatch, match="scaffold"):
        des.run(peg_a=pa, peg_b=pb, tail_len=2)  # the standard scaffold is expected by default


def test_a_guide_without_an_ngg_pam_is_caught(d):
    from dataclasses import replace

    with pytest.raises(PAMismatch, match="PAM"):
        d.run(a=replace(d.a, pam="TAC"))


@pytest.mark.parametrize("side", ["A", "B"])
def test_a_pam_that_is_not_ngg_is_refused_even_when_the_genome_agrees(d, side):
    """The genome-match check alone would pass this: the guide and the genome agree on a PAM of TAC."""
    from dataclasses import replace

    guide = d.a if side == "A" else d.b
    lo, hi = guide.pam_span
    sense = "TAC" if side == "A" else rc("TAC")  # a bottom guide's PAM is read reverse-complemented
    genome = d.genome[:lo] + sense + d.genome[hi:]
    odd = replace(guide, pam="TAC")
    args = dict(genome=genome, donor=d.donor, peg_a=d.peg("A"), peg_b=d.peg("B"))
    args.update(a=odd if side == "A" else d.a, b=odd if side == "B" else d.b)
    with pytest.raises(PAMismatch, match="NGG"):
        d.run(**args)


def test_the_simulator_never_changes_its_inputs(d):
    before = (d.genome, d.donor, d.peg("A"), d.peg("B"))
    d.run()
    assert before == (d.genome, d.donor, d.peg("A"), d.peg("B"))
