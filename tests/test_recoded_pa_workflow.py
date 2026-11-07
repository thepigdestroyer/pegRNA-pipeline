"""The Prime Assembly walk: a ranked search over PAM-in guide pairs (spec 4.2).

Fake scorers throughout (CI has no models); the real-model run is the acceptance test in
``scripts/analysis/acceptance_pa.py``. What is pinned here is the logic: which pairs exist, how they rank, what each
gate rejects and says, that nothing expensive runs twice, that the donor really restores the DNA between the nicks,
and that the answer is the same in every process.
"""

import dataclasses
import json
import subprocess
import sys
from pathlib import Path

import pytest
from pa_helpers import AllPassImportance, clean_offtarget, fake_scorers, make_cargo, make_locus, rc

from pegrna_pipeline.hpc.complex_dock import DockFake, DockingServiceError
from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.pa_geometry import find_guides, product_allele
from pegrna_pipeline.hpc.pa_simulator import PAMismatch, simulate_pa
from pegrna_pipeline.hpc.recoded_pa_workflow import (
    STANDING_CAVEAT,
    NoValidPADesign,
    RecodedPAResult,
    plan_pairs,
    rank_key,
    run_recoded_pa,
)
from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics

REPO = Path(__file__).resolve().parents[1]
S = 150
CARGO = make_cargo(120)

# tops: 123 (PAM AGG: Arg, whose only synonymous single change leaves NGA, so its PAM is weakened unless the
# recoder happened to break it) and 147; bottoms: 164 and 190. Coding up to S, so the
# 5' arm is recoded and the 3' arm, being non-coding, is copied exactly.
PLANTS = [("top", 123, "AGG"), ("top", 147), ("bottom", 164), ("bottom", 190)]


@pytest.fixture(scope="module")
def locus():
    return make_locus(PLANTS, n=300, site_s=S, seed=7, coding_end=S)


@pytest.fixture(scope="module")
def coding_locus():
    """Coding on BOTH sides of S: the 3' arm is recoded too, and a bottom guide's PAM can be broken."""
    return make_locus(
        [("top", 123, "AGG"), ("bottom", 164, "AGG")], n=300, site_s=S, seed=5, coding_end=300
    )


def inputs(loc, cargo=None, **over):
    base = dict(
        cargo=cargo or CARGO,
        genome=loc.genome,
        site_s=loc.site_s,
        cds_mask=loc.cds_mask,
        **fake_scorers(),
    )
    base.update(over)
    return base


def run(loc, cargo=None, **over):
    return run_recoded_pa(**inputs(loc, cargo, **over))


# ---- stand-in scorers ------------------------------------------------------------------------------


class RecordingDock(DockFake):
    def __init__(self, candidate_iptm=0.82):
        super().__init__(apo_iptm=0.85, candidate_iptm=candidate_iptm)
        self.shown = []

    def interface_iptm(self, pegrna, cas9):
        self.shown.append(pegrna)
        return super().interface_iptm(pegrna, cas9)


class RecordingFold(FullStructureFake):
    def __init__(self, pairs=None):
        super().__init__(pairs=pairs or {})
        self.seen = []

    def pair_probs(self, seq):
        self.seen.append(seq)
        return super().pair_probs(seq)


class CrossRegionFold:
    """Always calls one confident pair between the RTT and the PBS, whatever the pegRNA is."""

    def pair_probs(self, seq):
        return {(100, len(seq) - 6): 0.999}


class ScaffoldHuggingFold:
    """Pairs 16 scaffold nucleotides with the extension: no cross-region pair, so only the advisory fires."""

    def pair_probs(self, seq):
        return {(20 + k, 100 + k): 0.999 for k in range(16)}


class FlagBlockedImportance:
    def metrics(self, i):  # noqa: ARG002
        return ResidueMetrics(
            plddt=95.0,
            rsa=0.02,
            contact_number=30,
            is_catalytic=True,
            is_disulfide=False,
            is_positive_phi_gly=False,
            is_structural_pro=False,
        )


class AdvisorySplice:
    """A SpliceAI stand-in that says it had no real context (ADR 0028): reported, never rejecting."""

    binding = False
    context_note = "no genomic context: the design window padded with N"

    def __init__(self, delta):
        self._d = delta

    def delta(self, ref_seq, alt_seq):  # noqa: ARG002
        return self._d


def _dirty_for(spacers_to_fail):
    def scorer(protospacers):
        return {
            p: ((0.8, 4, True) if p in spacers_to_fail else (0.8, 0, False)) for p in protospacers
        }

    return scorer


# ---- refusals and argument checks ---------------------------------------------------------------------


def test_a_cargo_shorter_than_two_flaps_is_refused_with_the_way_out(locus):
    with pytest.raises(NoValidPADesign) as exc:
        run(locus, cargo=CARGO[:69])
    text = str(exc.value)
    assert "design-recoded" in text and "70" in text and "69" in text


def test_a_cargo_of_exactly_two_flaps_is_allowed(locus):
    assert len(CARGO[:70]) == 70 and run(locus, cargo=CARGO[:70]).donor


def test_the_refusal_is_a_no_valid_design_so_callers_that_catch_the_family_catch_it(locus):
    assert issubclass(NoValidPADesign, NoValidPegRNA)
    with pytest.raises(NoValidPegRNA):
        run(locus, cargo=make_cargo(40))


@pytest.mark.parametrize(
    "over, message",
    [
        (dict(flap_length=29), "flap_length"),
        (dict(flap_length=51), "flap_length"),
        (dict(flap_min=40, flap_max=30), "flap"),
        (dict(pbs_lengths=()), "pbs_lengths"),
        (dict(pbs_lengths=(13, 2)), "pbs_lengths"),
        (dict(min_span=22), "at least 23"),
        (dict(min_span=60, max_span=50), "max_span"),
        (dict(site_s=-1), "site_s"),
        (dict(site_s=301), "site_s"),
        (dict(tail="AAA"), "tail"),
        (dict(max_offtargets=-1), "max_offtargets"),
        (dict(max_pairs=0), "max_pairs"),
        (dict(max_run=0), "max_run"),
    ],
)
def test_bad_arguments_are_refused_with_the_name_of_the_argument(locus, over, message):
    with pytest.raises(ValueError, match=message):
        run(locus, **over)


def test_a_mask_that_does_not_fit_the_window_is_refused(locus):
    with pytest.raises(ValueError, match="cds_mask"):
        run(locus, cds_mask=locus.cds_mask[:-1])


def test_a_cargo_that_is_not_dna_is_refused(locus):
    with pytest.raises(ValueError, match="cargo"):
        run(locus, cargo="ACGT" * 20 + "N")


# ---- the happy path: the donor restores the DNA between the nicks ----------------------------------------


@pytest.fixture(scope="module")
def won(locus):
    return run(locus)


def test_a_design_is_found_and_has_the_documented_shape(won):
    assert isinstance(won, RecodedPAResult)
    assert won.guide_a.side == "A" and won.guide_b.side == "B"
    assert won.guide_a.strand == "top" and won.guide_b.strand == "bottom"
    assert won.guide_a.cut <= won.site_s <= won.guide_b.cut


def test_the_donor_is_arm5_cargo_arm3_and_the_product_is_the_donor_in_place_of_the_span(won, locus):
    assert won.donor == won.arm5 + won.cargo + won.arm3
    assert won.product == product_allele(locus.genome, *won.span, won.donor)
    assert [name for name, _ in won.donor_parts] == [
        n for n, seq in (("arm5", won.arm5), ("cargo", won.cargo), ("arm3", won.arm3)) if seq
    ]
    assert "".join(seq for _, seq in won.donor_parts) == won.donor


def test_nothing_is_lost_between_the_nicks_only_the_cargo_is_added(won, locus):
    """The owner's requirement (A+): the donor restores the DNA between the nicks. Net deletion is zero."""
    assert won.deletion_size == 0
    assert len(won.product) == len(locus.genome) + len(won.cargo)
    n1, n2 = won.span
    assert (
        won.product[:n1] == locus.genome[:n1]
        and won.product[n1 + len(won.donor) :] == locus.genome[n2:]
    )
    # remove the cargo again and what is left is the original window, apart from the silent recoding
    s = won.site_s
    without_cargo = won.product[:s] + won.product[s + len(won.cargo) :]
    diffs = [i for i, (x, y) in enumerate(zip(without_cargo, locus.genome)) if x != y]
    assert len(without_cargo) == len(locus.genome)
    assert len(diffs) == won.recode_edits and all(n1 <= i < n2 for i in diffs)


def test_the_non_coding_arm_is_copied_exactly(won, locus):
    """Decision 15: non-coding DNA is never recoded. Everything 3' of S is non-coding in this locus."""
    assert won.arm3 == locus.genome[won.site_s : won.span[1]]


def test_the_recoded_arm_keeps_the_protein_and_breaks_long_identity(won, locus):
    from Bio.Seq import Seq

    arm5_orig = locus.genome[won.span[0] : won.site_s]
    assert won.arm5 != arm5_orig or len(arm5_orig) <= 8
    frame = won.span[0] % 3  # the CDS starts at window position 0
    pad = (3 - frame) % 3
    a, b = arm5_orig[pad:], won.arm5[pad:]
    a, b = a[: len(a) // 3 * 3], b[: len(b) // 3 * 3]
    assert str(Seq(a).translate()) == str(Seq(b).translate())


def test_the_cargo_is_used_exactly_as_given(won):
    assert won.cargo == CARGO
    assert won.changes is not None and won.changes.cargo_unchanged
    assert not won.changes.cargo_changes


def test_both_pegrnas_are_within_the_calibrated_docking_range_and_at_most_166_nt(won):
    for g in (won.guide_a, won.guide_b):
        assert len(g.pegRNA_sequence) == 20 + 76 + g.flap_len + g.pbs_len + 3 <= 166
        assert 96 <= len(g.pegRNA_sequence) <= 232


def test_the_pegrnas_are_what_the_simulator_accepts(won, locus):
    """Re-run the self-check from the result alone, as a reader of the report could."""
    guides = {(g.strand, g.cut): g for g in find_guides(locus.genome)}
    a, b = guides[("top", won.guide_a.cut)], guides[("bottom", won.guide_b.cut)]
    out = simulate_pa(
        genome=locus.genome,
        a=a,
        b=b,
        donor=won.donor,
        peg_a=won.guide_a.pegRNA_sequence,
        peg_b=won.guide_b.pegRNA_sequence,
        flap_a=won.guide_a.flap_len,
        pbs_a=won.guide_a.pbs_len,
        flap_b=won.guide_b.flap_len,
        pbs_b=won.guide_b.pbs_len,
    )
    assert out == won.product and won.simulated


def test_the_flaps_are_the_ends_of_the_donor(won):
    assert rc(won.guide_a.rtt) == won.donor[: won.guide_a.flap_len]
    assert won.guide_b.rtt == won.donor[-won.guide_b.flap_len :]


def test_the_flap_lengths_are_in_the_paper_range_and_the_cargo_rule_holds(won):
    for g in (won.guide_a, won.guide_b):
        assert 30 <= g.flap_len <= 50
    assert 2 * max(won.guide_a.flap_len, won.guide_b.flap_len) <= len(won.cargo)


def test_the_standing_caveat_is_always_the_first_warning(won):
    assert won.warnings[0] == STANDING_CAVEAT
    assert "extension" in STANDING_CAVEAT and "untested" in STANDING_CAVEAT


def test_the_gate_results_name_every_step_in_pipeline_order(won):
    ids = [k for k, _ in won.gate_results]
    assert ids == [
        "offtarget",
        "recode",
        "inv0",
        "strand",
        "cryptic",
        "splice",
        "importance",
        "structure_a",
        "structure_b",
        "docking_a",
        "docking_b",
    ]
    info = dict(won.gate_results)
    assert "INV-0 held" in info["inv0"] and "simulate_pa" in info["strand"]


def test_counts_are_reported(won):
    assert won.pairs_considered == 3 and 1 <= won.candidates_tried <= won.pairs_considered


# ---- pairs and ranking ----------------------------------------------------------------------------------


def test_the_planning_finds_every_pam_in_pair_and_ranks_them_best_first(locus):
    planning = plan_pairs(**_plan_args(locus))
    assert planning.pairs_found == 3
    assert {(p.pair.a.cut, p.pair.b.cut) for p in planning.plans} == {
        (123, 164),
        (123, 190),
        (147, 190),
    }
    keys = [p.rank_key for p in planning.plans]
    assert keys == sorted(keys) and len(set(keys)) == len(keys), "a strict, total order"


def _plan_args(loc, **over):
    args = dict(
        cargo=CARGO,
        genome=loc.genome,
        site_s=loc.site_s,
        cds_mask=loc.cds_mask,
        offtarget_scorer=clean_offtarget,
        pbs_lengths=(13, 15, 17),
        flap_length=35,
        flap_min=30,
        flap_max=50,
        min_span=30,
        max_span=100,
        max_run=8,
        tail="UUU",
        max_offtargets=0,
        min_ontarget=0.5,
    )
    args.update(over)
    return args


def test_a_pair_that_breaks_a_pam_outranks_one_that_cannot():
    """Same locus and guides, only guide A's PAM differs. The arm is 6 nt, so no recoder runs and the PAM edit is
    the only change: CGG (Arg) can be broken to CGC or CGT; TGG (Trp) has no synonymous change at all.
    """

    def best(pam):
        loc = make_locus(
            [("top", 144, pam), ("bottom", 190)], n=300, site_s=S, seed=3, coding_end=S
        )
        return plan_pairs(**_plan_args(loc)).plans[0]

    breakable, regenerated = best("CGG"), best("TGG")
    assert breakable.pam_outcomes[0].status == "broken"
    assert regenerated.pam_outcomes[0].status == "intact-no-synonymous-change"
    assert len(regenerated.advisories) == len(breakable.advisories) + 1
    assert breakable.rank_key < regenerated.rank_key


def test_the_rank_key_orders_by_advisories_then_runs_then_span_then_on_target_then_cuts():
    base = dict(advisories=2, run_total=10, span=40, on_target_sum=1.6, a_cut=100, b_cut=140)

    def key(**over):
        return rank_key(**{**base, **over})

    assert key(advisories=1) < key()  # fewer advisories wins over everything else
    assert key(advisories=1, run_total=99, span=99) < key()
    assert key(run_total=9) < key() and key(run_total=9, span=99) < key()
    assert key(span=39) < key() and key(span=39, on_target_sum=0.0) < key()
    assert key(on_target_sum=1.9) < key()  # a HIGHER on-target sum is better
    assert key(a_cut=99) < key() and key(b_cut=139) < key()
    assert key(a_cut=99, b_cut=999) < key(
        a_cut=100, b_cut=1
    )  # a's cut before b's: the order is total


def test_the_span_is_the_number_of_restored_base_pairs_so_they_are_one_criterion(locus):
    for plan in plan_pairs(**_plan_args(locus)).plans:
        assert len(plan.arm5.sense) + len(plan.arm3.sense) == plan.pair.span


def test_a_shorter_span_ranks_above_a_higher_on_target_score(locus):
    """The span is a restored-base-pair count (fewer recoding edits), which outranks the on-target score."""
    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    favoured = guides[
        ("bottom", 190)
    ]  # every pair with it gets a better score, but it is the wider one

    def scorer(ps):
        return {p: (0.99 if p == favoured else 0.6, 0, False) for p in ps}

    plans = plan_pairs(**_plan_args(locus, offtarget_scorer=scorer)).plans
    assert [(p.pair.a.cut, p.pair.b.cut) for p in plans][0] == (123, 164)


def test_the_planning_screens_guides_in_one_batched_call_and_only_those_in_a_pair(locus):
    calls = []

    def spy(protospacers):
        calls.append(list(protospacers))
        return clean_offtarget(protospacers)

    plan_pairs(**_plan_args(locus, offtarget_scorer=spy))
    assert len(calls) == 1, "cas-offinder is ~100 s per invocation: one batched call"
    in_pairs = {
        g.spacer
        for g in find_guides(locus.genome)
        if (g.strand, g.cut) in {("top", 123), ("top", 147), ("bottom", 164), ("bottom", 190)}
    }
    assert set(calls[0]) <= in_pairs and len(set(calls[0])) == len(calls[0])


def test_a_guide_that_fails_the_off_target_screen_takes_its_pairs_with_it(locus):
    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    planning = plan_pairs(
        **_plan_args(locus, offtarget_scorer=_dirty_for({guides[("bottom", 164)]}))
    )
    assert {(p.pair.a.cut, p.pair.b.cut) for p in planning.plans} == {(123, 190), (147, 190)}
    assert [p.rank_key for p in planning.plans] == sorted(p.rank_key for p in planning.plans)
    assert any("failed the off-target screen" in d for d in planning.diagnostics)


def test_an_off_target_tolerance_lets_a_guide_through_and_is_counted(locus):
    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    dirty = _dirty_for({guides[("bottom", 164)]})
    strict = plan_pairs(**_plan_args(locus, offtarget_scorer=dirty))
    lax = plan_pairs(**_plan_args(locus, offtarget_scorer=dirty, max_offtargets=4))
    assert len(lax.plans) == len(strict.plans) + 1
    tolerated = next(p for p in lax.plans if p.pair.b.cut == 164)
    assert tolerated.offtargets == (0, 4)


def test_a_duplicate_of_a_guide_elsewhere_is_never_tolerated(locus):
    from pegrna_pipeline.hpc.offtarget_scorer import OffTargetScore

    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    target = guides[("bottom", 164)]

    def scorer(ps):
        return {
            p: OffTargetScore(0.8, 1, True, duplicates=1) if p == target else (0.8, 0, False)
            for p in ps
        }

    planning = plan_pairs(**_plan_args(locus, offtarget_scorer=scorer, max_offtargets=9))
    assert all(p.pair.b.cut != 164 for p in planning.plans)


def test_no_pam_in_pair_is_a_clear_failure_not_an_empty_result(locus):
    with pytest.raises(NoValidPADesign) as exc:
        run(locus, min_span=95, max_span=100)
    assert "0 PAM-in pairs" in str(exc.value) or "no PAM-in pair" in str(exc.value)


def test_every_guide_failing_the_screen_is_reported_with_the_reason(locus):
    with pytest.raises(NoValidPADesign) as exc:
        run(locus, offtarget_scorer=lambda ps: {p: (0.8, 9, True) for p in ps})
    text = " ".join(exc.value.diagnostics)
    assert "off-target" in text and "at most 0" in text


def test_a_span_floor_below_the_project_minimum_is_refused(locus):
    with pytest.raises(ValueError, match="at least 23"):
        plan_pairs(**_plan_args(locus, min_span=10))


# ---- what each gate rejects, and what it says ------------------------------------------------------------


def _fails(loc, needle, **over):
    with pytest.raises(NoValidPADesign) as exc:
        run(loc, **over)
    text = " ".join(exc.value.diagnostics)
    assert needle in text, text[:600]
    return exc.value


class DonorsInTheCargo:
    """Scores a strong donor only for the 9-mers the CARGO contains. The reference window has none of them, so
    the edit creates a strong site, which is what the gate rejects. (A scorer that scored every window strongly
    would not: the gate ignores a site no stronger than one the reference already had.)"""

    def __init__(self, cargo):
        self._nine = {cargo[i : i + 9] for i in range(len(cargo) - 8)}

    def score5(self, nine_mer):
        return 20.0 if nine_mer in self._nine else 0.0

    def score3(self, twentythree_mer):  # noqa: ARG002
        return 0.0


def test_the_cryptic_splice_gate_rejects_and_says_so(locus):
    _fails(locus, "cryptic splice", maxent_scorer=DonorsInTheCargo(CARGO))


def test_a_cryptic_site_no_stronger_than_one_the_reference_has_is_not_a_rejection(locus):
    assert run(locus, maxent_scorer=MaxEntScanFake(default=20.0)).donor


def test_the_binding_splice_gate_rejects_and_says_so(locus):
    _fails(
        locus, "splice (SpliceAI cross-check): splice delta 0.500", splice_scorer=SpliceFake(0.5)
    )


def test_the_structure_gate_rejects_naming_the_helix_and_the_pegrna(locus):
    err = _fails(
        locus, "cross-region pairing (CrossRegionFold)", full_structure_scorer=CrossRegionFold()
    )
    text = " ".join(err.diagnostics)
    assert "pegRNA A" in text or "A (" in text
    assert "RTT[" in text and "PBS1[" in text


def test_the_docking_gate_rejects_and_says_so(locus):
    _fails(locus, "docking", dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.10))


def test_the_importance_gate_rejects_once_for_the_whole_run_because_it_does_not_depend_on_the_pair(
    locus,
):
    err = _fails(locus, "insertion site", importance_scorer=FlagBlockedImportance())
    assert "no pair was tried" in " ".join(err.diagnostics)
    assert len([d for d in err.diagnostics if "insertion site" in d]) == 1


def test_a_scorer_that_crashes_is_a_diagnostic_not_a_crash(locus):
    class Boom:
        def pair_probs(self, seq):
            raise RuntimeError("folder exploded")

    _fails(locus, "structure scorer error: folder exploded", full_structure_scorer=Boom())

    class BoomDock(DockFake):
        def interface_iptm(self, pegrna, cas9):
            raise RuntimeError("dock exploded")

    _fails(locus, "docking scorer error: dock exploded", dock_scorer=BoomDock(0.85, 0.8))


def test_a_docking_outage_is_not_a_verdict_and_stops_the_run(locus):
    class Down(DockFake):
        def interface_iptm(self, pegrna, cas9):
            raise DockingServiceError("the service is down")

    with pytest.raises(DockingServiceError):
        run(locus, dock_scorer=Down(0.85, 0.8))


def test_a_failure_lists_the_pairs_that_were_tried_and_how_many_exist(locus):
    with pytest.raises(NoValidPADesign) as exc:
        run(locus, dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.10))
    assert "after 3 tried" in str(exc.value) and "3 PAM-in pair" in " ".join(exc.value.diagnostics)


def test_the_number_of_pairs_walked_is_bounded(locus):
    with pytest.raises(NoValidPADesign) as exc:
        run(locus, dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.10), max_pairs=2)
    assert "after 2 tried" in str(exc.value)


def test_the_walk_moves_on_to_the_next_pair_when_one_fails(locus):
    """Fail only the best pair's guide A: the walk must carry on and win with another pair."""
    best = plan_pairs(**_plan_args(locus)).plans[0]
    sacrificed = best.candidates_a[0].pegrna

    class OnePegRNAFails(DockFake):
        def interface_iptm(self, pegrna, cas9):
            return 0.1 if pegrna in {c.pegrna for c in best.candidates_a} else 0.82

    res = run(locus, dock_scorer=OnePegRNAFails(apo_iptm=0.85, candidate_iptm=0.82))
    assert (res.guide_a.cut, res.guide_b.cut) != (best.pair.a.cut, best.pair.b.cut)
    assert res.guide_a.pegRNA_sequence != sacrificed and res.candidates_tried >= 2


def test_the_pbs_length_is_a_lever_when_one_length_fails(locus):
    """A pegRNA that fails at PBS 13 is retried at 15 and 17 before the pair is given up on."""
    best = plan_pairs(**_plan_args(locus)).plans[0]
    first_a = best.candidates_a[0].pegrna

    class FirstLengthFails(DockFake):
        def interface_iptm(self, pegrna, cas9):
            return 0.1 if pegrna == first_a else 0.82

    res = run(locus, dock_scorer=FirstLengthFails(apo_iptm=0.85, candidate_iptm=0.82))
    assert (res.guide_a.cut, res.guide_b.cut) == (best.pair.a.cut, best.pair.b.cut)
    assert res.guide_a.pbs_len != best.candidates_a[0].extension.pbs_len


# ---- nothing expensive runs twice -------------------------------------------------------------------


def test_each_pegrna_is_docked_at_most_once_however_many_pairs_contain_it(locus):
    dock = RecordingDock(candidate_iptm=0.10)  # everything fails, so every pair is walked
    with pytest.raises(NoValidPADesign):
        run(locus, dock_scorer=dock)
    assert dock.shown, "the walk must reach docking for this to mean anything"
    assert len(dock.shown) == len(set(dock.shown)), "a pegRNA was docked more than once"


def test_a_pegrna_shared_by_two_pairs_is_folded_once_too(locus):
    fold = RecordingFold()
    with pytest.raises(NoValidPADesign):
        run(locus, full_structure_scorer=fold, dock_scorer=DockFake(0.85, 0.10))
    assert len(fold.seen) == len(set(fold.seen))


def test_the_importance_gate_runs_once_for_the_whole_walk(locus):
    calls = []

    class Counting(AllPassImportance):
        def metrics(self, i):
            calls.append(i)
            return super().metrics(i)

    with pytest.raises(NoValidPADesign):
        run(locus, importance_scorer=Counting(), dock_scorer=DockFake(0.85, 0.10))
    first_pass = len(calls)
    assert first_pass > 0
    calls.clear()
    with pytest.raises(NoValidPADesign):
        run(locus, importance_scorer=Counting(), dock_scorer=DockFake(0.85, 0.10), max_pairs=1)
    assert len(calls) == first_pass, "one pair or three, the insertion site is judged once"


# ---- warnings -------------------------------------------------------------------------------------------


def test_an_intact_pam_in_non_coding_dna_is_warned_about(won):
    assert any("PAM" in w and "non-coding" in w for w in won.warnings)


def test_a_flap_that_is_a_copy_of_the_genome_under_it_is_warned_about(won):
    assert any("re-anneal" in w and "identical" in w for w in won.warnings)


def test_an_advisory_splice_number_warns_and_never_rejects(locus):
    res = run(locus, splice_scorer=AdvisorySplice(0.34))
    assert any(w.startswith("splice cross-check:") and "does not reject" in w for w in res.warnings)
    assert "advisory only" in dict(res.gate_results)["splice"]


def test_a_binding_splice_verdict_names_where_it_was_computed(locus):
    class Bound(AdvisorySplice):
        binding = True
        context_note = "real genomic context (Ensembl test)"

    info = dict(run(locus, splice_scorer=Bound(0.05)).gate_results)["splice"]
    assert "real genomic context (Ensembl test)" in info and "advisory" not in info


def test_scaffold_pairing_is_an_advisory_not_a_rejection(locus):
    res = run(locus, full_structure_scorer=ScaffoldHuggingFold())
    assert any("scaffold nucleotides" in w and "pegRNA A" in w for w in res.warnings)
    assert any("pegRNA B" in w and "scaffold nucleotides" in w for w in res.warnings)


def test_a_homopolymer_in_a_spacer_is_reported_and_costs_rank(locus, monkeypatch):
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf

    real = wf.assess_spacer
    target = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}[("top", 123)]

    def patched(spacer):
        v = real(spacer)
        if spacer == target:
            return dataclasses.replace(
                v, advisories=("homopolymer run CCCCC (5 nt, positions 14-18)",)
            )
        return v

    monkeypatch.setattr(wf, "assess_spacer", patched)
    after = plan_pairs(**_plan_args(locus))
    monkeypatch.setattr(wf, "assess_spacer", real)
    before = plan_pairs(**_plan_args(locus))
    key_before = {(p.pair.a.cut, p.pair.b.cut): p.rank_key for p in before.plans}
    flagged = [p for p in after.plans if p.pair.a.cut == 123]
    assert flagged and all(any("CCCCC" in a for a in p.advisories) for p in flagged)
    for p in after.plans:
        cuts = (p.pair.a.cut, p.pair.b.cut)
        extra = 1 if p.pair.a.cut == 123 else 0
        assert (
            p.rank_key[0] == key_before[cuts][0] + extra
        ), "exactly the flagged guide's pairs lose one place"


def test_a_tolerated_off_target_is_disclosed_on_the_design(locus):
    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    dirty = _dirty_for({guides[("bottom", 164)], guides[("bottom", 190)]})
    res = run(locus, offtarget_scorer=dirty, max_offtargets=4)
    assert any("off-target site(s)" in w and "4" in w for w in res.warnings)
    assert "4 off-target" in dict(res.gate_results)["offtarget"]


# ---- the two things that restoring the span makes necessary ------------------------------------------------


@pytest.fixture(scope="module")
def short_arm_locus():
    """Guide A cuts 6 nt 5' of S, so its arm is too short to need any recoding (under either recoder): its CGG PAM
    is intact when the PAM is broken, which makes the dedicated PAM edit deterministic."""
    return make_locus([("top", 144, "CGG"), ("bottom", 190)], n=300, site_s=S, seed=3, coding_end=S)


@pytest.fixture(scope="module")
def weak_arm_locus():
    """The same locus with an AGG PAM: arginine, but the only synonymous change is AGG -> AGA, which is NGA."""
    return make_locus([("top", 144, "AGG"), ("bottom", 190)], n=300, site_s=S, seed=3, coding_end=S)


def test_a_pam_in_coding_dna_is_broken_by_one_synonymous_edit_and_reported(short_arm_locus):
    res = run(short_arm_locus)
    assert (res.guide_a.cut, res.guide_a.pam) == (144, "CGG")
    assert res.guide_a.pam_status == "broken"
    pam_edits = [c for c in res.arm_changes if "PAM" in c.part]
    assert len(pam_edits) == 1
    edit = pam_edits[0]
    assert (edit.old, edit.aa) == (
        "G",
        "R",
    ) and edit.new in "CT", "CGC and CGT are Arg and leave no PAM"
    assert res.guide_a.product_pam in {"CGC", "CGT"}
    assert "donor arm 5'" in edit.part and "breaks the PAM of guide A" in edit.part
    assert res.recode_edits == 1, "the PAM edit is the only change the short arm needs"
    assert not any("guide A's PAM" in w for w in res.warnings), "a broken PAM is not a warning"
    info = dict(res.gate_results)["recode"]
    assert "guide A's PAM broken at nt 150 by one synonymous change" in info


def test_a_pam_that_can_only_be_weakened_is_reported_as_that_and_warned_about(weak_arm_locus):
    """NAG and NGA are still cut by SpCas9 (about a fifth of NGG for NAG), so AGG -> AGA is not a break."""
    res = run(weak_arm_locus)
    assert (res.guide_a.cut, res.guide_a.pam) == (144, "AGG")
    assert res.guide_a.pam_status == "weakened" and res.guide_a.product_pam == "AGA"
    edit = next(c for c in res.arm_changes if "PAM" in c.part)
    assert (edit.old, edit.new, edit.aa) == ("G", "A", "R")
    assert "weakens the PAM of guide A" in edit.part and "breaks" not in edit.part
    warning = next(w for w in res.warnings if "guide A's PAM" in w)
    assert "can only be weakened, not broken" in warning and "reads AGA" in warning
    assert "weakened to AGA at nt 150" in dict(res.gate_results)["recode"]


def test_a_weakened_pam_ranks_below_a_pair_whose_pam_is_really_broken():
    """Same locus, same guides; only guide A's PAM differs (CGG can be broken, AGG only weakened)."""

    def best_key(pam):
        loc = make_locus(
            [("top", 144, pam), ("bottom", 190)], n=300, site_s=S, seed=3, coding_end=S
        )
        return plan_pairs(**_plan_args(loc)).plans[0]

    broken, weak = best_key("CGG"), best_key("AGG")
    assert broken.pam_outcomes[0].status == "broken" and weak.pam_outcomes[0].status == "weakened"
    assert len(weak.advisories) == len(broken.advisories) + 1
    assert weak.rank_key > broken.rank_key
    assert any("can only be weakened" in a for a in weak.advisories)


@pytest.mark.parametrize(
    "outcome, words",
    [
        (("broken", 150, "CGC"), "guide A's PAM broken at nt 151 by one synonymous change"),
        (("weakened", 150, "AAG"), "guide A's PAM weakened to AAG at nt 151"),
        (("weakened", None, "AAG"), "guide A's PAM left weak (AAG) by the recode"),
        (("broken-by-recode", None, "TTG"), "guide A's PAM broken by the recode"),
        (("intact-noncoding", None, "TGG"), "guide A's PAM regenerated"),
        (("intact-no-synonymous-change", None, "TGG"), "guide A's PAM regenerated"),
    ],
)
def test_the_recode_summary_says_what_became_of_each_pam(locus, outcome, words):
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf
    from pegrna_pipeline.hpc.pa_geometry import PamOutcome

    plan = plan_pairs(**_plan_args(locus)).plans[0]
    status, position, pam = outcome
    other = PamOutcome("B", "intact-noncoding", None, "AGG")
    shown = dataclasses.replace(plan, pam_outcomes=(PamOutcome("A", status, position, pam), other))
    text = wf._recode_info(shown)
    assert words in text and "guide B's PAM regenerated" in text
    assert f"(5' arm {plan.arm5.n_edits}, 3' arm {plan.arm3.n_edits})" in text


def test_the_pam_edit_is_in_the_donor_and_not_in_the_cargo_or_the_genome_window(short_arm_locus):
    res = run(short_arm_locus)
    assert res.arm5[-1] in "CT" and short_arm_locus.genome[res.span[0] : res.site_s][-1] == "G"
    assert res.cargo == CARGO and res.changes.cargo_unchanged


def test_the_two_donor_ends_must_not_pair_with_each_other(locus):
    """A 12-nt stretch shared by the donor's two ends lets the flaps anneal to each other, not to the donor."""
    shared = locus.genome[170:182]
    cargo = shared + make_cargo(108)
    planning = plan_pairs(**_plan_args(locus, cargo=cargo))
    assert all((p.pair.a.cut, p.pair.b.cut) != (147, 190) for p in planning.plans)
    assert any("pair with each other" in d and "147" in d for d in planning.diagnostics)


def test_a_terminator_in_every_flap_rejects_the_pair_and_the_others_carry_on(locus):
    cargo = (
        "ACGTACGTAC" + "AAAA" + make_cargo(106)
    )  # AAAA at 10: RC is TTTT in a flap that reaches it
    planning = plan_pairs(**_plan_args(locus, cargo=cargo))
    cuts = {(p.pair.a.cut, p.pair.b.cut) for p in planning.plans}
    assert (123, 164) in cuts and not any(c[0] == 147 for c in cuts)
    assert any("TTTT" in d and "147" in d for d in planning.diagnostics)


def test_a_guide_cutting_exactly_at_s_has_an_empty_arm_and_the_donor_starts_with_the_cargo():
    loc = make_locus([("top", S), ("bottom", 190)], n=300, site_s=S, seed=9, coding_end=S)
    planning = plan_pairs(**_plan_args(loc))
    assert len(planning.plans) == 1
    plan = planning.plans[0]
    assert plan.arm5.sense == "" and plan.donor.startswith(CARGO)
    res = run(loc)
    assert res.arm5 == "" and res.donor.startswith(CARGO) and res.deletion_size == 0


def test_an_arm_the_recoder_cannot_solve_is_a_diagnostic_not_a_crash(locus, monkeypatch):
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf

    real = wf.build_arm

    def failing(genome, lo, hi, mask, max_run):
        return None if lo == 123 else real(genome, lo, hi, mask, max_run)

    monkeypatch.setattr(wf, "build_arm", failing)
    planning = plan_pairs(**_plan_args(locus))
    assert all(p.pair.a.cut != 123 for p in planning.plans)
    assert any("arm partition failed" in d and "123" in d for d in planning.diagnostics)


def test_a_strand_mismatch_found_by_the_simulator_rejects_the_candidate(locus, monkeypatch):
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf

    def broken(**kw):
        raise PAMismatch("flap A does not match the 5' end of the donor")

    monkeypatch.setattr(wf, "simulate_side", broken)
    planning = plan_pairs(**_plan_args(locus))
    assert not planning.plans
    assert any("PA strand mismatch" in d for d in planning.diagnostics)


def test_a_coding_3_prime_arm_is_recoded_and_a_bottom_guides_pam_can_be_broken(coding_locus):
    planning = plan_pairs(**_plan_args(coding_locus))
    assert planning.plans, planning.diagnostics
    plan = planning.plans[0]
    assert plan.arm3.recode_span is not None and plan.arm3.n_edits > 0
    # it is in coding sequence now, so it is never classed as non-coding: it is broken or has no synonymous option
    assert plan.pam_outcomes[1].status != "intact-noncoding"


# ---- progress events -------------------------------------------------------------------------------------


def test_progress_events_carry_the_gate_the_verdict_and_the_numbers(locus):
    events = []
    run(
        locus,
        on_gate=lambda cand, gate, info, passed, metrics=None: events.append(
            (cand, gate, passed, info, metrics or {})
        ),
    )
    gates = [e[1] for e in events]
    for gate in (
        "cryptic",
        "splice",
        "importance",
        "structure_a",
        "structure_b",
        "docking_a",
        "docking_b",
    ):
        assert gate in gates, gate
    assert all(isinstance(e[0], int) and e[0] >= 1 for e in events)
    assert all(e[2] is True for e in events)
    cryptic = next(e for e in events if e[1] == "cryptic")
    assert cryptic[4].get("metric") == "MaxEnt"


def test_an_old_four_argument_callback_still_works(locus):
    seen = []
    run(locus, on_gate=lambda cand, gate, info, passed: seen.append(gate))
    assert "docking_a" in seen


def test_a_callback_that_raises_never_breaks_the_run(locus):
    def boom(*a, **k):
        raise RuntimeError("the UI went away")

    assert run(locus, on_gate=boom).donor


def test_a_rejection_is_streamed_as_a_failed_gate(locus):
    events = []
    with pytest.raises(NoValidPADesign):
        run(
            locus,
            dock_scorer=DockFake(0.85, 0.10),
            on_gate=lambda cand, gate, info, passed, metrics=None: events.append((gate, passed)),
        )
    assert ("docking_a", False) in events


# ---- the lock model: what may change, and what is reported -----------------------------------------------


def test_every_nucleotide_that_differs_from_the_window_is_listed(won, locus):
    n1, n2 = won.span
    listed = {c.position - 1 for c in won.arm_changes}
    actual = {
        n1 + i for i, (x, y) in enumerate(zip(won.arm5 + won.arm3, locus.genome[n1:n2])) if x != y
    }
    # arm5 + arm3 are the donor's two arms around the cargo; the cargo itself is never listed
    assert len(won.arm_changes) == won.recode_edits
    assert listed == actual
    assert all(c.old != c.new for c in won.arm_changes)


def test_the_report_object_says_what_was_kept_exactly(won):
    kept = " | ".join(won.changes.arm_kept_exactly)
    assert "non-coding" in kept and "arm" in kept
    assert not any(
        "arm" in line for line in won.changes.kept_exactly
    ), "an arm is not part of the cargo"


def test_user_locks_are_validated_and_reported_though_the_cargo_is_never_recoded(locus):
    res = run(locus, cargo_locked=[[1, 12]])
    assert any(
        line.startswith("nt 1-12 (12 nt): excluded by you") for line in res.changes.kept_exactly
    )
    parts = [
        {"name": "ORF", "dna": CARGO[:90], "locked": False},
        {"name": "tail", "dna": CARGO[90:], "locked": True},
    ]
    res = run(locus, cargo_parts=parts)
    assert any("tail" in line for line in res.changes.kept_exactly) and res.changes.cargo_unchanged
    with pytest.raises(ValueError, match="cargo_locked"):
        run(locus, cargo_locked=[[1, 999]])
    with pytest.raises(ValueError, match="do not add up"):
        run(locus, cargo_parts=[{"name": "x", "dna": CARGO[:-3], "locked": False}])


def test_a_t7_promoter_in_the_cargo_is_listed_as_kept_exactly_and_never_changed(locus):
    from pegrna_pipeline.hpc.insert_elements import T7G_PROMOTER

    cargo = make_cargo(111) + "TGA" + T7G_PROMOTER
    res = run(locus, cargo=cargo)
    assert res.cargo == cargo and res.changes.cargo_unchanged
    assert any("T7 promoter" in line for line in res.changes.kept_exactly)


# ---- wiring and invariants that nothing else would notice if they were removed ------------------------


def _corrupting(monkeypatch, lo_target, mutate):
    """Make ``build_arm`` hand back a corrupted arm for guide A's cut, as a recoder bug would."""
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf

    real = wf.build_arm

    def build(genome, lo, hi, mask, max_run):
        arm = real(genome, lo, hi, mask, max_run)
        if arm is not None and lo == lo_target and arm.sense:
            return dataclasses.replace(arm, sense=mutate(arm.sense))
        return arm

    monkeypatch.setattr(wf, "build_arm", build)


def _missense(sense):
    """Change the first codon's second base to one that changes the amino acid."""
    from Bio.Seq import Seq

    first = str(Seq(sense[:3]).translate())
    for base in "ACGT":
        if base != sense[1] and str(Seq(sense[0] + base + sense[2]).translate()) != first:
            return sense[0] + base + sense[2:]
    raise AssertionError("no missense found")


def test_an_arm_that_changes_the_protein_fails_inv0_and_the_pair_is_dropped(locus, monkeypatch):
    _corrupting(monkeypatch, 123, _missense)
    planning = plan_pairs(**_plan_args(locus))
    assert all(p.pair.a.cut != 123 for p in planning.plans) and planning.plans
    inv0 = [d for d in planning.diagnostics if "INV-0" in d]
    assert inv0 and all("123" in d and "not protein-silent" in d for d in inv0)
    assert not any("INV-0: INV-0" in d for d in inv0), "the message must not say INV-0 twice"


def test_an_arm_that_changes_the_length_is_caught_before_it_becomes_a_product(locus, monkeypatch):
    _corrupting(monkeypatch, 123, lambda sense: sense + "A")
    planning = plan_pairs(**_plan_args(locus))
    assert all(p.pair.a.cut != 123 for p in planning.plans)
    assert any("the product is" in d and "123" in d for d in planning.diagnostics)


def test_inv0_is_what_the_planning_calls_on_each_arm(locus, monkeypatch):
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf

    calls = []
    real = wf.assert_inv0
    monkeypatch.setattr(wf, "assert_inv0", lambda **kw: (calls.append(kw), real(**kw))[1])
    plan_pairs(**_plan_args(locus))
    assert len(calls) == 2 * 3, "two arms for each of the three pairs"
    assert all(
        c["allowed_silent"] is None for c in calls
    ), "no silent change is allowed outside the arms"


def test_the_strand_simulator_is_run_on_the_chosen_pair_and_a_mismatch_rejects_it(
    locus, monkeypatch
):
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf

    def broken(**kw):
        raise PAMismatch("flap B does not match the 3' end of the donor")

    monkeypatch.setattr(wf, "simulate_pa", broken)
    with pytest.raises(NoValidPADesign) as exc:
        run(locus)
    assert "PA strand mismatch" in " ".join(exc.value.diagnostics)


def test_the_simulator_is_given_the_expected_product(locus, monkeypatch):
    from pegrna_pipeline.hpc import recoded_pa_workflow as wf

    seen = []
    real = wf.simulate_pa
    monkeypatch.setattr(wf, "simulate_pa", lambda **kw: (seen.append(kw), real(**kw))[1])
    res = run(locus)
    assert len(seen) == 1 and seen[0]["expected_product"] == res.product


def test_every_candidate_of_a_side_has_the_same_flap_length(locus, coding_locus):
    """The reason the donor-ends check can be made once, when the pair is planned."""
    for loc in (locus, coding_locus):
        for plan in plan_pairs(**_plan_args(loc)).plans:
            assert len({c.extension.flap_len for c in plan.candidates_a}) == 1
            assert len({c.extension.flap_len for c in plan.candidates_b}) == 1
            assert len({c.extension.pbs_len for c in plan.candidates_a}) == len(plan.candidates_a)


def test_the_flap_never_exceeds_half_the_cargo(locus):
    """With a 70-nt cargo the longest flap is 35. Make flap 35 unattractive (its RTT would start with C) for
    both top guides, so a scan allowed past 35 would take 36 instead of stepping back to 34."""
    cargo = list(CARGO[:70])
    cargo[7], cargo[8] = "G", "A"  # guide 123: donor[34] = cargo[7]; donor[35] = cargo[8]
    cargo[30], cargo[31], cargo[32] = "A", "G", "A"  # guide 147: donor[33..35] = cargo[30..32]
    plans = plan_pairs(**_plan_args(locus, cargo="".join(cargo))).plans
    assert plans
    for plan in plans:
        for cand in plan.candidates_a + plan.candidates_b:
            assert cand.extension.flap_len <= 35
    # flap 35 starts with C for guide 123, so the scan steps BACK from it (which flap it finds depends on the cargo)
    assert {
        c.extension.flap_len for p in plans if p.pair.a.cut == 123 for c in p.candidates_a
    } < set(range(30, 35))


def test_a_failure_on_one_side_spares_the_expensive_gates_for_that_pair(locus):
    """If guide A's pegRNAs cannot fold, the pair is doomed: do not spend SpliceAI or docking on it."""
    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    doomed = guides[("top", 123)]

    class FoldsBadlyForGuide123:
        def pair_probs(self, seq):
            return {(100, len(seq) - 6): 0.999} if seq[:20] == doomed else {}

    class CountingSplice(SpliceFake):
        def __init__(self):
            super().__init__(0.0)
            self.calls = 0

        def delta(self, ref_seq, alt_seq):
            self.calls += 1
            return 0.0

    splice = CountingSplice()
    dock = RecordingDock()
    res = run(
        locus, full_structure_scorer=FoldsBadlyForGuide123(), splice_scorer=splice, dock_scorer=dock
    )
    assert res.guide_a.cut == 147, "the only pair without guide 123 that can win"
    assert (
        splice.calls == 1
    ), "SpliceAI ran for the one pair that survived folding, not for the two doomed ones"
    assert all(not p.startswith(doomed) for p in dock.shown), "the doomed pairs were never docked"


def test_a_failure_on_guide_a_does_not_dock_guide_b(locus):
    dock = RecordingDock(candidate_iptm=0.10)
    with pytest.raises(NoValidPADesign):
        run(locus, dock_scorer=dock)
    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    b_spacers = {guides[("bottom", 164)], guides[("bottom", 190)]}
    assert all(p[:20] not in b_spacers for p in dock.shown), "A fails first, so B is never docked"


def test_recoding_edits_cover_both_arms_and_every_one_is_listed(coding_locus):
    """The 3' arm is coding here, so it is recoded too: the count and the listing must include it."""
    plan = plan_pairs(**_plan_args(coding_locus)).plans[0]
    assert plan.arm3.n_edits > 0
    res = run(coding_locus)
    assert res.recode_edits == res.arm_changes.__len__() == plan.arm5.n_edits + plan.arm3.n_edits
    assert {c.part.split(",")[0] for c in res.arm_changes} == {"donor arm 3'"} | (
        {"donor arm 5'"} if plan.arm5.n_edits else set()
    )


def test_the_results_guide_records_are_the_walks_own(won):
    for g, side in ((won.guide_a, "A"), (won.guide_b, "B")):
        assert (
            g.side == side
            and g.tail == "UUU"
            and g.pbs_len == len(g.pbs)
            and g.flap_len == len(g.rtt)
        )
        assert g.pegRNA_sequence.endswith("UUU") and g.spacer == g.pegRNA_sequence[:20]
        assert 0 <= g.genomic_run <= g.flap_len


def test_both_guides_off_target_counts_are_recorded_separately(locus):
    guides = {(g.strand, g.cut): g.spacer for g in find_guides(locus.genome)}
    scorer = _dirty_for({guides[("bottom", 164)]})
    res = run(locus, offtarget_scorer=scorer, max_offtargets=4)
    assert (res.guide_a.offtargets, res.guide_b.offtargets) == (0, 4)


# ---- determinism ---------------------------------------------------------------------------------------


def _digest(res):
    return json.dumps(
        [
            res.guide_a.pegRNA_sequence,
            res.guide_b.pegRNA_sequence,
            res.donor,
            res.product,
            res.warnings,
            res.gate_results,
            [(c.position, c.old, c.new) for c in res.arm_changes],
            res.candidates_tried,
            res.pairs_considered,
        ]
    )


def test_two_runs_in_one_process_agree(locus):
    assert _digest(run(locus)) == _digest(run(locus))


_PROBE = """
import json, sys
sys.path.insert(0, sys.argv[1])
from pa_helpers import fake_scorers, make_cargo, make_locus
from pegrna_pipeline.hpc.recoded_pa_workflow import run_recoded_pa
loc = make_locus([("top", 123, "AGG"), ("top", 147), ("bottom", 164), ("bottom", 190)], n=300, site_s=150, seed=7, coding_end=150)
res = run_recoded_pa(cargo=make_cargo(120), genome=loc.genome, site_s=loc.site_s, cds_mask=loc.cds_mask, **fake_scorers())
print(json.dumps([res.guide_a.pegRNA_sequence, res.guide_b.pegRNA_sequence, res.donor, res.product,
                  res.warnings, res.gate_results, [(c.position, c.old, c.new) for c in res.arm_changes],
                  res.candidates_tried, res.pairs_considered]))
"""


@pytest.mark.parametrize("hash_seeds", [(0, 1, 2, 5, 9)])
def test_the_design_is_identical_in_every_process(hash_seeds):
    """String hashing is randomised per process; a design must not depend on it (this is what #30 was)."""
    outputs = set()
    for seed in hash_seeds:
        proc = subprocess.run(
            [sys.executable, "-c", _PROBE, str(REPO / "tests")],
            capture_output=True,
            text=True,
            cwd=str(REPO),
            env={
                "PYTHONHASHSEED": str(seed),
                "PATH": "/usr/bin:/bin",
                "PYTHONPATH": str(REPO / "src"),
            },
            timeout=300,
        )
        assert proc.returncode == 0, proc.stderr[-800:]
        outputs.add(proc.stdout.strip().splitlines()[-1])
    assert len(outputs) == 1, "the design changed with PYTHONHASHSEED"


def test_the_inputs_are_never_modified(locus):
    genome, mask, cargo = locus.genome, list(locus.cds_mask), CARGO
    run(locus)
    assert (locus.genome, list(locus.cds_mask), CARGO) == (genome, mask, cargo)
