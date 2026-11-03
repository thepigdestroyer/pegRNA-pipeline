"""The second nick's off-target tolerance is explicit, configurable, and disclosed.

ADR 0010 chose zero tolerance for the ngRNA, and that stays the default. But the knob that was
supposed to relax it, ``max_offtargets``, did nothing above 0: the screens rejected on
``has_off or off > max_offtargets``, and the scorers set ``has_off = (off > 0)``, so any guide
with a single off-target was rejected whatever the limit. These tests pin the repaired rule on
both candidate paths, through the driver, and in what the report says about the chosen guide.
"""

import dataclasses
import json
from pathlib import Path

import pytest

# Reuse the workflow fixtures of the PRIDICT wiring tests (same fake genome / scorers).
import test_pridict_ngrna_wiring as W  # noqa: E402  (tests/ is on sys.path under pytest)

from pegrna_pipeline.hpc import recoded_driver
from pegrna_pipeline.hpc import recoded_tjpe_workflow as wf
from pegrna_pipeline.hpc.nick_search import (
    NickCandidate,
    offtarget_rejection,
    screen_nick_candidates,
)
from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer, OffTargetScore
from pegrna_pipeline.hpc.pridict_ngrna import build_pridict_ngrnas
from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA, run_recoded_tjpe
from pegrna_pipeline.hpc.scorer_factory import make_scorers

CLEAN = (1.0, 0, False)


def near(n: int):
    """A scorer reporting ``n`` near-matches for every guide (``n=0`` is clean)."""
    return lambda protospacers: {p: (1.0, n, n > 0) for p in protospacers}


# ---- the rule itself ----------------------------------------------------------------


def test_the_count_is_compared_with_the_limit():
    assert offtarget_rejection((1.0, 0, False), max_offtargets=0, min_ontarget=0.5) is None
    assert offtarget_rejection((1.0, 1, True), max_offtargets=0, min_ontarget=0.5)
    assert offtarget_rejection((1.0, 1, True), max_offtargets=1, min_ontarget=0.5) is None
    assert offtarget_rejection((1.0, 2, True), max_offtargets=1, min_ontarget=0.5)


def test_the_flag_alone_no_longer_overrides_the_limit():
    """The bug: ``has_off`` is true for any count > 0, so it defeated every limit above 0."""
    assert offtarget_rejection((1.0, 3, True), max_offtargets=5, min_ontarget=0.5) is None


def test_a_flag_with_no_count_is_contradictory_and_rejected():
    assert offtarget_rejection((1.0, 0, True), max_offtargets=3, min_ontarget=0.5)


def test_a_missing_on_target_is_rejected_whatever_the_limit():
    reason = offtarget_rejection((0.0, 0, False), max_offtargets=9, min_ontarget=0.5)
    assert reason and "on-target" in reason


def test_the_reason_names_the_count_and_the_limit():
    reason = offtarget_rejection((1.0, 4, True), max_offtargets=1, min_ontarget=0.5)
    assert "4" in reason and "1" in reason


# ---- exact duplicates are never a "tolerated" off-target ----------------------------------


class _StubCasOffinder(CasOffinderScorer):
    """CasOffinderScorer whose binary is replaced by canned cas-offinder rows."""

    def __init__(self, tiers, **kw):
        super().__init__(**kw)
        self._tiers = tiers

    def _run(self, input_path, output_path):
        with open(input_path) as fh:
            guide = fh.read().split()[2]  # first query line: <guide+NNN> <mismatches>
        guide = guide[:-3]
        with open(output_path, "w") as fh:
            for mm, n in self._tiers.items():
                for _ in range(n):
                    fh.write(
                        "\t".join(
                            str(x) for x in (guide + "NGG", "chr1", 1, guide + "TGG", "+", mm)
                        )
                        + "\n"
                    )


GUIDE = "ACGTACGTACGTACGTACGT"


def test_the_score_is_still_a_plain_triple():
    s = OffTargetScore(1.0, 3, True, duplicates=1)
    assert s == (1.0, 3, True) and tuple(s) == (1.0, 3, True)
    assert (lambda a, b, c: (a, b, c))(*s) == (1.0, 3, True)
    assert s.duplicates == 1 and getattr((1.0, 3, True), "duplicates", 0) == 0


def test_the_scorer_reports_how_many_hits_are_exact_copies():
    sc = _StubCasOffinder({0: 3, 1: 2}, genome="/x.fa")
    score = sc([GUIDE])[GUIDE]
    assert score == (1.0, 4, True)  # 2 extra exact copies + 2 near matches
    assert score.duplicates == 2


def test_near_matches_alone_report_no_duplicates():
    score = _StubCasOffinder({0: 1, 2: 3}, genome="/x.fa")([GUIDE])[GUIDE]
    assert score == (1.0, 3, True) and score.duplicates == 0


def test_a_duplicate_is_rejected_at_any_tolerance():
    s = OffTargetScore(1.0, 1, True, duplicates=1)
    reason = offtarget_rejection(s, max_offtargets=10, min_ontarget=0.5)
    assert reason and "exact copy" in reason and "never tolerated" in reason


def test_a_tolerance_admits_near_matches_but_not_a_duplicate_hiding_among_them():
    near_only = OffTargetScore(1.0, 2, True, duplicates=0)
    with_dup = OffTargetScore(1.0, 2, True, duplicates=1)
    assert offtarget_rejection(near_only, max_offtargets=2, min_ontarget=0.5) is None
    assert offtarget_rejection(with_dup, max_offtargets=2, min_ontarget=0.5)


def test_the_real_scorer_through_the_screen_never_lets_a_duplicate_through():
    cand = NickCandidate(GUIDE, "TGG", 205, 55, "Fw")
    sc = _StubCasOffinder({0: 2}, genome="/x.fa")  # one extra perfect copy
    assert screen_nick_candidates([cand], scorer=sc, max_offtargets=5) == []
    near = _StubCasOffinder({0: 1, 1: 1, 2: 1}, genome="/x.fa")
    (kept,) = screen_nick_candidates([cand], scorer=near, max_offtargets=2)
    assert kept.offtargets == 2


def test_the_pridict_path_has_the_same_protection():
    sc = lambda ps: {p: OffTargetScore(1.0, 1, True, duplicates=1) for p in ps}  # noqa: E731
    assert _build(sc, max_offtargets=5) == []


# ---- genome-enumerated path ---------------------------------------------------------

CAND = NickCandidate("TGGCATGACTGAACGATGCA", "TGG", 205, 55, "Fw")


def test_screen_default_is_still_zero_tolerance():
    assert screen_nick_candidates([CAND], scorer=near(1)) == []
    assert len(screen_nick_candidates([CAND], scorer=near(0))) == 1


def test_screen_honours_a_tolerance_and_reports_the_count():
    kept = screen_nick_candidates([CAND], scorer=near(2), max_offtargets=2)
    assert [c.offtargets for c in kept] == [2]
    assert screen_nick_candidates([CAND], scorer=near(3), max_offtargets=2) == []


def test_an_unscreened_candidate_has_no_count():
    assert CAND.offtargets is None


# ---- PRIDICT path -------------------------------------------------------------------


def _build(scorer, **kw):
    rows = [W._nrow(W._MARKER_A, "50", "Fw", "0.9")]
    return build_pridict_ngrnas(
        rows,
        genome=W.GENOME,
        site_s=W.SITE_S,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        offtarget_scorer=scorer,
        **kw,
    )


def test_pridict_path_default_rejects_a_single_near_match():
    assert _build(near(1)) == []


def test_pridict_path_honours_the_tolerance_and_carries_the_count():
    (inj,) = _build(near(1), max_offtargets=1)
    assert inj.offtargets == 1


def test_a_clean_pridict_guide_carries_zero():
    (inj,) = _build(near(0))
    assert inj.offtargets == 0


# ---- the walk: default, tolerance, disclosure ----------------------------------------


def test_default_workflow_rejects_guides_with_near_matches_and_says_why():
    with pytest.raises(NoValidPegRNA) as exc:
        run_recoded_tjpe(**W._wf_inputs(offtarget_scorer=near(1)))
    text = " ".join(exc.value.diagnostics)
    assert "failed the off-target screen" in text and "at most 0" in text


def test_a_granted_tolerance_lets_the_design_through_and_discloses_it():
    res = run_recoded_tjpe(**W._wf_inputs(offtarget_scorer=near(1)), ngrna_max_offtargets=1)
    off = dict(res.gate_results)["offtarget"]
    assert "1 off-target site(s)" in off and "tolerance 1" in off
    assert any("second-nick guide has 1 off-target" in w for w in res.warnings)


def test_a_clean_guide_reports_zero_and_raises_no_off_target_warning():
    res = run_recoded_tjpe(**W._wf_inputs(offtarget_scorer=near(0)), ngrna_max_offtargets=3)
    assert "0 off-target sites" in dict(res.gate_results)["offtarget"]
    assert not any("off-target" in w for w in res.warnings)


def test_the_tolerance_is_a_ceiling_not_a_target():
    """Guides with 2 near-matches still fail when only 1 is allowed."""
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(**W._wf_inputs(offtarget_scorer=near(2)), ngrna_max_offtargets=1)


@pytest.mark.parametrize("bad", [-1, 1.5])
def test_a_nonsensical_tolerance_is_refused(bad):
    with pytest.raises(ValueError, match="ngrna_max_offtargets"):
        run_recoded_tjpe(**W._wf_inputs(), ngrna_max_offtargets=bad)


def test_unusable_sequences_are_never_sent_to_the_genome_scan(monkeypatch):
    """The expensive step is a whole-genome scan; a TTTT guide cannot be used, so it must not
    reach it (and, before the shared check existed, it would have been chosen)."""
    ok = NickCandidate("TGGCATGACTGAACGATGCA", "TGG", 205, 55, "Fw")
    poly_t = NickCandidate("TTCTCCCTTCTTTCTTTTTT", "TGG", 195, 45, "Fw")
    monkeypatch.setattr(wf, "enumerate_nick_sites", lambda *a, **k: [poly_t, ok])
    scanned = []

    def spy(protospacers):
        scanned.extend(protospacers)
        return {p: CLEAN for p in protospacers}

    res = run_recoded_tjpe(**W._wf_inputs(offtarget_scorer=spy))
    assert scanned == [ok.protospacer]
    assert res.nick_position == 205


def test_when_every_candidate_is_unusable_the_reason_is_the_sequence(monkeypatch):
    poly_t = NickCandidate("TTCTCCCTTCTTTCTTTTTT", "TGG", 205, 55, "Fw")
    monkeypatch.setattr(wf, "enumerate_nick_sites", lambda *a, **k: [poly_t])
    with pytest.raises(NoValidPegRNA) as exc:
        run_recoded_tjpe(**W._wf_inputs())
    text = " ".join(exc.value.diagnostics)
    assert "failed the sequence rules" in text and "TTTT" in text


def test_homopolymer_advisories_surface_on_the_chosen_guide(monkeypatch):
    runs = NickCandidate("GCTGCCACAGATGCCCCCAA", "TGG", 205, 55, "Fw")  # CCCCC
    monkeypatch.setattr(wf, "enumerate_nick_sites", lambda *a, **k: [runs])
    res = run_recoded_tjpe(**W._wf_inputs(offtarget_scorer=near(0)))
    assert res.ngrna_protospacer == runs.protospacer, "a homopolymer must not disqualify it"
    assert any("ngRNA spacer advisory" in w and "CCCCC" in w for w in res.warnings)


def test_a_pegrna_spacer_that_breaks_a_hard_rule_is_reported_not_silently_used():
    res = run_recoded_tjpe(**W._wf_inputs(spacer="ACGTTTTACGCAGGCATCGA"))
    assert any(w.startswith("pegRNA spacer:") and "TTTT" in w for w in res.warnings)


# ---- through the driver (the path the CLI and the web app take) -----------------------


def _with_offtarget(monkeypatch, scorer):
    real = recoded_driver.make_scorers

    def patched(**kw):
        return dataclasses.replace(real(**kw), offtarget_scorer=scorer)

    monkeypatch.setattr(recoded_driver, "make_scorers", patched)


def _runner():
    return W._NickRunner(
        spacer="TGCATGCATGCATGCATGCA",
        nicking=[W._nrow(W._MARKER_B, "50", "Fw", "0.7")],
    )


def test_driver_default_is_zero_tolerance(monkeypatch):
    _with_offtarget(monkeypatch, near(1))
    with pytest.raises(NoValidPegRNA):
        recoded_driver.design_recoded(**W._auto_inputs(), mode="fake", pridict_runner=_runner())


def test_driver_passes_the_tolerance_to_both_screens(monkeypatch):
    _with_offtarget(monkeypatch, near(1))
    out = recoded_driver.design_recoded(
        **W._auto_inputs(), mode="fake", pridict_runner=_runner(), ngrna_max_offtargets=1
    )
    # PRIDICT's own nicking guide was usable under the tolerance, so it is the one chosen.
    assert out.result.ngrna_protospacer == W._MARKER_B
    assert "1 off-target site(s)" in dict(out.result.gate_results)["offtarget"]
    assert any("second-nick guide has 1 off-target" in w for w in out.result.warnings)


def test_the_cli_accepts_the_setting_from_a_config_file():
    """The CLI keeps only config keys that are design_recoded parameters, so being in the
    signature is what makes ``ngrna_max_offtargets`` settable from the web UI and a JSON config."""
    import inspect

    assert "ngrna_max_offtargets" in inspect.signature(recoded_driver.design_recoded).parameters
    cfg = json.loads(Path("examples/ctcf_cterm/config.json").read_text())
    assert "ngrna_max_offtargets" not in cfg, "the shipped examples keep the strict default"


def test_fake_mode_stays_clean_by_default():
    scorers = make_scorers(mode="fake")
    assert scorers.offtarget_scorer(["ACGT" * 5])["ACGT" * 5][1] == 0
