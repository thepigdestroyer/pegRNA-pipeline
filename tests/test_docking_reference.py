"""The docking reference table is a function of the measurements, and the gate is tested on them.

``DOCKING_REFERENCE`` is meant to be derived from ``scripts/analysis/data/docking_calibration.jsonl``
(by ``scripts/analysis/docking_reference.py``), never typed. These tests fail when the table drifts
from its evidence, and they pin the figures that ADR 0024 quotes, so a sentence in the ADR cannot
go on being true only because nobody re-checked it.
"""

import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))
import docking_reference as dr  # noqa: E402

from pegrna_pipeline.hpc.gate_thresholds import (  # noqa: E402
    DOCKING_CALIBRATED_RANGE,
    DOCKING_REFERENCE,
    THRESHOLDS,
)
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD  # noqa: E402

RECORDS = dr.load()
FRACTION = THRESHOLDS.docking_pass_fraction
KNOWN_CLASSES = {
    "neutral",
    "neutral-shuffled",
    "nonbinder",
    "sequestering",
    "designed",
    "control-canonical",
    "control-scrambled",
    "control-scaffold-only",
    "candidate",  # the pipeline's own pool, measured by sample_pool_candidates.py
}
RANGE_TOP = DOCKING_CALIBRATED_RANGE[1]


# ---- the table is the data --------------------------------------------------------------


def test_the_table_in_the_source_is_the_one_the_data_implies():
    assert dr.stale(RECORDS) == [], (
        "gate_thresholds.DOCKING_REFERENCE no longer matches the measurements. Regenerate it: "
        "python scripts/analysis/docking_reference.py --emit"
    )


def test_every_row_rests_on_enough_independent_sequences():
    for (kind, ext), values in dr.cells(dr.reference_records(RECORDS)).items():
        assert len(values) >= 5, f"{kind} at {96 + ext} nt rests on only {len(values)} sequence(s)"


def test_the_table_is_ordered_and_neutral_stays_above_non_binder():
    lengths = [row[0] for row in DOCKING_REFERENCE]
    assert lengths == sorted(lengths)
    for length, neutral, nonbinder in DOCKING_REFERENCE:
        assert neutral > nonbinder
        if length <= RANGE_TOP:
            assert neutral - nonbinder > 0.15, "inside the trusted range the classes stay apart"
    neutrals = [row[1] for row in DOCKING_REFERENCE]
    assert neutrals == sorted(neutrals, reverse=True), "ipTM falls with RNA length"


def test_the_gap_between_the_classes_closes_with_length():
    """The reason the range ends where it does: the span shrinks until it is inside the noise."""
    spans = [n - b for _, n, b in DOCKING_REFERENCE]
    assert spans == sorted(spans, reverse=True)
    assert spans[-1] < 0.5 * spans[0]


def test_the_calibrated_range_is_where_the_data_stop_supporting_a_verdict():
    """Every row up to the top is trusted, and the next row is not (or there is none)."""
    quality = dr.row_quality(RECORDS)
    top = DOCKING_CALIBRATED_RANGE[1]
    for length, _span, _sd, separation, n_neutral, n_nonbinder in quality:
        if length <= top:
            assert separation >= dr.MIN_SEPARATION_SD, (length, separation)
            assert min(n_neutral, n_nonbinder) >= dr.MIN_SEQUENCES_PER_CLASS
    beyond = [q for q in quality if q[0] > top]
    assert beyond, "the range is only meaningful if longer rows were measured"
    first = beyond[0]
    errors = sum(
        ok != (dr.kind_of(r) == "neutral")
        for r, _f, ok in dr.leave_one_out(RECORDS, FRACTION)
        if r["len"] == first[0]
    )
    assert (
        first[3] < dr.MIN_SEPARATION_SD or errors > 0
    ), f"the {first[0]} nt row looks trustworthy: the range could extend"


# ---- the records are sound --------------------------------------------------------------


def test_records_are_well_formed():
    for r in RECORDS:
        assert r["class"] in KNOWN_CLASSES, r["class"]
        assert len(r["samples"]) == 5, "every measurement is the mean of 5 diffusion samples"
        assert r["mean"] == pytest.approx(st.fmean(r["samples"]), abs=1e-3)
        assert all(0.0 <= v <= 1.0 for v in r["samples"])
        if r.get("ext_len") is not None:
            assert r["len"] == 96 + r["ext_len"], r


def test_no_sequence_is_counted_twice():
    seqs = [r["rna"] for r in RECORDS if r.get("rna")]
    assert len(seqs) == len(set(seqs)), "a duplicated RNA would count as two independent sequences"


def test_recorded_sequences_are_what_their_class_says():
    scaffold = STANDARD_SGRNA_SCAFFOLD
    for r in RECORDS:
        rna = r.get("rna")
        if not rna:
            continue
        # The pipeline's own pegRNAs end in a UUU tail written as RNA letters (the NIM client
        # transcribes T to U anyway); everything built for calibration is plain DNA letters.
        alphabet = set("ACGTU") if r["class"] in ("designed", "candidate") else set("ACGT")
        assert set(rna) <= alphabet and len(rna) == r["len"]
        backbone = r["spacer"] + scaffold
        if r["class"] in ("neutral", "sequestering", "control-canonical"):
            assert rna[:96] == backbone, "built on the functional backbone"
        if r["class"] in ("nonbinder", "control-scrambled"):
            assert rna[:96] != backbone and sorted(rna[:96]) == sorted(
                backbone
            ), "a scrambled backbone has the same composition in a different order"


# ---- the gate against the evidence ------------------------------------------------------


def _in_range(loo):
    return [x for x in loo if x[0]["len"] <= RANGE_TOP]


def test_leave_one_out_neutral_accepted_and_non_binders_rejected_in_the_trusted_range():
    loo = _in_range(dr.leave_one_out(RECORDS, FRACTION))
    neutral = [ok for r, _, ok in loo if dr.kind_of(r) == "neutral"]
    nonbinder = [ok for r, _, ok in loo if dr.kind_of(r) == "nonbinder"]
    assert len(neutral) >= 60 and len(nonbinder) >= 45, "enough sequences to mean something"
    assert all(neutral), "a normal design was rejected: too harsh"
    assert not any(nonbinder), "a non-binder was accepted"


def test_beyond_the_range_the_verdicts_really_are_unreliable():
    """The other half of the claim: past the range leave-one-out does make errors."""
    beyond = [x for x in dr.leave_one_out(RECORDS, FRACTION) if x[0]["len"] > RANGE_TOP]
    assert len(beyond) >= 20
    assert any(ok != (dr.kind_of(r) == "neutral") for r, _f, ok in beyond)


def test_a_held_out_spacer_is_judged_correctly():
    held = [r for r in RECORDS if r["source"] in dr.HELD_OUT_SOURCES and dr.kind_of(r)]
    assert len(held) >= 10
    for kind, want in (("neutral", True), ("nonbinder", False)):
        verdicts = [
            dr.evaluate(r, DOCKING_REFERENCE, FRACTION)[1] for r in held if dr.kind_of(r) == kind
        ]
        assert verdicts and sum(v == want for v in verdicts) / len(verdicts) >= 0.9, kind


def test_the_held_out_spacer_really_is_a_different_spacer():
    primary = {r["spacer"] for r in RECORDS if r["source"] not in dr.HELD_OUT_SOURCES}
    held = {r["spacer"] for r in RECORDS if r["source"] in dr.HELD_OUT_SOURCES}
    assert held and not (held & primary)


def test_the_pipelines_own_designs_pass():
    designed = [r for r in RECORDS if r["class"] == "designed"]
    assert len(designed) >= 2
    for r in designed:
        assert dr.evaluate(r, DOCKING_REFERENCE, FRACTION)[1], r["len"]


def test_the_cutoff_sits_between_the_closest_calls_and_the_band_contains_them():
    """The citation for the borderline margin says the closest calls on either side lie inside
    it. Pin that, so the margin cannot silently stop covering them."""
    loo = _in_range(dr.leave_one_out(RECORDS, FRACTION))
    lowest_neutral = min(f for r, f, _ in loo if dr.kind_of(r) == "neutral")
    highest_nonbinder = max(f for r, f, _ in loo if dr.kind_of(r) == "nonbinder")
    assert highest_nonbinder < FRACTION < lowest_neutral, "no overlap, so the cutoff separates"
    margin = THRESHOLDS.docking_borderline_margin
    assert abs(lowest_neutral - FRACTION) < margin
    assert abs(highest_nonbinder - FRACTION) < margin


# ---- what the gate does and does not see (the limits ADR 0024 states) ------------------


def _sequestering_by_window():
    out = defaultdict(list)
    for r in RECORDS:
        if r["class"] == "sequestering":
            a, b = map(int, r["segment"].split(":"))
            out[(a, b)].append(dr.evaluate(r, DOCKING_REFERENCE, FRACTION))
    return out


def test_extensive_sequestration_of_the_five_prime_scaffold_is_rejected():
    for (a, b), results in _sequestering_by_window().items():
        if b - a >= 40 and a < 30:
            assert all(not ok for _, ok in results), f"complement of scaffold[{a}:{b}]"


def test_most_localized_24nt_pairings_are_not_detected():
    """The limitation, pinned: if the gate ever becomes sensitive to these, the docstring in
    complex_dock.py, the citation and ADR 0024 all need to change with it."""
    windows = [v for (a, b), v in _sequestering_by_window().items() if b - a == 24]
    results = [ok for rs in windows for _, ok in rs]
    assert len(results) >= 10
    assert sum(results) / len(results) >= 0.5, "the gate now rejects most 24-nt pairings"


def test_the_three_prime_hairpin_is_invisible_to_the_gate():
    for f, ok in _sequestering_by_window()[(52, 76)]:
        assert ok and f > 0.8


def test_a_folding_check_would_separate_what_the_gate_cannot():
    """Evidence for the follow-up (extension-to-scaffold pairing is checked by nothing): in
    ViennaRNA's MFE structure neutral designs pair at most a dozen scaffold nucleotides with the
    extension and every sweep construct pairs at least twenty, though the gate accepts most.
    (Within the trusted range. Past it, longer random extensions pair more by chance under MFE,
    up to 24 nt at 296 nt, which is why the advisory itself uses the stricter p >= 0.99 measure;
    that one does not grow with length, see test_scaffold_pairing_advisory.)"""
    neutral = [
        r["scaffold_paired"]
        for r in RECORDS
        if r["class"] == "neutral"
        and r.get("scaffold_paired") is not None
        and r["len"] <= RANGE_TOP
    ]
    sweep = [
        r["scaffold_paired"]
        for r in RECORDS
        if r["class"] == "sequestering" and r.get("scaffold_paired") is not None
    ]
    assert len(neutral) >= 25 and len(sweep) >= 10
    assert max(neutral) < min(sweep)


# ---- the pipeline's own candidates: a test of the gate, not part of its table ------------------
#
# The reference was built from random and scrambled extensions. The question that matters for "is
# the gate too harsh?" is what happens to the pegRNAs the pipeline actually produces, so 59 of
# them (CTCF: 8 cargos, 2 spacers, all four 5' nicks) were measured by
# scripts/analysis/sample_pool_candidates.py.


def _candidates():
    return [r for r in RECORDS if r["class"] == "candidate"]


def _summary():
    return dr.candidate_summary(_candidates(), DOCKING_REFERENCE, FRACTION, RANGE_TOP)


def test_the_real_pool_sample_is_large_and_spans_the_lengths_users_will_make():
    c = _candidates()
    assert len(c) >= 50
    assert min(r["len"] for r in c) <= 200 and max(r["len"] for r in c) >= 300


def test_real_candidates_are_not_part_of_the_reference_they_test():
    assert not [r for r in dr.reference_records(RECORDS) if r["class"] == "candidate"]


def test_every_real_candidate_inside_the_range_is_accepted():
    s = _summary()
    assert s["inside"] >= 15, "enough of them inside the range to mean something"
    assert s["accepted"] == s["inside"], "the gate rejected a design the pipeline itself produced"


def test_real_designs_score_at_least_as_well_as_neutral_ones():
    """They sit ABOVE the random-extension baseline, so that baseline is conservative for real
    designs rather than flattering."""
    assert st.median(_summary()["f_values"]) >= 1.0


def test_the_old_floor_would_have_rejected_nearly_every_real_candidate():
    """The measured answer to "was the 0.6 floor too harsh?": 1 of 59 real candidates clears it."""
    s = _summary()
    assert s["old_floor_pass"] <= 0.1 * s["n"], s["old_floor_pass"]
    assert s["accepted"] == s["inside"]


def test_the_first_candidates_production_would_try_pass():
    head = _summary()["head"]
    assert len(head) >= 5
    assert all(ok for _r, _f, ok in head)


def test_candidates_beyond_the_range_are_reported_not_judged():
    """41 of the 59 are longer than the range. They pass as advisory-only, with the reason."""
    from pegrna_pipeline.hpc.complex_dock import passes_calibrated_docking_gate

    beyond = [r for r in _candidates() if r["len"] > RANGE_TOP]
    assert len(beyond) >= 20

    class Recorded:
        def __init__(self, samples):
            self._s = samples

        def interface_iptm_samples(self, pegrna, cas9):  # noqa: ARG002
            return self._s

    for r in beyond:
        metrics = {}
        ok, why = passes_calibrated_docking_gate(
            "A" * r["len"], "nCas9", scorer=Recorded(r["samples"]), metrics=metrics
        )
        assert ok and why.startswith("ADVISORY ONLY") and metrics["advisory"]
        assert "converge" in why, "the message must say why it is not judged"


def test_candidate_records_say_where_they_came_from():
    for r in _candidates():
        assert r["locus"] and r["cargo"] and r["cargo_len"] > 0
        assert 0 <= r["pool_rank"] < r["pool_size"]
        assert r["rna"].startswith(r["spacer"]) and r["len"] == len(r["rna"])
