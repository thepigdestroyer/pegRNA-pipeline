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

from pegrna_pipeline.hpc.gate_thresholds import DOCKING_REFERENCE, THRESHOLDS  # noqa: E402
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
}


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
    for _, neutral, nonbinder in DOCKING_REFERENCE:
        assert neutral - nonbinder > 0.15, "the classes must stay well separated"
    neutrals = [row[1] for row in DOCKING_REFERENCE]
    assert neutrals == sorted(neutrals, reverse=True), "ipTM falls with RNA length"


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
        alphabet = set("ACGTU") if r["class"] == "designed" else set("ACGT")
        assert set(rna) <= alphabet and len(rna) == r["len"]
        backbone = r["spacer"] + scaffold
        if r["class"] in ("neutral", "sequestering", "control-canonical"):
            assert rna[:96] == backbone, "built on the functional backbone"
        if r["class"] in ("nonbinder", "control-scrambled"):
            assert rna[:96] != backbone and sorted(rna[:96]) == sorted(
                backbone
            ), "a scrambled backbone has the same composition in a different order"


# ---- the gate against the evidence ------------------------------------------------------


def test_leave_one_out_neutral_accepted_and_non_binders_rejected():
    loo = dr.leave_one_out(RECORDS, FRACTION)
    neutral = [ok for r, _, ok in loo if dr.kind_of(r) == "neutral"]
    nonbinder = [ok for r, _, ok in loo if dr.kind_of(r) == "nonbinder"]
    assert len(neutral) >= 50 and len(nonbinder) >= 35, "enough sequences to mean something"
    assert sum(neutral) / len(neutral) >= 0.95, "too many normal designs rejected (too harsh)"
    assert 1 - sum(nonbinder) / len(nonbinder) >= 0.95, "too many non-binders accepted"


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
    loo = dr.leave_one_out(RECORDS, FRACTION)
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
    extension and every sweep construct pairs at least twenty, though the gate accepts most."""
    neutral = [
        r["scaffold_paired"]
        for r in RECORDS
        if r["class"] == "neutral" and r.get("scaffold_paired") is not None
    ]
    sweep = [
        r["scaffold_paired"]
        for r in RECORDS
        if r["class"] == "sequestering" and r.get("scaffold_paired") is not None
    ]
    assert len(neutral) >= 25 and len(sweep) >= 10
    assert max(neutral) < min(sweep)
