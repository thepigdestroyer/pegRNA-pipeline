"""Only protein-coding bases are ever recoded, and every change is reported (ADR 0029).

The CTCF cassette ends in a T7 promoter read in frame 0 after the ORF's stop codon. The recoder's guard counted
the promoter's own stop codons as "internal" and switched the lever off for the whole cargo; relaxing it naively
would have started recoding the promoter, whose consensus cannot change without breaking it. These tests pin the
lock model: what is locked, why, that nothing locked can move, and that every nucleotide that does move is listed.
"""

import dataclasses

import pytest
import test_pridict_ngrna_wiring as W

from pegrna_pipeline.hpc import recoded_driver
from pegrna_pipeline.hpc import recoded_tjpe_workflow as rw
from pegrna_pipeline.hpc.cargo_spec import (
    CODE,
    MIN_ORF_CODONS,
    arm_changes,
    assert_locks_hold,
    build_cargo_spec,
    build_change_report,
    cargo_variants,
    diff_cargo,
    free_codons,
    library_part_is_coding,
)
from pegrna_pipeline.hpc.insert_elements import RC_T7G_PROMOTER, T7G_PROMOTER
from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA, run_recoded_tjpe
from pegrna_pipeline.hpc.report_writer import render_changes_csv

ORF = "CTGAAGGAGCTG" * 10  # 40 codons, LKEL x 10: easy residues, so the CAI generator can recode it
STOP = "TGA"
OLD_T7 = "TAATACGACTCACTATAGGGG"  # the 4-G tail of the shipped cassette


def _aa(dna):
    return "".join(CODE[dna[i : i + 3]] for i in range(0, len(dna) - len(dna) % 3, 3))


def _rc(dna):
    return dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _part(spec, name):
    return next(p for p in spec.parts if p.name == name)


# ---- what is locked, with no annotation ----------------------------------------------------------------


def test_the_ctcf_shape_orf_then_stop_then_t7_locks_the_stop_and_the_promoter():
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER)
    orf, stop, t7 = spec.parts
    assert (orf.name, orf.locked, orf.start, orf.end) == ("protein-coding ORF", False, 0, 120)
    assert (stop.name, stop.locked, stop.start, stop.end) == ("stop codon", True, 120, 123)
    assert (t7.name, t7.locked, t7.start, t7.end) == ("T7 promoter", True, 123, 142)
    assert len(free_codons(spec)) == 40


def test_the_old_four_g_tail_is_locked_whole_including_the_g_run():
    spec = build_cargo_spec(ORF + STOP + OLD_T7)
    t7 = _part(spec, "T7 promoter")
    assert (t7.start, t7.end) == (123, 123 + len(OLD_T7)) and t7.locked


def test_the_reverse_promoter_is_found_with_its_c_run():
    spec = build_cargo_spec(ORF + STOP + RC_T7G_PROMOTER)
    rev = _part(spec, "T7 promoter (reverse)")
    assert (rev.start, rev.end) == (123, 123 + len(RC_T7G_PROMOTER)) and rev.locked


def test_a_promoter_at_the_start_leaves_nothing_to_recode_and_says_so():
    spec = build_cargo_spec(RC_T7G_PROMOTER + ORF)
    assert free_codons(spec) == []
    assert all(p.locked for p in spec.parts)
    assert any("Nothing is recoded" in n for n in spec.notes)


def test_a_promoter_inside_the_reading_frame_clips_the_orf_at_its_start():
    spec = build_cargo_spec(ORF + T7G_PROMOTER + ORF)
    orf = _part(spec, "protein-coding ORF")
    assert (orf.start, orf.end) == (0, 120)
    assert _part(spec, "T7 promoter").locked
    assert all(p.locked for p in spec.parts if p.start >= 120)  # everything after it: kept exactly


def test_a_reverse_promoter_inside_the_frame_is_never_recoded():
    """The reverse T7 has no stop codon in frame 0, so only recognising it protects it: without that the ORF
    would run straight through it and into the next one, and the lever would rewrite the promoter.
    """
    cargo = ORF + RC_T7G_PROMOTER + ORF + STOP
    spec = build_cargo_spec(cargo)
    assert (_part(spec, "protein-coding ORF").end, _part(spec, "T7 promoter (reverse)").locked) == (
        120,
        True,
    )
    variants = cargo_variants(spec, 6)
    assert len(variants) > 1
    start = 120
    assert all(v[start : start + len(RC_T7G_PROMOTER)] == RC_T7G_PROMOTER for v in variants)
    assert all(
        v[start + len(RC_T7G_PROMOTER) :] == cargo[start + len(RC_T7G_PROMOTER) :] for v in variants
    )


def test_a_clean_cds_of_any_length_is_coding_as_before():
    flag = "GACTACAAGGACGACGATGACAAG"
    spec = build_cargo_spec(flag)
    assert len(free_codons(spec)) == 8 and not spec.parts[0].locked


def test_a_cds_with_one_terminal_stop_is_coding_and_keeps_its_stop():
    spec = build_cargo_spec("GACTACAAGGACGACGATGACAAG" + "TAA")
    assert len(free_codons(spec)) == 8
    assert _part(spec, "stop codon").locked


def test_a_short_orf_followed_by_junk_is_not_taken_for_protein():
    short = "CTGAAGGAGCTG" * 3  # 12 codons
    spec = build_cargo_spec(short + STOP + "ACGTACGTAC")
    assert free_codons(spec) == [] and any(str(MIN_ORF_CODONS) in n for n in spec.notes)


def test_a_long_orf_with_a_trailing_tail_of_any_length_is_coding():
    spec = build_cargo_spec(ORF + STOP + "ACGTACGTACGT" + "A")  # 1 nt over a multiple of 3
    assert len(free_codons(spec)) == 40 and _part(spec, "sequence outside the ORF").locked


def test_the_library_rule_for_what_is_protein():
    assert library_part_is_coding(aa="DYKDDDDK", category="epitope_tag")
    assert not library_part_is_coding(aa="", category="regulatory")  # Kozak
    assert not library_part_is_coding(aa="*", category="stop")
    assert not library_part_is_coding(aa="STOP", category="stop")  # a stop entry is never recoded
    assert not library_part_is_coding(aa="", category="amplification")  # the T7 elements


# ---- explicit parts and excluded stretches ---------------------------------------------------------------


def test_marked_parts_are_taken_at_their_word_and_unmarked_means_untouched():
    parts = [
        {"name": "FLAG", "dna": "GACTACAAGGACGACGATGACAAG", "locked": False},
        {"name": "my promoter", "dna": "TTGACAATTAATCATCGG"},  # no flag: locked
        {"name": "ZeoR", "dna": ORF, "locked": False},
    ]
    spec = build_cargo_spec("".join(p["dna"] for p in parts), parts=parts)
    assert [p.locked for p in spec.parts] == [False, True, False]
    assert len(free_codons(spec)) == 8 + 40


def test_a_part_marked_coding_that_is_not_a_clean_cds_is_kept_exactly_with_a_note():
    parts = [
        {"name": "odd", "dna": "ACGTA", "locked": False},
        {"name": "stopped", "dna": "ATGTAAATG", "locked": False},
    ]
    spec = build_cargo_spec("ACGTAATGTAAATG", parts=parts)
    assert all(p.locked for p in spec.parts) and len(spec.notes) == 2
    assert "not a multiple of 3" in spec.notes[0] and "stop codon before its last" in spec.notes[1]


def test_parts_that_do_not_add_up_to_the_cargo_are_refused():
    with pytest.raises(ValueError, match="do not add up"):
        build_cargo_spec(ORF, parts=[{"name": "a", "dna": ORF[:-3], "locked": False}])


def test_two_coding_parts_in_different_frames_get_a_note():
    parts = [
        {"name": "A", "dna": ORF, "locked": False},
        {"name": "T7-G", "dna": T7G_PROMOTER, "locked": True},  # 19 nt: moves the frame by 1
        {"name": "B", "dna": ORF, "locked": False},
    ]
    spec = build_cargo_spec("".join(p["dna"] for p in parts), parts=parts)
    assert any("read in a different frame from A" in n for n in spec.notes)


def test_excluded_ranges_are_one_based_inclusive_and_merge():
    spec = build_cargo_spec(ORF, locked=[[1, 3], [3, 6], [10, 12]])
    assert spec.excluded == ((0, 6), (9, 12))


@pytest.mark.parametrize("bad", [[[0, 3]], [[5, 4]], [[1, 999]], [["a", 2]], [[1]]])
def test_a_bad_excluded_range_is_refused_with_the_reason(bad):
    with pytest.raises(ValueError, match="cargo_locked"):
        build_cargo_spec(ORF, locked=bad)


def test_a_codon_that_touches_an_excluded_base_is_kept_whole_and_its_neighbours_are_free():
    spec = build_cargo_spec(ORF, locked=[[5, 5]])  # nt 5 is in codon 2 (nt 4-6)
    free = free_codons(spec)
    assert 3 not in free and 0 in free and 6 in free and len(free) == 39


def test_an_excluded_stretch_is_reported_as_kept_exactly():
    spec = build_cargo_spec(ORF + STOP, locked=[[7, 15]])
    assert any(
        line.startswith("nt 7-15 (9 nt): excluded by you")
        for line in build_change_report(spec, ORF + STOP, 0).kept_exactly
    )


# ---- nothing locked can move -------------------------------------------------------------------------------


def _check_every_variant(spec):
    variants = cargo_variants(spec, 6)
    assert (
        variants[0] == spec.sequence and len(variants) > 1
    ), "the generator should find variants here"
    mask = spec.locked_mask()
    for v in variants[1:]:
        assert len(v) == len(spec.sequence) and v != spec.sequence
        assert all(a == b for a, b, m in zip(spec.sequence, v, mask) if m), "a locked base changed"
        assert _aa(v) == _aa(spec.sequence), "an amino acid changed"
    for part in spec.parts:  # a codon that touches a locked base is kept WHOLE, not just that base
        for s in range(part.start, part.end - 2, 3):
            if not part.locked and any(mask[s : s + 3]):
                assert all(v[s : s + 3] == spec.sequence[s : s + 3] for v in variants), s
    return variants


def test_no_variant_ever_touches_the_stop_or_the_promoter():
    for tail in (T7G_PROMOTER, OLD_T7, RC_T7G_PROMOTER):
        spec = build_cargo_spec(ORF + STOP + tail)
        variants = _check_every_variant(spec)
        assert all(v.endswith(STOP + tail) for v in variants)


def test_no_variant_touches_a_user_excluded_stretch_or_the_codons_it_touches():
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER, locked=[[10, 20], [61, 62]])
    for v in _check_every_variant(spec):
        assert v[9:21] == spec.sequence[9:21]  # nt 10-20 plus the codons around them kept whole
        assert v[57:63] == spec.sequence[57:63]


def test_a_worst_case_recoder_still_cannot_move_a_locked_codon(monkeypatch):
    """The real generator seldom touches these codons (the test ORF is already CAI-optimal), so use one that
    rewrites EVERY codon it is handed. The expectation is hard-coded, not read back from free_codons.
    """
    from pegrna_pipeline.hpc import codon_expander

    swap = {"L": "TTA", "K": "AAA", "E": "GAA"}  # none of these is the codon the ORF uses
    monkeypatch.setattr(
        codon_expander,
        "rank_insert_variants_by_cai",
        lambda aa, **kw: [dataclasses.make_dataclass("V", ["dna"])("".join(swap[a] for a in aa))],
    )
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER, locked=[[10, 20], [61, 62]])
    original, variant = cargo_variants(spec, 6)
    kept_whole = {9, 12, 15, 18, 60}  # the codons that touch nt 10-20 and nt 61-62 (0-based starts)
    for start in range(0, 120, 3):
        same = variant[start : start + 3] == original[start : start + 3]
        assert same == (start in kept_whole), start
    assert variant[120:] == original[120:]  # the stop and the T7: untouched


def test_no_variant_when_there_is_no_free_codon():
    spec = build_cargo_spec(ORF, parts=[{"name": "mine", "dna": ORF, "locked": True}])
    assert cargo_variants(spec, 6) == [ORF]


def test_a_cap_of_one_disables_the_lever():
    assert cargo_variants(build_cargo_spec(ORF), 1) == [ORF]


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda s: s[:123] + "A" + s[124:], "locked base 124"),  # inside the promoter
        (lambda s: s[:120] + "TAA" + s[123:], "locked base"),  # the stop codon
        (lambda s: s[:3] + "GAA" + s[6:], "changed amino acid"),  # AAG (K) -> GAA (E)
        (lambda s: s + "A", "keep its length"),
    ],
)
def test_the_guard_refuses_any_change_to_a_locked_base_a_stop_an_amino_acid_or_the_length(
    mutate, message
):
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER)
    with pytest.raises(AssertionError, match=message):
        assert_locks_hold(spec, mutate(spec.sequence))


# ---- exactly what changed ----------------------------------------------------------------------------------------


def test_the_diff_lists_every_changed_nucleotide_with_its_part_position_and_codon():
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER)
    new = ORF[:3] + "AAA" + ORF[6:] + STOP + T7G_PROMOTER  # codon 2: AAG -> AAA (Lys)
    (c,) = diff_cargo(spec, new)
    assert (c.scope, c.part, c.position, c.codon) == ("cargo", "protein-coding ORF", 6, 2)
    assert (c.old_codon, c.new_codon, c.aa, c.old, c.new) == ("AAG", "AAA", "K", "G", "A")
    assert "protein-coding ORF, nt 6 (codon 2, K): AAG -> AAA  [G->A]" == c.line()


def test_the_number_of_listed_changes_is_exactly_the_hamming_distance():
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER)
    for v in cargo_variants(spec, 6)[1:]:
        assert len(diff_cargo(spec, v)) == sum(a != b for a, b in zip(spec.sequence, v))


def test_a_cargo_used_as_given_says_so_and_still_lists_what_was_protected():
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER)
    report = build_change_report(spec, spec.sequence, 0)
    text = "\n".join(report.lines())
    assert report.cargo_unchanged and "used exactly as you gave it (142 nt, 0 changed)" in text
    assert (
        "kept exactly: T7 promoter (nt 124-142, 19 nt)" in text
        and "kept exactly: stop codon (nt 121-123" in text
    )
    assert text.endswith("Everything not listed above is exactly as you specified.")


def test_a_recoded_cargo_names_the_parts_and_every_nucleotide():
    spec = build_cargo_spec(ORF + STOP + T7G_PROMOTER)
    variant = cargo_variants(spec, 6)[1]
    report = build_change_report(spec, variant, 1)
    text = "\n".join(report.lines(limit=500))
    assert (
        not report.cargo_unchanged and "RECODED, not as you gave it" in text and "variant 1" in text
    )
    assert text.count("protein-coding ORF, nt ") == len(report.cargo_changes)
    assert "T7 promoter, nt" not in text  # the promoter is listed as kept, never as changed


def test_the_rows_and_the_dict_carry_the_same_changes():
    spec = build_cargo_spec(ORF + STOP)
    report = build_change_report(spec, cargo_variants(spec, 6)[1], 1)
    assert [r["position"] for r in report.rows()] == [c.position for c in report.cargo_changes]
    d = report.to_dict()
    assert d["cargo_unchanged"] is False and d["cargo_variant"] == 1
    assert len(d["cargo_changes"]) == len(report.cargo_changes) and d["arm_changes"] == []


def test_arm_changes_report_window_position_codon_and_amino_acid():
    genome = "AAA" + "CTG" * 6 + "CCC"  # CDS from nt 4 (1-based), Leu x 6
    mask = [False] * 3 + [True] * 18 + [False] * 3
    arm = "CTC" + "CTG" * 2  # replaces genome[3:12): the first codon CTG -> CTC, still Leu
    (c,) = arm_changes(genome, 3, 12, arm, mask)
    assert (c.scope, c.part, c.position, c.codon) == (
        "genomic arm",
        "replacement arm",
        6,
        1,
    )  # 3rd base
    assert (c.old_codon, c.new_codon, c.aa, c.old, c.new) == ("CTG", "CTC", "L", "G", "C")


def test_a_non_coding_arm_base_has_no_codon():
    genome = "ACGTACGTAC"
    (c,) = arm_changes(
        genome, 2, 6, "GTAT", [False] * 10
    )  # GTAC -> GTAT: one change, non-coding DNA
    assert (c.position, c.codon, c.old_codon, c.aa) == (6, None, "", "")


def test_an_arm_of_the_wrong_length_is_refused():
    with pytest.raises(ValueError, match="as long as"):
        arm_changes("ACGTACGT", 0, 4, "ACG", [True] * 8)


# ---- through the walk ---------------------------------------------------------------------------------------------


def _walk(cargo, **over):
    return run_recoded_tjpe(**W._wf_inputs(cargo=cargo, **over))


class _FailsOnTheOriginal:
    """A splice stand-in that rejects exactly the cargo as given, so only a recoded variant can pass."""

    binding = True
    context_note = ""

    def __init__(self, original):
        self._original = original

    def delta(self, ref_seq, alt_seq):  # noqa: ARG002
        return 0.5 if self._original in alt_seq else 0.0


CARGO = ORF + STOP + T7G_PROMOTER


def test_the_walk_attaches_the_report_and_the_arm_changes_are_exactly_the_recode_edits():
    result = _walk(CARGO)
    assert result.changes is not None and result.changes.cargo_unchanged
    assert len(result.arm_changes) == result.recode_edits == len(result.changes.arm_changes)
    assert all(c.old != c.new for c in result.arm_changes)


def test_a_failing_cargo_is_recoded_in_its_orf_only_and_the_report_lists_every_change():
    result = _walk(CARGO, splice_scorer=_FailsOnTheOriginal(CARGO))
    ch = result.changes
    assert ch.cargo_variant >= 1 and not ch.cargo_unchanged
    assert {c.part for c in ch.cargo_changes} == {"protein-coding ORF"}
    assert all(c.position <= len(ORF) for c in ch.cargo_changes)  # nothing at or after the stop
    # the RTT of an Rv pegRNA carries the reverse complement of the cargo: the stop and T7 reach it intact
    assert _rc(result.insert_seq).endswith(STOP + T7G_PROMOTER)
    assert any("cargo was RECODED" in w for w in result.warnings)


def test_a_recoded_variant_is_never_reported_as_the_original():
    result = _walk(CARGO, splice_scorer=_FailsOnTheOriginal(CARGO))
    written = _rc(result.insert_seq)  # the cargo the pegRNA writes, back in sense orientation
    assert written != CARGO and _aa(written[: len(ORF)]) == _aa(ORF)


def test_marking_the_orf_locked_stops_the_lever_and_says_why():
    parts = [
        {"name": "my ORF", "dna": ORF, "locked": True},
        {"name": "rest", "dna": STOP + T7G_PROMOTER},
    ]
    with pytest.raises(NoValidPegRNA):
        _walk(CARGO, splice_scorer=_FailsOnTheOriginal(CARGO), cargo_parts=parts)


def test_an_excluded_range_holds_in_the_designed_pegrna():
    result = _walk(CARGO, splice_scorer=_FailsOnTheOriginal(CARGO), cargo_locked=[[1, 30]])
    assert _rc(result.insert_seq)[:30] == CARGO[:30]
    assert all(c.position > 30 for c in result.changes.cargo_changes)


def test_a_lever_that_cannot_generate_a_variant_says_so_in_the_failure(monkeypatch):
    monkeypatch.setattr(rw, "cargo_variants", lambda spec, n: [spec.sequence])
    with pytest.raises(NoValidPegRNA) as exc:
        _walk(CARGO, splice_scorer=_FailsOnTheOriginal(CARGO))
    assert any("no synonymous cargo variant could be generated" in d for d in exc.value.diagnostics)
    assert "no synonymous cargo variant could be generated" in str(
        exc.value
    )  # the text people read


def test_a_bad_lock_input_fails_before_any_design_work():
    with pytest.raises(ValueError, match="cargo_locked"):
        _walk(CARGO, cargo_locked=[[1, 9999]])


# ---- the files and the JSON the UI reads ---------------------------------------------------------------------------


def _design(tmp_path, **kw):
    return recoded_driver.design_recoded(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo=CARGO,
        genome=W.GENOME,
        site_s=W.SITE_S,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
        mode="fake",
        output_dir=str(tmp_path),
        **kw,
    )


def test_the_report_has_a_changes_section(tmp_path):
    out = _design(tmp_path)
    assert "Changes from your input:" in out.report_text
    assert "Cargo: used exactly as you gave it (142 nt, 0 changed)." in out.report_text
    assert "kept exactly: T7 promoter (nt 124-142, 19 nt)" in out.report_text
    assert "Everything not listed above is exactly as you specified." in out.report_text


def test_the_changes_csv_has_one_row_per_changed_nucleotide(tmp_path):
    out = _design(tmp_path)
    rows = (tmp_path / "recoded_tjpe_changes.csv").read_text().strip().splitlines()
    assert rows[0] == "scope,part,position,codon,old_codon,new_codon,amino_acid,old_base,new_base"
    assert len(rows) - 1 == out.result.recode_edits and all(
        r.startswith("genomic arm,") for r in rows[1:]
    )


def test_the_gates_json_carries_the_changes_for_the_ui(tmp_path):
    import json

    _design(tmp_path)
    gates = json.loads((tmp_path / "recoded_tjpe_gates.json").read_text())
    assert gates["changes"]["cargo_unchanged"] is True
    assert len(gates["changes"]["arm_changes"]) == gates["recode_edits"]
    assert any("T7 promoter" in line for line in gates["changes"]["kept_exactly"])


def test_the_driver_passes_the_lock_inputs_through(tmp_path):
    out = _design(tmp_path, cargo_locked=[[1, 12]])
    assert any(
        line.startswith("nt 1-12 (12 nt): excluded by you")
        for line in out.result.changes.kept_exactly
    )


def test_a_result_built_by_hand_still_renders_without_a_changes_object():
    result = dataclasses.replace(_walk(CARGO), changes=None)
    assert (
        render_changes_csv(result).strip()
        == "scope,part,position,codon,old_codon,new_codon,amino_acid,old_base,new_base"
    )
