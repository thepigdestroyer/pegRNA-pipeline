"""Reports, files, driver and command line for ``design-pa`` (spec 4.4).

The report is what a reader orders DNA from, so what it must NEVER do is pass a placeholder off as a verdict, hide
that a donor restoring the span is this project's extension of the method, or make the reader work out which
sequence to order. These tests pin the words as well as the numbers.
"""

import csv
import dataclasses
import io
import json

import pytest
from click.testing import CliRunner
from pa_helpers import fake_scorers, make_cargo, make_locus

from pegrna_pipeline.agent.cli import cli
from pegrna_pipeline.hpc import recoded_pa_driver as driver
from pegrna_pipeline.hpc.pa_report import (
    pa_gates_payload,
    render_pa_csv,
    render_pa_fasta,
    render_pa_report,
)
from pegrna_pipeline.hpc.recoded_pa_driver import PA_FILES, design_pa, pa_stand_in_gates
from pegrna_pipeline.hpc.recoded_pa_workflow import (
    STANDING_CAVEAT,
    NoValidPADesign,
    run_recoded_pa,
)
from pegrna_pipeline.hpc.scorer_factory import GateScorers

S = 150
CARGO = make_cargo(120)


@pytest.fixture(scope="module")
def locus():
    return make_locus(
        [("top", 123, "AGG"), ("top", 147), ("bottom", 164), ("bottom", 190)],
        n=300,
        site_s=S,
        seed=7,
        coding_end=S,
    )


@pytest.fixture(scope="module")
def result(locus):
    return run_recoded_pa(
        cargo=CARGO,
        genome=locus.genome,
        site_s=S,
        cds_mask=locus.cds_mask,
        **fake_scorers(),
    )


MANIFEST = {"scorer_mode": "api", "uniprot": "P49711", "gene": "CTCF"}


# ---- the text report ------------------------------------------------------------------------------------


def test_the_report_names_the_design_and_says_what_to_order(result):
    text = render_pa_report(result, MANIFEST)
    assert text.startswith("=== Prime Assembly design: two pegRNAs and one donor ===")
    assert "WHAT TO ORDER" in text
    for needle in ("pegRNA A", "pegRNA B", "donor DNA"):
        assert needle in text


def test_both_full_pegrnas_and_the_full_donor_are_in_the_report_verbatim(result):
    text = render_pa_report(result, MANIFEST)
    assert result.guide_a.pegRNA_sequence in text
    assert result.guide_b.pegRNA_sequence in text
    assert result.donor in text


def test_each_pegrna_is_broken_into_its_parts(result):
    text = render_pa_report(result, MANIFEST)
    for g in (result.guide_a, result.guide_b):
        assert f"spacer  : {g.spacer}" in text
        assert f"RTT     : {g.rtt}" in text
        assert f"PBS     : {g.pbs}" in text
    assert "writes flap A" in text and "writes flap B" in text
    assert "UUU" in text


def test_the_donor_is_broken_into_its_arms_and_the_cargo_with_the_flap_ends_marked(result):
    text = render_pa_report(result, MANIFEST)
    assert f"cargo  : {result.cargo}" in text
    if result.arm5:
        assert f"arm5   : {result.arm5}" in text
    if result.arm3:
        assert f"arm3   : {result.arm3}" in text
    assert f"5' {result.guide_a.flap_len} nt" in text and f"3' {result.guide_b.flap_len} nt" in text
    assert "match flap A" in text and "match flap B" in text


def test_the_standing_caveat_is_in_the_report_and_cannot_be_missed(result):
    text = render_pa_report(result, MANIFEST)
    assert STANDING_CAVEAT in text
    first_third = text[: len(text) // 3]
    assert (
        "extension of Prime Assembly" in first_third
    ), "it belongs near the top, not buried in the warnings"


def test_a_pegrna_is_rna_and_the_report_says_to_read_t_as_u(result):
    assert "read each T as U" in render_pa_report(result, MANIFEST)


def _prep_section(text):
    return text[text.index("HOW TO PREPARE THE DONOR") : text.index("Coordinates & invariants:")]


def test_how_to_prepare_the_donor_is_given_with_what_it_rests_on(result):
    text = render_pa_report(result, MANIFEST)
    assert "HOW TO PREPARE THE DONOR" in text
    prep = _prep_section(text)
    for needle in (
        "PCR",
        "lambda exonuclease",
        "phosphorothioate",
        "AZD-7648",
        "PE6c",
        "HEK293T",
        "AAVS1",
    ):
        assert needle in prep, needle


def test_the_donor_notes_say_what_the_paper_showed_and_what_it_did_not(result):
    prep = _prep_section(render_pa_report(result, MANIFEST))
    assert "about 50% insertion at AAVS1" in prep and "20-40% insertion" in prep
    assert (
        "slightly less than the plain PCR product" in prep
    ), "the overhang donor was a little LESS efficient"
    assert "restores the DNA between the nicks" in prep and "has to be measured" in prep


def test_the_donor_notes_do_not_make_claims_the_paper_does_not(result):
    """These two were in an earlier draft and are wrong: the phosphorothioate linkages are internal, next to the
    region matching the flap (not at the primers' 5' ends), and the overhang donor was not the more efficient one.
    """
    prep = _prep_section(render_pa_report(result, MANIFEST))
    assert "5' end of each primer" not in prep
    assert "higher efficiency" not in prep and "more efficient" not in prep
    assert (
        "one locus in one cell line" not in prep
    ), "the paper shows several loci and two cell types"


def test_the_coordinates_and_invariants_are_stated(result):
    text = render_pa_report(result, MANIFEST)
    n1, n2 = result.span
    assert f"insert site S      : {result.site_s}" in text
    assert f"{n2 - n1} bp" in text and "deletion size      : 0 (INV-0 asserted zero)" in text
    assert f"recode edits       : {result.recode_edits}" in text
    assert f"pairs considered   : {result.pairs_considered}" in text
    assert f"pairs walked       : {result.candidates_tried}" in text
    assert f"{len(result.product)} nt" in text


def test_the_nicks_are_described_between_numbered_bases_not_as_a_bare_index(result):
    text = render_pa_report(result, MANIFEST)
    n1, n2 = result.span
    assert f"between window nt {n1} and {n1 + 1}" in text
    assert f"between window nt {n2} and {n2 + 1}" in text


def test_the_pam_status_of_each_guide_is_stated_in_words(result):
    text = render_pa_report(result, MANIFEST)
    assert "PAM of guide A" in text and "PAM of guide B" in text


def _block(text, start, stop):
    return text[text.index(start) : text.index(stop)]


def test_each_pegrna_block_says_which_end_of_the_donor_its_flap_is(result):
    text = render_pa_report(result, MANIFEST)
    block_a = _block(text, "  pegRNA A  (", "  pegRNA B  (")
    block_b = _block(text, "  pegRNA B  (", "  donor DNA  (")
    assert f"writes flap A, {result.guide_a.flap_len} nt = the donor's 5' end" in block_a
    assert f"writes flap B, {result.guide_b.flap_len} nt = the donor's 3' end" in block_b
    assert "3' end)" not in block_a.replace("3' end of the nicked", "") and "5' end)" not in block_b


def test_the_coordinates_lines_give_both_nicks_and_the_window_range_of_the_span(result):
    """The per-guide blocks print each nick too, so these are the exact coordinate lines, not a substring."""
    text = render_pa_report(result, MANIFEST)
    n1, n2 = result.span
    assert (
        f"  nicks              : A between window nt {n1} and {n1 + 1}; "
        f"B between window nt {n2} and {n2 + 1}\n"
    ) in text
    assert f"  span replaced      : {n2 - n1} bp (window nt {n1 + 1}-{n2}), restored by" in text


def test_each_donor_arm_says_which_window_range_it_restores(result):
    text = render_pa_report(result, MANIFEST)
    n1, n2 = result.span
    assert result.arm5 and result.arm3, "the fixture must have both arms for this to mean anything"
    assert (
        f"{len(result.arm5)} nt: the restored copy of window nt {n1 + 1}-{result.site_s})" in text
    )
    assert (
        f"{len(result.arm3)} nt: the restored copy of window nt {result.site_s + 1}-{n2})" in text
    )


_PAM_FRAGMENT = {
    "broken": "broken in the donor by one synonymous change, so the guide cannot cut the product",
    "broken-by-recode": "broken by the recoding of the arm, so the guide cannot cut the product",
    "weakened": "WEAKENED: it reads {pam} in the product, which SpCas9 still cuts at a fraction of its NGG rate",
    "intact-noncoding": "REGENERATED: it lies in non-coding DNA, which is kept exactly, so the guide can cut the product again",
    "intact-no-synonymous-change": "REGENERATED: it lies in coding sequence with no synonymous change that breaks it, so the guide can cut the product again",
}
_PAM_OF = {
    "A": "AAG",
    "B": "AGA",
}  # told apart on purpose: each line must carry its own guide's PAM


@pytest.mark.parametrize("status", sorted(_PAM_FRAGMENT))
@pytest.mark.parametrize("side", ["A", "B"])
def test_each_guides_pam_line_describes_that_guides_pam_in_words(result, side, status):
    """Whether the guide can cut the product again is the thing a reader has to act on, so each outcome has its
    own wording and each guide's line carries its own guide's outcome and PAM, not the other's."""
    other = "B" if side == "A" else "A"
    other_status = "broken" if status != "broken" else "intact-noncoding"
    ga = dataclasses.replace(
        result.guide_a,
        pam_status=status if side == "A" else other_status,
        product_pam=_PAM_OF["A"],
    )
    gb = dataclasses.replace(
        result.guide_b,
        pam_status=status if side == "B" else other_status,
        product_pam=_PAM_OF["B"],
    )
    text = render_pa_report(dataclasses.replace(result, guide_a=ga, guide_b=gb), MANIFEST)
    lines = {
        ln.split(":")[0].strip(): ln for ln in text.splitlines() if ln.startswith("  PAM of guide")
    }
    want = _PAM_FRAGMENT[status].replace("{pam}", _PAM_OF[side])
    other_want = _PAM_FRAGMENT[other_status].replace("{pam}", _PAM_OF[other])
    assert want in lines[f"PAM of guide {side}"]
    assert other_want in lines[f"PAM of guide {other}"]
    assert want not in lines[f"PAM of guide {other}"]


def test_every_pam_outcome_the_geometry_can_produce_has_words():
    """The statuses documented on ``PamOutcome``; one added there without words here would print raw."""
    from pegrna_pipeline.hpc import pa_report

    assert (
        set(pa_report._PAM_WORDS)
        == set(_PAM_FRAGMENT)
        == {
            "broken-by-recode",
            "broken",
            "weakened",
            "intact-noncoding",
            "intact-no-synonymous-change",
        }
    )


def test_the_arms_kept_exactly_line_is_under_the_arm_not_under_the_cargo(result):
    text = render_pa_report(result, MANIFEST)
    section = _block(text, "Changes from your input:", "Advisory warnings:")
    lines = section.splitlines()
    cargo_at = next(i for i, ln in enumerate(lines) if ln.strip().startswith("Cargo:"))
    arm_at = next(i for i, ln in enumerate(lines) if ln.strip().startswith("Genomic arm:"))
    kept_arm = [i for i, ln in enumerate(lines) if "kept exactly: donor arm" in ln]
    assert kept_arm and all(i > arm_at > cargo_at for i in kept_arm)


def test_a_pegrna_block_says_when_its_rtt_end_pairs_with_its_spacer_end_and_why_that_is_fine(
    result,
):
    ga = dataclasses.replace(result.guide_a, nick_pairing_bp=3)
    gb = dataclasses.replace(result.guide_b, nick_pairing_bp=0)
    text = render_pa_report(dataclasses.replace(result, guide_a=ga, guide_b=gb), MANIFEST)
    block_a = _block(text, "  pegRNA A  (", "  pegRNA B  (")
    block_b = _block(text, "  pegRNA B  (", "  donor DNA  (")
    assert "the RTT's last 3 nt pair with the spacer's last 3 nt" in block_a
    assert (
        "the donor copies the genome at this nick: expected, and not counted as a misfold"
        in block_a
    )
    assert (
        block_a.index("PBS     :")
        < block_a.index("the RTT's last 3 nt")
        < block_a.index("PolyU   :")
    )
    assert "RTT's last" not in block_b, "a guide with no such pairing gets no note"


def test_the_note_names_the_length_of_the_pairing_it_is_about(result):
    gb = dataclasses.replace(result.guide_b, nick_pairing_bp=2)
    text = render_pa_report(dataclasses.replace(result, guide_b=gb), MANIFEST)
    assert "the RTT's last 2 nt pair with the spacer's last 2 nt" in _block(
        text, "  pegRNA B  (", "  donor DNA  ("
    )


def test_without_such_a_pairing_the_report_has_no_note(result):
    assert result.guide_a.nick_pairing_bp == 0 and result.guide_b.nick_pairing_bp == 0
    assert "RTT's last" not in render_pa_report(result, MANIFEST)


def test_the_payload_carries_each_guides_pairing_length(result):
    ga = dataclasses.replace(result.guide_a, nick_pairing_bp=3)
    gb = dataclasses.replace(result.guide_b, nick_pairing_bp=1)
    guides = pa_gates_payload(dataclasses.replace(result, guide_a=ga, guide_b=gb), [], "api")["pa"][
        "guides"
    ]
    assert [g["nick_pairing_bp"] for g in guides] == [3, 1]


def test_the_changes_section_and_the_warnings_are_present(result):
    text = render_pa_report(result, MANIFEST)
    assert "Changes from your input:" in text
    assert "Cargo: used exactly as you gave it" in text
    assert "Everything not listed above is exactly as you specified." in text
    assert "Advisory warnings:" in text
    for w in result.warnings:
        assert w in text


def test_the_provenance_manifest_is_listed_sorted(result):
    text = render_pa_report(result, MANIFEST)
    tail = text[text.index("Provenance manifest:") :]
    assert tail.index("gene: CTCF") < tail.index("scorer_mode: api") < tail.index("uniprot: P49711")


def test_a_run_with_stand_in_gates_cannot_be_mistaken_for_a_validated_design(result):
    text = render_pa_report(result, {"scorer_mode": "api", "stand_in_gates": "docking"})
    assert text.index("NOT A VALIDATED DESIGN") < text.index("WHAT TO ORDER")
    assert "did NOT run a real model: docking" in text
    assert "Re-run with the missing credentials/resources before ordering DNA." in text


def test_a_fake_mode_run_says_it_is_a_dry_run(result):
    text = render_pa_report(result, {"scorer_mode": "fake", "stand_in_gates": "splice, docking"})
    assert "scorer_mode=fake: this whole run is a dry run of the design logic." in text


def test_a_real_run_has_no_banner(result):
    assert "NOT A VALIDATED DESIGN" not in render_pa_report(result, MANIFEST)


# ---- the components CSV ---------------------------------------------------------------------------------


def _rows(text):
    return [r for r in csv.reader(io.StringIO(text)) if r and not r[0].startswith("#")]


def test_the_csv_has_the_same_columns_as_the_tj_pe_one(result):
    assert _rows(render_pa_csv(result))[0] == [
        "index",
        "component",
        "role",
        "part_of",
        "sequence",
        "length",
    ]


def test_the_csv_rows_are_in_construct_order_with_correct_lengths(result):
    rows = _rows(render_pa_csv(result))[1:]
    assert [int(r[0]) for r in rows] == list(range(1, len(rows) + 1))
    assert all(int(r[5]) == len(r[4]) for r in rows)
    parts = [(r[3], r[1]) for r in rows]
    assert [p for p, _ in parts if p.startswith("pegRNA A")][0] == "pegRNA A"
    order = [p for p, _ in parts]
    assert order.index("pegRNA A") < order.index("pegRNA B") < order.index("donor")


def test_each_pegrna_reassembles_from_its_csv_rows(result):
    rows = _rows(render_pa_csv(result))[1:]
    for side, g in (("A", result.guide_a), ("B", result.guide_b)):
        joined = "".join(r[4] for r in rows if r[3] == f"pegRNA {side}")
        assert joined == g.pegRNA_sequence


def test_the_donor_reassembles_from_its_csv_rows(result):
    rows = _rows(render_pa_csv(result))[1:]
    assert "".join(r[4] for r in rows if r[3] == "donor") == result.donor


def test_the_csv_names_the_donor_parts_and_carries_metadata_rows(result):
    text = render_pa_csv(result, {"gene": "CTCF"})
    names = [r[1] for r in _rows(text)[1:]]
    assert (
        "cargo" in names
        and any(n.startswith("RTT A") for n in names)
        and any(n.startswith("PBS B") for n in names)
    )
    assert "# full_pegRNA_A" in text and "# donor_length" in text and "# gene" in text


# ---- the FASTA -----------------------------------------------------------------------------------------


def test_the_fasta_has_three_records_pegrna_a_pegrna_b_and_the_donor(result):
    records = [r for r in render_pa_fasta(result).split(">") if r.strip()]
    assert len(records) == 3
    heads = [r.splitlines()[0] for r in records]
    assert (
        heads[0].startswith("pegRNA_A")
        and heads[1].startswith("pegRNA_B")
        and heads[2].startswith("donor")
    )
    seqs = ["".join(r.splitlines()[1:]) for r in records]
    assert seqs == [result.guide_a.pegRNA_sequence, result.guide_b.pegRNA_sequence, result.donor]


def test_the_fasta_headers_carry_what_a_reader_needs_to_tell_them_apart(result):
    heads = [r.splitlines()[0] for r in render_pa_fasta(result).split(">") if r.strip()]
    assert f"{len(result.guide_a.pegRNA_sequence)} nt" in heads[0] and "top strand" in heads[0]
    assert f"{len(result.guide_b.pegRNA_sequence)} nt" in heads[1] and "bottom strand" in heads[1]
    assert f"{len(result.donor)} nt" in heads[2]


def test_the_fasta_wraps_long_sequences_without_changing_them(result):
    text = render_pa_fasta(result)
    assert max(len(line) for line in text.splitlines() if not line.startswith(">")) <= 80


# ---- the JSON the web UI reads ---------------------------------------------------------------------------


def test_the_payload_carries_the_generic_fields_the_ui_already_reads(result):
    payload = pa_gates_payload(result, stand_in_gates=["docking"], mode="api")
    assert payload["gate_results"] == [list(x) for x in result.gate_results]
    assert payload["recode_edits"] == result.recode_edits and payload["deletion_size"] == 0
    assert payload["candidates_tried"] == result.candidates_tried
    assert payload["warnings"] == result.warnings
    assert payload["stand_in_gates"] == ["docking"] and payload["scorer_mode"] == "api"
    assert payload["changes"]["cargo_unchanged"] is True
    assert payload["strategy"] == "pa"


def test_the_payload_describes_both_guides_and_the_donor_for_the_result_view(result):
    pa = pa_gates_payload(result, stand_in_gates=[], mode="api")["pa"]
    assert pa["site_s"] == result.site_s and pa["span"] == list(result.span)
    assert pa["cargo_length"] == len(result.cargo) and pa["donor_length"] == len(result.donor)
    assert pa["product_length"] == len(result.product)
    a, b = pa["guides"]
    assert (a["side"], a["strand"], b["side"], b["strand"]) == ("A", "top", "B", "bottom")
    for g, src in ((a, result.guide_a), (b, result.guide_b)):
        assert g["pegRNA"] == src.pegRNA_sequence and g["pegRNA_length"] == len(src.pegRNA_sequence)
        assert (g["flap_len"], g["pbs_len"], g["cut"]) == (src.flap_len, src.pbs_len, src.cut)
        assert g["pam_status"] == src.pam_status and g["spacer"] == src.spacer
        assert g["product_pam"] == src.product_pam and len(g["product_pam"]) == 3
    donor = pa["donor"]
    assert donor["sequence"] == result.donor
    assert [p["name"] for p in donor["parts"]] == [n for n, _ in result.donor_parts]
    assert donor["parts"][0]["start"] == 0 and donor["parts"][-1]["end"] == len(result.donor)
    assert all(p["end"] - p["start"] == len(p["sequence"]) for p in donor["parts"])
    assert donor["flap_a"] == [0, result.guide_a.flap_len]
    assert donor["flap_b"] == [len(result.donor) - result.guide_b.flap_len, len(result.donor)]


def test_the_payload_gives_the_window_length_and_the_donor_notes_the_page_shows(result, locus):
    pa = pa_gates_payload(result, stand_in_gates=[], mode="api")["pa"]
    assert pa["window_length"] == len(locus.genome)
    assert pa["window_length"] == pa["product_length"] - pa["cargo_length"]
    from pegrna_pipeline.hpc import pa_report

    assert pa["donor_prep"] == list(pa_report._DONOR_PREP) and len(pa["donor_prep"]) == 4
    assert all(note in render_pa_report(result, MANIFEST) for note in pa["donor_prep"])


def test_the_payload_is_plain_json(result):
    payload = pa_gates_payload(result, stand_in_gates=[], mode="fake")
    assert json.loads(json.dumps(payload)) == payload


def test_the_payload_lists_the_recoded_and_the_exactly_copied_spans(result, locus):
    pa = pa_gates_payload(result, stand_in_gates=[], mode="fake")["pa"]
    assert pa["locked_spans"] == [list(s) for s in result.locked_spans]
    assert pa["recode_spans"] == [list(s) for s in result.recode_spans]


# ---- the driver ------------------------------------------------------------------------------------------


def _inputs(loc, **over):
    base = dict(cargo=CARGO, genome=loc.genome, site_s=S, cds_mask=list(loc.cds_mask))
    base.update(over)
    return base


def test_the_fake_mode_driver_designs_and_returns_every_rendering(locus):
    out = design_pa(**_inputs(locus), mode="fake")
    assert out.report_text.startswith("=== Prime Assembly design: two pegRNAs and one donor ===")
    assert out.result.donor in out.report_text and out.csv_text and out.fasta_text.count(">") == 3
    assert (out.report_path, out.csv_path, out.fasta_path) == (
        None,
        None,
        None,
    ), "no output_dir, no files"


def test_the_driver_writes_the_five_files_and_they_hold_what_it_returns(locus, tmp_path):
    where = tmp_path / "nested" / "run"
    out = design_pa(**_inputs(locus), mode="fake", output_dir=where)
    assert sorted(p.name for p in where.iterdir()) == sorted(PA_FILES.values())
    assert (where / PA_FILES["report"]).read_text() == out.report_text
    # the bytes, not read_text(): the csv module ends its rows with CRLF, which reading as text would translate
    assert (where / PA_FILES["components"]).read_bytes().decode() == out.csv_text
    assert (where / PA_FILES["fasta"]).read_bytes().decode() == out.fasta_text
    assert (out.report_path, out.csv_path, out.fasta_path) == (
        str(where / PA_FILES["report"]),
        str(where / PA_FILES["components"]),
        str(where / PA_FILES["fasta"]),
    )
    rows = list(csv.DictReader(io.StringIO((where / PA_FILES["changes"]).read_text())))
    assert len(rows) == out.result.recode_edits
    assert {r["scope"] for r in rows} <= {"genomic arm"} and all(
        r["old_base"] != r["new_base"] for r in rows
    )


def test_the_gates_json_is_what_the_web_reads(locus, tmp_path):
    out = design_pa(**_inputs(locus), mode="fake", output_dir=tmp_path)
    gates = json.loads((tmp_path / PA_FILES["gates"]).read_text())
    assert gates == pa_gates_payload(out.result, stand_in_gates=out.stand_in_gates, mode="fake")
    assert gates["strategy"] == "pa" and gates["scorer_mode"] == "fake"
    assert [k for k, _ in gates["gate_results"]][-2:] == ["docking_a", "docking_b"]


def test_a_fake_run_declares_the_gates_a_pa_run_uses_and_not_pridict(locus):
    """The fake bundle also lists ``pridict`` (tj-PE's auto-spacer source), which PA never calls."""
    out = design_pa(**_inputs(locus), mode="fake")
    assert set(out.stand_in_gates) == {
        "splice",
        "structure",
        "importance",
        "docking",
        "offtarget",
        "cryptic",
    }
    assert "pridict" not in out.report_text.split("Provenance manifest:")[1]


def test_pridict_is_dropped_from_the_stand_in_gates_and_nothing_else():
    class Bundle:
        stub_gates = ("splice", "pridict", "docking")

    assert pa_stand_in_gates(Bundle()) == ("splice", "docking")


def test_in_api_mode_the_staged_resources_reach_the_scorer_factory(locus, monkeypatch, tmp_path):
    seen = {}

    def factory(*, mode, **kw):
        seen["mode"], seen["kw"] = mode, kw
        return GateScorers(**fake_scorers(), stub_gates=("docking",))

    monkeypatch.setattr(driver, "make_scorers", factory)
    out = design_pa(
        **_inputs(locus),
        mode="api",
        uniprot_acc="P49711",
        genome_fasta="/data/genome.fa",
        boltz_cache="/cache/boltz",
        max_mismatches=2,
        output_dir=tmp_path,
    )
    gates = json.loads((tmp_path / PA_FILES["gates"]).read_text())
    assert gates["scorer_mode"] == "api" and gates["stand_in_gates"] == ["docking"]
    assert seen["mode"] == "api"
    assert seen["kw"]["uniprot_acc"] == "P49711" and seen["kw"]["genome"] == "/data/genome.fa"
    assert seen["kw"]["boltz_cache"] == "/cache/boltz" and seen["kw"]["max_mismatches"] == 2
    assert out.stand_in_gates == ("docking",)
    banner = out.report_text.split("WHAT TO ORDER")[0]
    assert "NOT A VALIDATED DESIGN" in banner and "did NOT run a real model: docking" in banner


def test_in_fake_mode_the_scorer_factory_gets_no_resource_arguments(locus, monkeypatch):
    seen = {}

    def factory(*, mode, **kw):
        seen["kw"] = kw
        return GateScorers(**fake_scorers(), stub_gates=())

    monkeypatch.setattr(driver, "make_scorers", factory)
    out = design_pa(
        **_inputs(locus), mode="fake", genome_fasta="/data/genome.fa", uniprot_acc="P49711"
    )
    assert seen["kw"] == {}
    assert (
        "NOT A VALIDATED DESIGN" not in out.report_text
    ), "no stand-ins were declared by this bundle"


def test_a_real_bundle_without_stand_ins_gives_a_report_without_a_banner(locus, monkeypatch):
    monkeypatch.setattr(
        driver, "make_scorers", lambda **kw: GateScorers(**fake_scorers(), stub_gates=())
    )
    out = design_pa(**_inputs(locus), mode="api")
    assert out.stand_in_gates == () and "NOT A VALIDATED DESIGN" not in out.report_text


def test_the_transcript_and_the_window_are_handed_to_the_splice_context(locus, monkeypatch):
    calls = []

    def attach(scorer, *, transcript, genome):
        calls.append((transcript, genome))
        return "SpliceAI cross-check: real flanks (test)"

    monkeypatch.setattr(driver, "_attach_splice_context", attach)
    out = design_pa(**_inputs(locus), mode="fake", transcript="ENST00000264010")
    assert calls == [("ENST00000264010", locus.genome)]
    assert "splice_context: SpliceAI cross-check: real flanks (test)" in out.report_text


def test_the_callers_manifest_is_reported_and_never_modified(locus):
    manifest = {"gene": "CTCF", "uniprot": "P49711"}
    out = design_pa(**_inputs(locus), mode="fake", manifest=manifest)
    tail = out.report_text.split("Provenance manifest:")[1]
    assert "gene: CTCF" in tail and "uniprot: P49711" in tail and "scorer_mode: fake" in tail
    assert manifest == {"gene": "CTCF", "uniprot": "P49711"}


def test_the_progress_callback_receives_each_gate_with_its_numbers(locus):
    events = []
    design_pa(
        **_inputs(locus),
        mode="fake",
        progress_cb=lambda cand, gate, info, passed, metrics=None: events.append(
            (gate, passed, metrics)
        ),
    )
    gates = [g for g, _, _ in events]
    assert {
        "importance",
        "cryptic",
        "splice",
        "structure_a",
        "structure_b",
        "docking_a",
        "docking_b",
    } <= set(gates)
    assert gates.index("structure_a") < gates.index("docking_a") < gates.index("docking_b")
    assert all(passed for _, passed, _ in events)


def test_the_locks_and_parts_reach_the_walk(locus):
    out = design_pa(**_inputs(locus), mode="fake", cargo_locked=[[1, 12]])
    assert any("excluded by you" in line for line in out.result.changes.kept_exactly)


def test_every_walk_setting_is_forwarded_to_the_walk_unchanged(locus, monkeypatch):
    seen = {}
    real = driver.run_recoded_pa

    def spy(**kw):
        seen.update(kw)
        return real(**kw)

    monkeypatch.setattr(driver, "run_recoded_pa", spy)
    parts = [{"name": "all", "dna": CARGO, "locked": False}]
    settings = dict(
        pbs_lengths=(13, 17),
        flap_length=34,
        flap_min=31,
        flap_max=49,
        min_span=31,
        max_span=99,
        max_run=7,
        tail="UU",
        cas9="nCas9-test",
        max_offtargets=2,
        max_pairs=5,
        cargo_locked=[[1, 6]],
        cargo_parts=parts,
        protein_length=700,
        structure_threshold=0.9,
    )
    design_pa(**_inputs(locus), mode="fake", **settings)
    for key, value in settings.items():
        assert seen[key] == value, key
    assert (seen["cargo"], seen["genome"], seen["site_s"]) == (CARGO, locus.genome, S)


def test_the_driver_accepts_only_spcas9_ngg(locus):
    with pytest.raises(ValueError, match="NGG"):
        design_pa(**_inputs(locus), mode="fake", pam="NG")


def test_bad_arguments_are_refused_before_anything_runs(locus):
    with pytest.raises(ValueError, match="flap_length"):
        design_pa(**_inputs(locus), mode="fake", flap_length=20)


def test_a_refused_design_raises_the_walks_error_and_writes_nothing(locus, tmp_path):
    with pytest.raises(NoValidPADesign, match="design-recoded"):
        design_pa(**_inputs(locus, cargo=CARGO[:60]), mode="fake", output_dir=tmp_path / "o")
    assert not (tmp_path / "o").exists()


# ---- the command line -------------------------------------------------------------------------------------


def _config(loc, **over):
    cfg = dict(
        gene="TEST",
        uniprot_acc="P49711",
        cargo=CARGO,
        genome=loc.genome,
        site_s=S,
        cds_mask=list(loc.cds_mask),
    )
    cfg.update(over)
    return cfg


def _invoke(tmp_path, cfg, *extra, mode="fake"):
    path = tmp_path / "pa.json"
    path.write_text(json.dumps(cfg))
    out = tmp_path / "out"
    result = CliRunner().invoke(
        cli, ["design-pa", "--config", str(path), "--mode", mode, "--output-dir", str(out), *extra]
    )
    return result, out


def test_the_command_is_listed():
    result = CliRunner().invoke(cli, ["--help"])
    assert "design-pa" in result.output and "design-recoded" in result.output


def test_run_mode_pa_points_to_design_pa_and_submits_nothing(monkeypatch):
    """The cluster job behind ``run --mode pa`` was the prototype that could never pair real PRIDICT output."""
    from pegrna_pipeline.agent import launcher

    def submitted(*_a, **_k):
        raise AssertionError("a Prime Assembly run must not be queued on the cluster")

    monkeypatch.setattr(launcher.HPCLauncher, "launch_run", submitted)
    monkeypatch.setattr(launcher.HPCLauncher, "__init__", submitted)
    result = CliRunner().invoke(
        cli, ["run", "--gene", "CTCF", "--construct", "FLAG", "--mode", "pa", "--non-interactive"]
    )
    assert result.exit_code == 2
    assert "design-pa" in result.output and "examples/ctcf_cterm/config_pa.json" in result.output
    assert (
        "prototype" not in result.output.split("design-pa")[0]
    ), "it must not say it is merely a prototype"


def test_the_command_designs_prints_the_report_and_writes_the_files(locus, tmp_path):
    result, out = _invoke(tmp_path, _config(locus))
    assert result.exit_code == 0, result.output
    assert "=== Prime Assembly design" in result.output
    assert f"Report written to {out / PA_FILES['report']}" in result.output
    assert sorted(p.name for p in out.iterdir()) == sorted(PA_FILES.values())


def test_config_keys_that_are_not_parameters_become_provenance(locus, tmp_path):
    result, out = _invoke(tmp_path, _config(locus, _comment="a note", strategy="pa"))
    assert result.exit_code == 0, result.output
    tail = (out / PA_FILES["report"]).read_text().split("Provenance manifest:")[1]
    for line in (
        "comment: a note",
        "gene: TEST",
        "strategy: pa",
        "uniprot: P49711",
        "scorer_mode: fake",
    ):
        assert line in tail, line


def test_the_command_line_decides_mode_output_and_progress_not_the_config(locus, tmp_path):
    cfg = _config(
        locus, mode="api", output_dir="/nonexistent/elsewhere", progress_cb=1, manifest={"k": "v"}
    )
    result, out = _invoke(tmp_path, cfg)
    assert result.exit_code == 0, result.output
    assert json.loads((out / PA_FILES["gates"]).read_text())["scorer_mode"] == "fake"
    assert "k: v" in (out / PA_FILES["report"]).read_text()
    assert "@@P" not in result.output


def test_progress_prints_one_json_line_per_gate_and_then_the_report(locus, tmp_path):
    result, _ = _invoke(tmp_path, _config(locus), "--progress")
    assert result.exit_code == 0, result.output
    events = [json.loads(ln[3:]) for ln in result.output.splitlines() if ln.startswith("@@P ")]
    assert events and all(set(e) == {"cand", "gate", "info", "passed", "metrics"} for e in events)
    assert {"splice", "structure_a", "structure_b", "docking_a", "docking_b"} <= {
        e["gate"] for e in events
    }
    assert result.output.index("@@P ") < result.output.index("=== Prime Assembly design")


def test_a_refused_design_exits_nonzero_says_why_and_records_the_stand_ins(locus, tmp_path):
    result, out = _invoke(tmp_path, _config(locus, cargo=CARGO[:60]))
    assert result.exit_code == 1
    assert (
        "design-pa failed: NoValidPADesign:" in result.output and "design-recoded" in result.output
    )
    gates = json.loads((out / PA_FILES["gates"]).read_text())
    assert gates["failed"] is True and gates["strategy"] == "pa" and gates["scorer_mode"] == "fake"
    assert gates["error"].startswith("NoValidPADesign: ") and gates["gate_results"] == []
    assert "docking" in gates["stand_in_gates"] and "pridict" not in gates["stand_in_gates"]
    assert gates["diagnostics"] and "cargo 60 nt" in gates["diagnostics"][0]
    assert not (out / PA_FILES["report"]).exists(), "a refused design has no report to order from"


def test_a_walk_that_finds_no_pair_records_every_diagnostic(locus, tmp_path):
    result, out = _invoke(tmp_path, _config(locus, min_span=95, max_span=100))
    assert result.exit_code == 1
    gates = json.loads((out / PA_FILES["gates"]).read_text())
    assert any("PAM-in" in d for d in gates["diagnostics"]) and "NoValidPADesign" in gates["error"]


def test_a_bad_argument_fails_loudly_through_the_command(locus, tmp_path):
    result, out = _invoke(tmp_path, _config(locus, flap_length=12))
    assert result.exit_code == 1 and "design-pa failed: ValueError:" in result.output
    assert json.loads((out / PA_FILES["gates"]).read_text())["failed"] is True


def test_the_failure_record_never_masks_the_failure(locus, tmp_path, monkeypatch):
    import pegrna_pipeline.hpc.scorer_factory as factory

    def explode(**_kw):
        raise RuntimeError("probe failed")

    monkeypatch.setattr(factory, "make_scorers", explode)
    result, out = _invoke(tmp_path, _config(locus, cargo=CARGO[:60]))
    assert result.exit_code == 1 and "design-pa failed: NoValidPADesign" in result.output
    assert not (out / PA_FILES["gates"]).exists()


def test_in_api_mode_a_failure_probes_the_scorers_the_driver_would_have_used(
    locus, tmp_path, monkeypatch
):
    import pegrna_pipeline.hpc.scorer_factory as factory

    probes = []

    def make(*, mode, **kw):
        probes.append((mode, kw))
        return GateScorers(**fake_scorers(), stub_gates=("docking", "pridict"))

    monkeypatch.setattr(factory, "make_scorers", make)
    monkeypatch.setattr(driver, "make_scorers", make)
    monkeypatch.setenv("PEGRNA_GENOME", "/data/genome.fa")
    result, out = _invoke(tmp_path, _config(locus, cargo=CARGO[:60]), mode="api")
    assert result.exit_code == 1
    assert probes[-1] == ("api", {"uniprot_acc": "P49711", "genome": "/data/genome.fa"})
    gates = json.loads((out / PA_FILES["gates"]).read_text())
    assert gates["stand_in_gates"] == ["docking"], "the PA list: the bundle's, without PRIDICT"
    assert gates["scorer_mode"] == "api"


def test_design_recoded_still_writes_its_own_failure_record_through_the_shared_helper(tmp_path):
    """The failure record is shared code now: tj-PE keeps its file name, keys and stand-in list."""
    cfg = dict(
        rtt_aa="L",
        cargo="GACTACAAGGACGACGATGACAAG",
        genome="ACGT" * 20,
        site_s=40,
        cds_mask=[True] * 40 + [False] * 40,
        pbs_lengths=[13],
    )
    path = tmp_path / "tj.json"
    path.write_text(json.dumps(cfg))
    out = tmp_path / "tj_out"
    result = CliRunner().invoke(
        cli, ["design-recoded", "--config", str(path), "--mode", "fake", "--output-dir", str(out)]
    )
    assert result.exit_code == 1
    gates = json.loads((out / "recoded_tjpe_gates.json").read_text())
    assert set(gates) == {"gate_results", "failed", "error", "stand_in_gates", "scorer_mode"}
    assert "pridict" in gates["stand_in_gates"], "tj-PE does use PRIDICT, so it stays listed there"
    assert not (out / PA_FILES["gates"]).exists()


# ---- the shipped example: the user's CTCF cassette ---------------------------------------------------------

from pathlib import Path  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "ctcf_cterm" / "config_pa.json"


def test_the_ctcf_cassette_example_designs_in_fake_mode(tmp_path):
    cfg = json.loads(EXAMPLE.read_text())
    path = tmp_path / "example.json"
    path.write_text(json.dumps(cfg))
    out = tmp_path / "out"
    result = CliRunner().invoke(
        cli, ["design-pa", "--config", str(path), "--mode", "fake", "--output-dir", str(out)]
    )
    assert result.exit_code == 0, result.output
    gates = json.loads((out / PA_FILES["gates"]).read_text())
    pa = gates["pa"]
    assert len(cfg["cargo"]) == 499 and cfg["cargo"].endswith(
        "TAATACGACTCACTATAGG"
    ), "the two-G T7 cassette"
    assert pa["cargo_length"] == 499 and gates["deletion_size"] == 0
    assert [g["pegRNA_length"] <= 166 for g in pa["guides"]] == [True, True]
    assert pa["product_length"] == len(cfg["genome"]) + len(cfg["cargo"])
    donor = pa["donor"]
    assert cfg["cargo"] in donor["sequence"] and donor["sequence"] == "".join(
        part["sequence"] for part in donor["parts"]
    )
    cargo_part = next(part for part in donor["parts"] if part["name"] == "cargo")
    assert cargo_part["sequence"] == cfg["cargo"], "the cassette is used exactly as given"
    assert gates["changes"]["cargo_unchanged"] is True
    assert any("T7 promoter" in line for line in gates["changes"]["kept_exactly"])
