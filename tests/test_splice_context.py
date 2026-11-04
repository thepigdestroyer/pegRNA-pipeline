"""SpliceAI binds only on real genomic context; without it the number is advisory (ADR 0028).

A real run of the CTCF example with the 501-bp GGGGS-FLAG-2A-resistance cassette passed the MaxEntScan gate and
was then rejected on every candidate by the SpliceAI cross-check: "splice delta 0.336 >= 0.2". ADR 0011 had made
SpliceAI an ADVISORY cross-check that binds only when a caller supplies real >= 10 kb genomic context. In api mode
the 297-nt window is padded with N to reach 10 kb, so a number came back where the original code raised, and the
walk applied it as a hard gate.

Measured on the CTCF locus (scripts/analysis/probe_splice_gate.py), padded with N, all of 10 random synonymous
recodings of the cassette and all of 10 random 501-nt sequences scored >= 0.2 (0.286 to 0.774), and this cassette
scored 0.336. With the real +/-5 kb flanks the same cassette scored 0.159, and the strongest site SpliceAI saw
inside the cassette fell from 0.321 to 0.007. So the padded number was a length effect, not evidence, and real
context is informative (5 of 6 random recodings still failed it).

The pure logic is tested here with stand-ins; the numbers above are from the real models.
"""

import dataclasses
import os
import random

import pytest
import test_pridict_ngrna_wiring as W

from pegrna_pipeline.hpc import recoded_driver
from pegrna_pipeline.hpc.genome_context import MissingGenomicContextError
from pegrna_pipeline.hpc.genome_ingest import Flanks, bounded_ensembl, fetch_flanking_sequence
from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA, run_recoded_tjpe
from pegrna_pipeline.hpc.scorer_factory import _PaddedSpliceAIScorer
from pegrna_pipeline.hpc.splice_screen import SpliceFake

HIGH, LOW = 0.34, 0.15
NO_CONTEXT = "no genomic context: the design window padded with N"


class _Splice:
    """A SpliceAI stand-in with a fixed number and a stated context."""

    def __init__(self, delta, *, binding, note=NO_CONTEXT):
        self._d, self.binding, self.context_note = delta, binding, note

    def delta(self, ref_seq, alt_seq):  # noqa: ARG002
        return self._d


def _walk(scorer, **over):
    return run_recoded_tjpe(**W._wf_inputs(splice_scorer=scorer, **over))


# ---- the walk ------------------------------------------------------------------------------------


def test_a_high_number_from_a_scorer_without_real_context_does_not_reject():
    res = _walk(_Splice(HIGH, binding=False))
    info = dict(res.gate_results)["splice"]
    assert "advisory only" in info and "0.340" in info and "padded with N" in info


def test_the_same_number_from_a_scorer_with_real_context_rejects():
    with pytest.raises(NoValidPegRNA) as exc:
        _walk(_Splice(HIGH, binding=True, note="real genomic context (Ensembl)"))
    assert "splice delta 0.340 >= 0.2" in " ".join(exc.value.diagnostics)


def test_a_scorer_that_does_not_say_it_is_advisory_binds_as_it_always_did():
    with pytest.raises(NoValidPegRNA) as exc:
        _walk(SpliceFake(delta=HIGH))
    assert "splice delta 0.340 >= 0.2" in " ".join(exc.value.diagnostics)


def test_the_advisory_is_a_warning_and_says_what_the_primary_gate_found():
    res = _walk(_Splice(HIGH, binding=False))
    (warning,) = [w for w in res.warnings if w.startswith("splice cross-check:")]
    assert "does not reject" in warning and "MaxEntScan) passed" in warning
    assert "cannot tell a clean large insert from a bad one" in warning


def test_a_number_under_the_cutoff_is_reported_with_its_context_and_raises_no_warning():
    res = _walk(_Splice(LOW, binding=False))
    assert "splice delta 0.150 ok (advisory only: " + NO_CONTEXT in dict(res.gate_results)["splice"]
    assert not any(w.startswith("splice cross-check:") for w in res.warnings)


def test_a_binding_verdict_names_the_context_it_was_computed_in():
    res = _walk(_Splice(LOW, binding=True, note="real genomic context (Ensembl GRCh38 chr16)"))
    info = dict(res.gate_results)["splice"]
    assert info.startswith("splice delta 0.150 ok") and "real genomic context" in info


def test_the_advisory_reaches_the_progress_stream_as_a_pass_with_no_cutoff_to_contradict_it():
    """The web flowchart draws a cutoff comparison when a gate event carries one; an advisory verdict that
    claimed 'outside cutoff' under a green tick would contradict itself."""
    events = []
    _walk(
        _Splice(HIGH, binding=False),
        on_gate=lambda cand, gate, info, passed, metrics=None: events.append(
            (gate, passed, metrics or {})
        ),
    )
    gate, passed, metrics = next(e for e in events if e[0] == "splice")
    assert passed is True
    assert "cutoff" not in metrics and metrics["value"] == HIGH
    assert metrics["note"].startswith("SpliceAI Δ 0.34: advisory, not binding")


def test_a_context_free_number_cannot_stop_a_design_the_primary_gate_passed():
    """The regression itself: every candidate used to fail with the same number."""
    res = _walk(_Splice(HIGH, binding=False), pbs_lengths=(13, 14, 15, 17, 19))
    assert res.pegRNA_sequence and res.candidates_tried >= 1


# ---- the scorer ------------------------------------------------------------------------------------


class _Inner:
    _context = 5000

    def __init__(self):
        self.calls = []

    def delta(self, ref_seq, alt_seq):
        self.calls.append((ref_seq, alt_seq))
        return 0.1


def _scorer():
    s = _PaddedSpliceAIScorer()
    s._inner = _Inner()
    return s


REF = "ACGT" * 20
ALT = REF + "TTGCA" * 100
LEFT, RIGHT = "AC" * 2500, "GT" * 2500


def test_without_context_both_alleles_are_padded_with_n_and_it_does_not_bind():
    s = _scorer()
    s.delta(REF, ALT)
    ref, alt = s._inner.calls[0]
    assert ref == "N" * 5000 + REF + "N" * 5000 and alt == "N" * 5000 + ALT + "N" * 5000
    assert s.binding is False and "no genomic context" in s.context_note


def test_with_complete_real_flanks_they_surround_both_alleles_and_it_binds():
    s = _scorer()
    s.set_context(LEFT, RIGHT, complete=True, source="Ensembl X")
    s.delta(REF, ALT)
    ref, alt = s._inner.calls[0]
    assert ref == LEFT + REF + RIGHT and alt == LEFT + ALT + RIGHT
    assert s.binding is True and "real genomic context (Ensembl X)" == s.context_note


def test_flanks_that_run_out_are_used_but_do_not_bind():
    s = _scorer()
    s.set_context("N" * 3000 + LEFT[:2000], RIGHT, complete=False, source="Ensembl X")
    s.delta(REF, ALT)
    assert s.binding is False and "partial genomic context" in s.context_note


def test_a_window_the_caller_already_made_ten_kb_long_is_scored_as_given_and_binds():
    s = _scorer()
    big_ref, big_alt = "ACGT" * 2600, "ACGT" * 2600 + "TTGCA" * 100
    s.delta(big_ref, big_alt)
    assert s._inner.calls == [(big_ref, big_alt)]
    assert s.binding is True and "supplied by the caller" in s.context_note


def test_a_flank_of_the_wrong_length_is_refused():
    with pytest.raises(ValueError, match="exactly 5000"):
        _scorer().set_context("A" * 10, RIGHT, complete=True, source="x")


def test_the_flanks_are_remembered_across_calls():
    s = _scorer()
    s.set_context(LEFT, RIGHT, complete=True, source="x")
    s.delta(REF, ALT)
    s.delta(REF, REF + "GATTACA" * 50)
    assert all(r.startswith(LEFT) and r.endswith(RIGHT) for pair in s._inner.calls for r in pair)


# ---- fetching the flanks ---------------------------------------------------------------------------


class _Ensembl:
    """Stands in for the REST API: one transcript, one region whose content the test controls."""

    def __init__(self, region, *, strand=1, start=100_000, end=190_000, chrom="16"):
        self.region, self.urls = region, []
        self.lookup = {"seq_region_name": chrom, "start": start, "end": end, "strand": strand}

    def __call__(self, url, **params):
        self.urls.append(url)
        if "/lookup/id/" in url:
            return self.lookup
        return {"seq": self.region}


def _dna(n, seed):
    rng = random.Random(seed)  # one generator for the whole sequence, not one per base
    return "".join(rng.choice("ACGT") for _ in range(n))


def _region(left_len, window, right_len, seed=1):
    return _dna(left_len, seed) + window + _dna(right_len, seed + 1000)


WINDOW = _dna(120, 9)


def test_the_flanks_are_the_real_sequence_either_side_of_the_one_exact_match():
    region = _region(7000, WINDOW, 7000)
    flanks = fetch_flanking_sequence("ENST00000264010", WINDOW, http_get=_Ensembl(region))
    assert flanks.left == region[2000:7000] and flanks.right == region[7120:12120]
    assert flanks.complete and "ENST00000264010" in flanks.source and "(+)" in flanks.source


def test_the_request_asks_for_the_transcripts_own_strand_with_room_for_the_flanks():
    api = _Ensembl(_region(7000, WINDOW, 7000), strand=-1, start=100_000, end=190_000)
    fetch_flanking_sequence("ENST1", WINDOW, http_get=api)
    assert api.urls[-1].endswith("/sequence/region/human/16:94000..196000:-1")


def test_a_gene_near_the_start_of_the_chromosome_does_not_ask_for_a_negative_coordinate():
    api = _Ensembl(_region(7000, WINDOW, 7000), start=2000, end=90_000)
    fetch_flanking_sequence("ENST1", WINDOW, http_get=api)
    assert "/16:1..96000:1" in api.urls[-1]


def test_where_the_genome_runs_out_the_flank_is_n_filled_and_flagged_incomplete():
    region = _region(1200, WINDOW, 7000)
    flanks = fetch_flanking_sequence("ENST1", WINDOW, http_get=_Ensembl(region))
    assert flanks.left == "N" * 3800 + region[:1200] and len(flanks.left) == 5000
    assert flanks.right == region[1320:6320] and not flanks.complete


def test_a_window_that_is_not_in_the_transcript_is_refused_not_guessed():
    other = _region(7000, _dna(120, 5), 7000)
    with pytest.raises(MissingGenomicContextError, match="not in the fetched region"):
        fetch_flanking_sequence("ENST1", WINDOW, http_get=_Ensembl(other))


def test_a_window_that_occurs_twice_has_ambiguous_flanks_and_is_refused():
    region = _region(3000, WINDOW, 3000) + _region(1000, WINDOW, 3000, seed=2)
    with pytest.raises(MissingGenomicContextError, match="more than once"):
        fetch_flanking_sequence("ENST1", WINDOW, http_get=_Ensembl(region))


def test_a_lookup_without_coordinates_is_refused():
    api = _Ensembl("ACGT")
    api.lookup = {"error": "ID not found"}
    with pytest.raises(MissingGenomicContextError, match="no usable coordinates"):
        fetch_flanking_sequence("ENST1", WINDOW, http_get=api)


def test_a_bounded_budget_is_set_for_the_call_and_restored(monkeypatch):
    monkeypatch.delenv("PEGRNA_ENSEMBL_RETRIES", raising=False)
    monkeypatch.delenv("PEGRNA_ENSEMBL_TIMEOUT", raising=False)
    with bounded_ensembl(2, 45):
        assert os.environ["PEGRNA_ENSEMBL_RETRIES"] == "2"
        assert os.environ["PEGRNA_ENSEMBL_TIMEOUT"] == "45"
    assert "PEGRNA_ENSEMBL_RETRIES" not in os.environ and "PEGRNA_ENSEMBL_TIMEOUT" not in os.environ


def test_an_operators_own_budget_is_left_alone(monkeypatch):
    monkeypatch.setenv("PEGRNA_ENSEMBL_RETRIES", "7")
    with bounded_ensembl(2, 45):
        assert os.environ["PEGRNA_ENSEMBL_RETRIES"] == "7"
    assert os.environ["PEGRNA_ENSEMBL_RETRIES"] == "7"


# ---- the driver ------------------------------------------------------------------------------------


class _ContextScorer:
    """A splice stand-in that takes flanks like the api scorer, and binds only once it has complete ones."""

    def __init__(self, delta):
        self._d, self.flanks = delta, None

    def set_context(self, left, right, *, complete, source):
        self.flanks = (left, right, complete, source)

    @property
    def binding(self):
        return bool(self.flanks and self.flanks[2])

    @property
    def context_note(self):
        return f"real genomic context ({self.flanks[3]})" if self.binding else NO_CONTEXT

    def delta(self, ref_seq, alt_seq):  # noqa: ARG002
        return self._d


def _fetched(left="L", right="R", complete=True, source="Ensembl X"):
    return lambda transcript, window: Flanks(left * 5000, right * 5000, complete, source)


def _failing(transcript, window):
    raise RuntimeError("Ensembl request failed after 2 attempts")


def test_the_driver_attaches_real_flanks_to_a_scorer_that_takes_them(monkeypatch):
    monkeypatch.setattr(recoded_driver, "fetch_flanking_sequence", _fetched())
    scorer = _ContextScorer(LOW)
    note = recoded_driver._attach_splice_context(scorer, transcript="ENST1", genome="ACGT")
    assert scorer.flanks == ("L" * 5000, "R" * 5000, True, "Ensembl X")
    assert note == "SpliceAI cross-check: Ensembl X"


def test_a_partial_context_is_attached_and_labelled_partial(monkeypatch):
    monkeypatch.setattr(recoded_driver, "fetch_flanking_sequence", _fetched(complete=False))
    note = recoded_driver._attach_splice_context(
        _ContextScorer(LOW), transcript="ENST1", genome="A"
    )
    assert note.endswith("(partial)")


def test_a_failed_fetch_never_fails_the_design_and_says_so(monkeypatch):
    monkeypatch.setattr(recoded_driver, "fetch_flanking_sequence", _failing)
    scorer = _ContextScorer(LOW)
    note = recoded_driver._attach_splice_context(scorer, transcript="ENST1", genome="ACGT")
    assert scorer.flanks is None
    assert "unavailable" in note and "RuntimeError" in note and "advisory only" in note


def test_no_transcript_means_no_fetch(monkeypatch):
    monkeypatch.setattr(
        recoded_driver, "fetch_flanking_sequence", _failing
    )  # would raise if called
    scorer = _ContextScorer(LOW)
    note = recoded_driver._attach_splice_context(scorer, transcript=None, genome="ACGT")
    assert scorer.flanks is None and "no transcript given" in note


def test_a_scorer_that_takes_no_context_is_left_alone_and_nothing_is_fetched(monkeypatch):
    monkeypatch.setattr(recoded_driver, "fetch_flanking_sequence", _failing)
    assert (
        recoded_driver._attach_splice_context(SpliceFake(0.0), transcript="ENST1", genome="A") == ""
    )


def _design(monkeypatch, scorer, fetch, **extra):
    real = recoded_driver.make_scorers
    monkeypatch.setattr(
        recoded_driver,
        "make_scorers",
        lambda **kw: dataclasses.replace(real(**kw), splice_scorer=scorer),
    )
    monkeypatch.setattr(recoded_driver, "fetch_flanking_sequence", fetch)
    return recoded_driver.design_recoded(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome=W.GENOME,
        site_s=W.SITE_S,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
        mode="fake",
        transcript="ENST1",
        **extra,
    )


def test_with_real_context_a_high_number_rejects_the_design(monkeypatch):
    with pytest.raises(NoValidPegRNA) as exc:
        _design(monkeypatch, _ContextScorer(HIGH), _fetched())
    assert "splice delta 0.340 >= 0.2" in " ".join(exc.value.diagnostics)


def test_without_real_context_the_same_number_is_a_warning_and_the_design_goes_through(monkeypatch):
    """The situation the user hit: the same number, but nothing real to judge it against."""
    out = _design(monkeypatch, _ContextScorer(HIGH), _failing)
    assert out.result.pegRNA_sequence
    assert any(w.startswith("splice cross-check:") for w in out.result.warnings)
    assert "advisory only" in dict(out.result.gate_results)["splice"]


def test_the_run_manifest_records_where_the_splice_context_came_from(monkeypatch):
    out = _design(monkeypatch, _ContextScorer(LOW), _fetched(source="Ensembl GRCh38 chr16:1-9 (+)"))
    assert "splice_context: SpliceAI cross-check: Ensembl GRCh38 chr16:1-9 (+)" in out.report_text
    assert (
        "real genomic context (Ensembl GRCh38 chr16:1-9 (+))"
        in dict(out.result.gate_results)["splice"]
    )
