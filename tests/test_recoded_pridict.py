"""Stage C (C1/C3/C4): two-pegRNA PRIDICT adapter + orientation + ngRNA lifecycle.

The tj-pegRNA is assembled from TWO PRIDICT designs (ADR 0014/0016):
- pegRNA#1 (designed on the target strand) -> PBS1 + first insert end, used UNTOUCHED
  (PRIDICT already emits them as RC(TS) via PBSrevcomp/RTrevcomp).
- pegRNA#2 (designed on the OPPOSITE strand) -> rcPBS2 = RC(its PBSrevcomp),
  last_end = RC(its RTrevcomp); its spacer+scaffold BECOMES the ngRNA, bound to the design.
So every 3'-extension component ends up RC(TS).
"""

from __future__ import annotations

import dataclasses

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import PegRNADesign
from pegrna_pipeline.hpc.recoded_pridict import (
    TwoPegDesign,
    TwoPegPridictFake,
    combine_two_peg,
    exclude_rtgdna_tail,
    split_rt_genomic_tail,
)

SCAF = "GTTTTAGAGCTAGAAATAGCAAG"  # stand-in scaffold for tests


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def _peg(name: str, spacer: str, pbs: str, rt: str, ts: str = "Rv") -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer=spacer,
        pbs=pbs,
        rt=rt,
        pbs_length=len(pbs),
        rt_length=len(rt),
        target_strand=ts,
        editing_position=5,
        pe2_hek293t=50.0,
        pe2_k562=50.0,
        pegrna_full="",
    )


# --------------------------------------------------------------- C1: dataclass + fake
def test_twopegdesign_offtarget_risk_defaults_none():
    d = TwoPegDesign(
        pbs1="AAA",
        first_end="CCC",
        rcpbs2="GGG",
        last_end="TTT",
        ngrna="ACGT" + SCAF,
        pbs1_len=3,
        rt_len=3,
        editing_position=5,
        target_strand="Rv",
    )
    # PRIDICT2 exposes no per-spacer off-target-risk field (ADR 0016) -> None.
    assert d.offtarget_risk is None


def test_fake_emits_pridict_orientation():
    # All window args are target-strand. Real PRIDICT emits RC(TS): pegRNA#1 (on TS) ->
    # RC(TS) windows; pegRNA#2 (on OS) -> PBSrevcomp=RC(OS)=TS literally. The fake mirrors that.
    peg1, peg2 = TwoPegPridictFake().design_pair(
        first_end="ACGTACGT",
        last_end="TTTTGGGG",
        pbs1_window="CCCCAAAA",
        pbs2_window="GGGGTTTT",
        spacer1="G" * 20,
        spacer2="A" * 20,
        target_strand="Rv",
    )
    assert isinstance(peg1, PegRNADesign) and isinstance(peg2, PegRNADesign)
    assert peg1.pbs == _rc("CCCCAAAA") and peg1.rt == _rc("ACGTACGT")  # peg1 on TS -> RC(TS)
    assert peg2.pbs == "GGGGTTTT" and peg2.rt == "TTTTGGGG"  # peg2 on OS -> TS literally


def test_fake_then_combine_is_rc_of_target_strand_end_to_end():
    # The CI path (fake -> combine) must yield RC(TS) for every component starting from
    # target-strand windows -- guards a future fake regression that silently inverts a strand.
    kw = dict(
        first_end="ACGTACGTAC",
        last_end="TTTTGGGGAC",
        pbs1_window="CCCCAAAAGT",
        pbs2_window="GGGGTTTTCA",
        spacer1="G" * 20,
        spacer2="A" * 20,
        target_strand="Rv",
    )
    peg1, peg2 = TwoPegPridictFake().design_pair(**kw)
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.pbs1 == _rc(kw["pbs1_window"])
    assert d.first_end == _rc(kw["first_end"])
    assert d.rcpbs2 == _rc(kw["pbs2_window"])
    assert d.last_end == _rc(kw["last_end"])
    assert d.rcpbs2 != kw["pbs2_window"]  # the RC is not a no-op (non-palindromic window)


# ---------------------------------------------------- C3: orientation = RC(TS) end-to-end
def test_combine_two_peg_is_rc_of_target_strand_end_to_end():
    ts_pbs1, ts_first, ts_pbs2, ts_last = "CCCCAAAAGT", "ACGTACGTAC", "GGGGTTTTCA", "TTTTGGGGAC"
    # pegRNA#1 designed on TS: PRIDICT emits RC(TS windows).
    peg1 = _peg("p1", "G" * 20, pbs=_rc(ts_pbs1), rt=_rc(ts_first))
    # pegRNA#2 designed on the OPPOSITE strand: its revcomp fields == the TS windows literally.
    peg2 = _peg("p2", "A" * 20, pbs=ts_pbs2, rt=ts_last)
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.pbs1 == _rc(ts_pbs1)  # untouched (already RC(TS))
    assert d.first_end == _rc(ts_first)  # untouched (already RC(TS))
    assert d.rcpbs2 == _rc(ts_pbs2)  # RC of peg2's PBSrevcomp -> RC(TS)
    assert d.last_end == _rc(ts_last)  # RC of peg2's RTrevcomp -> RC(TS)
    # Every extension component is RC of its target-strand window:
    for comp, ts_win in [
        (d.pbs1, ts_pbs1),
        (d.first_end, ts_first),
        (d.rcpbs2, ts_pbs2),
        (d.last_end, ts_last),
    ]:
        assert comp == _rc(ts_win)
    # Non-tautology guard: the RC on peg2 is actually applied (would fail if combine
    # dropped the RC, since the windows are non-palindromic).
    assert d.rcpbs2 != peg2.pbs
    assert d.last_end != peg2.rt


def test_combine_pbs1_and_first_end_are_verbatim_from_peg1():
    peg1 = _peg("p1", "G" * 20, pbs="ACGTACGT", rt="TTTTCCCCAA")
    peg2 = _peg("p2", "A" * 20, pbs="GGGGAAAA", rt="CCCCGGGG")
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.pbs1 == peg1.pbs and d.first_end == peg1.rt  # verbatim, not re-oriented
    assert d.pbs1_len == peg1.pbs_length and d.rt_len == peg1.rt_length


# ----------------------------------------------------------------- C4: ngRNA lifecycle
def test_ngrna_bound_as_peg2_spacer_plus_scaffold():
    peg1 = _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT")
    peg2 = _peg("p2", "ACGTACGTACGTACGTACGT", pbs="GGGG", rt="CCCC")
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.ngrna == peg2.spacer + SCAF  # ngRNA = pegRNA#2 spacer + scaffold, untouched


def test_ngrna_rebinds_on_new_peg2():
    # The ngRNA is bound to (never independent of) its tj-pegRNA: a new peg2 -> a new design
    # with a new ngRNA; discarding the design discards the ngRNA (they are one object).
    peg1 = _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT")
    da = combine_two_peg(peg1, _peg("p2a", "A" * 20, pbs="GGGG", rt="CCCC"), scaffold=SCAF)
    db = combine_two_peg(peg1, _peg("p2b", "C" * 20, pbs="GGGG", rt="CCCC"), scaffold=SCAF)
    assert da.ngrna != db.ngrna
    assert db.ngrna == "C" * 20 + SCAF


# ---------------------------------------------------------------- misc invariants
def test_combine_propagates_peg1_edit_metadata_and_offtarget_risk():
    peg1 = _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT", ts="Fw")
    peg2 = _peg("p2", "A" * 20, pbs="GGGG", rt="CCCC")
    d = combine_two_peg(peg1, peg2, scaffold=SCAF, offtarget_risk=0.7)
    assert d.editing_position == peg1.editing_position  # peg1 authoritative for edit metadata
    assert d.target_strand == peg1.target_strand == "Fw"
    assert d.offtarget_risk == 0.7  # passthrough when a caller supplies one


def test_twopegdesign_is_frozen():
    d = combine_two_peg(
        _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT"),
        _peg("p2", "A" * 20, pbs="GGGG", rt="CCCC"),
        scaffold=SCAF,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        d.rcpbs2 = "AAAA"  # type: ignore[misc]


# ---------------------------------------------------- C5: RTgDNA / RT-tail exclusion
def test_split_rt_genomic_tail_separates_insert_and_tail():
    insert_part = "ACGTACGTAC"
    tail = "GGGGCCCCAA"  # genomic-homology tail PRIDICT appended past the insert
    ins, got_tail = split_rt_genomic_tail(insert_part + tail, insert_part)
    assert ins == insert_part and got_tail == tail


def test_split_rt_genomic_tail_empty_when_no_tail():
    insert_part = "ACGTACGTAC"
    ins, tail = split_rt_genomic_tail(insert_part, insert_part)
    assert ins == insert_part and tail == ""


def test_split_rt_genomic_tail_raises_when_misaligned():
    # The RT must start with the intended insert part; otherwise the RT and insert are
    # misaligned -> fail loud rather than silently mis-trim.
    with pytest.raises(ValueError):
        split_rt_genomic_tail("TTTTACGTACGT", "ACGTACGTAC")


def test_split_rt_genomic_tail_rejects_empty_insert_part():
    # Regression: an empty intended insert starts every string and would silently mark the
    # whole RT as genomic tail, deleting the insert end instead of failing loud.
    with pytest.raises(ValueError, match="insert_part"):
        split_rt_genomic_tail("ACGTACGT", "")


def test_exclude_rtgdna_tail_drops_the_genomic_tail():
    # C5: the recoded between-nicks arm authoritatively covers the genomic span, so the
    # RT's genomic-homology tail is excluded -- only the insert part remains.
    insert_part = "ACGTACGTAC"
    rt = insert_part + "GGGGCCCCAA"
    assert exclude_rtgdna_tail(rt, insert_part) == insert_part


def test_exclude_rtgdna_tail_noop_when_no_tail():
    insert_part = "ACGTACGTAC"
    assert exclude_rtgdna_tail(insert_part, insert_part) == insert_part


def test_exclude_rtgdna_tail_keeps_non_overlapping_tail_after_arm_overlap():
    # C5 excludes only the portion of PRIDICT's genomic tail that overlaps the separately
    # recoded between-nicks arm. If the tail extends past that arm, the non-overlapping suffix
    # must stay so the RT template is not over-trimmed.
    insert_part = "ACGTACGTAC"
    overlap = "GGGG"
    outside_arm = "CCCC"
    rt = insert_part + overlap + outside_arm
    assert (
        exclude_rtgdna_tail(rt, insert_part, overlap_len=len(overlap)) == insert_part + outside_arm
    )


def test_exclude_rtgdna_tail_overlap_len_zero_keeps_full_tail():
    # overlap_len=0 is the explicit inverse of the default (None -> drop whole tail): keep it all.
    # Guards a "overlap_len = overlap_len or len(tail)" refactor that would treat 0 as falsy.
    insert_part = "ACGTACGTAC"
    tail = "GGGGCCCC"
    assert exclude_rtgdna_tail(insert_part + tail, insert_part, overlap_len=0) == insert_part + tail


def test_exclude_rtgdna_tail_overlap_len_equals_tail_len_drops_whole_tail():
    insert_part = "ACGTACGTAC"
    tail = "GGGGCCCC"
    assert (
        exclude_rtgdna_tail(insert_part + tail, insert_part, overlap_len=len(tail)) == insert_part
    )


def test_exclude_rtgdna_tail_rejects_out_of_range_overlap_len():
    insert_part = "ACGTACGTAC"
    rt = insert_part + "GGGG"  # tail length 4
    with pytest.raises(ValueError, match="overlap_len"):
        exclude_rtgdna_tail(rt, insert_part, overlap_len=-1)
    with pytest.raises(ValueError, match="overlap_len"):
        exclude_rtgdna_tail(rt, insert_part, overlap_len=5)


# ---------------------------------------------------------- C2: real adapter (mocked I/O)
class _FakeSSH:
    def __init__(self):
        self.downloads = []

    def download(self, remote, local):
        self.downloads.append((remote, local))

    def exec(self, *a, **k):
        return ("", "", 0)

    def upload(self, *a, **k):
        pass


def _pegp(name: str, pe2: float, *, k562: float | None = None) -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer="G" * 20,
        pbs="ACGT",
        rt="TTTT",
        pbs_length=4,
        rt_length=4,
        target_strand="Rv",
        editing_position=5,
        pe2_hek293t=pe2,
        pe2_k562=pe2 if k562 is None else k562,
        pegrna_full="",
    )


def test_real_adapter_submits_and_picks_top_pe2_per_end(monkeypatch):
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    monkeypatch.setattr(
        ps,
        "submit_pridict_batches",
        lambda **kw: BatchResult(submitted=[object()], failed=[], output_dirs=["/out"]),
    )
    loaded = {
        "e1": [_pegp("e1", 10.0), _pegp("e1", 90.0)],
        "e2": [_pegp("e2", 30.0), _pegp("e2", 70.0)],
    }
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    monkeypatch.setattr(adapter, "_load_designs", lambda output_dir, name: loaded[name])
    peg1, peg2 = adapter.design_pair(
        editseq1="ACGT(+GG)AC", editseq2="TTTT(+CC)AC", name1="e1", name2="e2"
    )
    assert peg1.pe2_hek293t == 90.0 and peg2.pe2_hek293t == 70.0  # top PE2 per end


def test_real_adapter_defaults_to_project_score_cell_line_k562(monkeypatch):
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    monkeypatch.setattr(
        ps,
        "submit_pridict_batches",
        lambda **kw: BatchResult(submitted=[object()], failed=[], output_dirs=["/out"]),
    )
    loaded = {
        # HEK would pick the first row; K562 must pick the second row by default.
        "e1": [_pegp("e1", 90.0, k562=10.0), _pegp("e1", 20.0, k562=99.0)],
        "e2": [_pegp("e2", 80.0, k562=15.0), _pegp("e2", 30.0, k562=88.0)],
    }
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    monkeypatch.setattr(adapter, "_load_designs", lambda output_dir, name: loaded[name])
    peg1, peg2 = adapter.design_pair(editseq1="x", editseq2="y", name1="e1", name2="e2")
    assert peg1.pe2_k562 == 99.0 and peg2.pe2_k562 == 88.0


def test_real_adapter_can_rank_by_hek293t_when_explicit(monkeypatch):
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    monkeypatch.setattr(
        ps,
        "submit_pridict_batches",
        lambda **kw: BatchResult(submitted=[object()], failed=[], output_dirs=["/out"]),
    )
    loaded = {
        "e1": [_pegp("e1", 91.0, k562=10.0), _pegp("e1", 20.0, k562=99.0)],
        "e2": [_pegp("e2", 81.0, k562=15.0), _pegp("e2", 30.0, k562=88.0)],
    }
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(),
        remote_run_dir="/run",
        pridict_install_path="/PRIDICT2",
        score_cell_line="hek293t",
    )
    monkeypatch.setattr(adapter, "_load_designs", lambda output_dir, name: loaded[name])
    peg1, peg2 = adapter.design_pair(editseq1="x", editseq2="y", name1="e1", name2="e2")
    assert peg1.pe2_hek293t == 91.0 and peg2.pe2_hek293t == 81.0


def test_real_adapter_rejects_invalid_score_cell_line_at_construction():
    # Fail fast: an invalid score column must raise when the adapter is built, NOT after a
    # multi-hour CX3 batch has already run (reviewers flagged the late-validation trap).
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    with pytest.raises(ValueError, match="score_cell_line"):
        TwoPegPridictReal(
            ssh=_FakeSSH(),
            remote_run_dir="/run",
            pridict_install_path="/PRIDICT2",
            score_cell_line="cho",
        )


def test_real_adapter_revalidates_score_cell_line_before_submitting(monkeypatch):
    # __post_init__ fails fast, but the dataclass is mutable: reassigning score_cell_line after
    # construction must not slip an invalid column past submit and only fail post-batch.
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    calls = {"n": 0}

    def _spy(**kw):
        calls["n"] += 1
        return BatchResult(submitted=[object()], failed=[], output_dirs=["/out"])

    monkeypatch.setattr(ps, "submit_pridict_batches", _spy)
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    adapter.score_cell_line = "cho"  # mutate past the __post_init__ guard
    with pytest.raises(ValueError, match="score_cell_line"):
        adapter.design_pair(editseq1="x", editseq2="y", name1="e1", name2="e2")
    assert calls["n"] == 0  # never submitted a batch with an invalid score column


def test_real_adapter_validates_names_before_submitting(monkeypatch):
    # design_pair must reject bad/duplicate names BEFORE the expensive submit runs.
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    calls = {"n": 0}

    def _spy(**kw):
        calls["n"] += 1
        return BatchResult(submitted=[object()], failed=[], output_dirs=["/out"])

    monkeypatch.setattr(ps, "submit_pridict_batches", _spy)
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    # Equal names would collide on one PRIDICT CSV -> identical pair.
    with pytest.raises(ValueError, match="must differ"):
        adapter.design_pair(editseq1="x", editseq2="y", name1="dup", name2="dup")
    # Path-traversal / shell-metachar names are rejected.
    with pytest.raises(ValueError, match=r"\[A-Za-z0-9_-\]"):
        adapter.design_pair(editseq1="x", editseq2="y", name1="../evil", name2="e2")
    assert calls["n"] == 0  # never submitted a batch for an invalid request


def test_real_adapter_rejects_submit_kwargs_that_override_fixed_fields():
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    with pytest.raises(ValueError, match="submit_kwargs.*remote_run_dir"):
        adapter.design_pair(
            editseq1="x",
            editseq2="y",
            name1="e1",
            name2="e2",
            remote_run_dir="/clobber",
        )
    # `variants` is the most dangerous collision (it would silently retarget the whole batch)
    # and is a distinct guarded field, so cover it too.
    with pytest.raises(ValueError, match="submit_kwargs.*variants"):
        adapter.design_pair(
            editseq1="x",
            editseq2="y",
            name1="e1",
            name2="e2",
            variants=[{"sequence_name": "z", "editseq": "z"}],
        )


def test_real_adapter_forwards_extra_submit_kwargs(monkeypatch):
    # The exclusion guard must NOT block legitimate passthrough kwargs (e.g. a walltime override):
    # they should reach submit_pridict_batches unchanged.
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    seen = {}

    def _spy(**kw):
        seen.update(kw)
        return BatchResult(submitted=[object()], failed=[], output_dirs=["/out"])

    monkeypatch.setattr(ps, "submit_pridict_batches", _spy)
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    monkeypatch.setattr(adapter, "_load_designs", lambda output_dir, name: [_pegp(name, 1.0)])
    adapter.design_pair(editseq1="x", editseq2="y", name1="e1", name2="e2", walltime_hours=12)
    assert seen.get("walltime_hours") == 12


def test_real_adapter_raises_on_failed_batch(monkeypatch):
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    monkeypatch.setattr(
        ps,
        "submit_pridict_batches",
        lambda **kw: BatchResult(submitted=[], failed=[object()], output_dirs=["/out"]),
    )
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    with pytest.raises(RuntimeError, match="PRIDICT batch failed"):
        adapter.design_pair(editseq1="x", editseq2="y", name1="e1", name2="e2")


def test_real_adapter_raises_when_no_designs_found(monkeypatch):
    from pegrna_pipeline.hpc import pridict_submitter as ps
    from pegrna_pipeline.hpc.pridict_submitter import BatchResult
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    monkeypatch.setattr(
        ps,
        "submit_pridict_batches",
        lambda **kw: BatchResult(submitted=[object()], failed=[], output_dirs=["/out"]),
    )
    adapter = TwoPegPridictReal(
        ssh=_FakeSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    monkeypatch.setattr(adapter, "_load_designs", lambda output_dir, name: [])
    with pytest.raises(RuntimeError, match="no PRIDICT pegRNA designs"):
        adapter.design_pair(editseq1="x", editseq2="y", name1="e1", name2="e2")


def test_real_adapter_load_designs_returns_empty_when_output_csv_missing():
    from pegrna_pipeline.hpc.recoded_pridict import TwoPegPridictReal

    class MissingCSVSSH(_FakeSSH):
        def download(self, remote, local):
            raise FileNotFoundError(remote)

    adapter = TwoPegPridictReal(
        ssh=MissingCSVSSH(), remote_run_dir="/run", pridict_install_path="/PRIDICT2"
    )
    assert adapter._load_designs("/out", "e1") == []
