"""Stage F: the INSERT (synonymous-cargo) lever + the every-gate-invoked guard (F2).

``run_recoded_tjpe`` is a thin wrapper over ``_run_recoded_tjpe_single`` that retries with
the next CAI-ranked SYNONYMOUS cargo variant ONLY on an insert-side failure (splice /
cryptic-splice / RNA structure -- the gates a cargo SEQUENCE affects), never on a
cargo-independent one (importance = genomic site, off-target, INV-0, nick). The routing is
pinned by monkeypatching the core (so it does not depend on cargo-dependent fake scorers).
F2 pins that a full pass actually invokes every gate scorer.
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc import recoded_tjpe_workflow as rw
from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA, run_recoded_tjpe
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics

_CLEAN_CARGO = "ATGAAACCCGGG"  # ATG AAA CCC GGG -> MKPG, clean in-frame CDS, no stop


# --------------------------------------------------------------- INSERT-lever routing


def test_insert_lever_retries_on_insert_side_failure(monkeypatch):
    calls = []

    def _fail_then_win(*, cargo, **_kw):
        calls.append(cargo)
        if len(calls) == 1:
            raise NoValidPegRNA(
                "no candidate ... cand1: RiNALMo cross-region pairing",
                diagnostics=["cand1: RiNALMo cross-region pairing"],
            )
        return "WINNER"

    monkeypatch.setattr(rw, "_run_recoded_tjpe_single", _fail_then_win)
    out = run_recoded_tjpe(cargo=_CLEAN_CARGO, max_insert_variants=6)
    assert out == "WINNER"
    assert len(calls) == 2  # retried with a second cargo variant
    assert calls[0] != calls[1]  # a DIFFERENT (synonymous) cargo
    assert str(Seq(calls[0]).translate()) == str(Seq(calls[1]).translate())  # synonymous


def test_insert_lever_does_not_retry_on_cargo_independent_failure(monkeypatch):
    calls = []

    def _fail_importance(*, cargo, **_kw):
        calls.append(cargo)
        raise NoValidPegRNA(
            "no candidate ... importance",
            diagnostics=["cand1: insertion site residue(s) [7] structurally important (buried)"],
        )

    monkeypatch.setattr(rw, "_run_recoded_tjpe_single", _fail_importance)
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(cargo=_CLEAN_CARGO, max_insert_variants=6)
    assert len(calls) == 1  # a synonymous swap cannot rescue a buried genomic site -> no retry


def test_insert_lever_exhausts_variants_then_raises(monkeypatch):
    calls = []

    def _always_insert_fail(*, cargo, **_kw):
        calls.append(cargo)
        raise NoValidPegRNA("... cryptic", diagnostics=["cand1: cryptic splice: created"])

    monkeypatch.setattr(rw, "_run_recoded_tjpe_single", _always_insert_fail)
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(cargo=_CLEAN_CARGO, max_insert_variants=4)
    assert 1 < len(calls) <= 4  # tried several variants, bounded by the cap


def test_cargo_variants_noop_for_non_clean_cds():
    assert rw._cargo_variants("ATGAA", 6) == ["ATGAA"]  # not a multiple of 3 -> lever off
    vs = rw._cargo_variants(_CLEAN_CARGO, 5)
    prot = str(Seq(_CLEAN_CARGO).translate())
    assert vs[0] == _CLEAN_CARGO and len(vs) >= 2
    assert all(str(Seq(v).translate()) == prot for v in vs)  # every variant is synonymous
    assert len(set(vs)) == len(vs)  # distinct
    assert rw._cargo_variants(_CLEAN_CARGO, 1) == [_CLEAN_CARGO]  # cap 1 disables the lever


# ------------------------------------------------------------------ F2: every gate invoked


class _Rec:
    """Wrap a scorer, recording which gate methods get called."""

    def __init__(self, inner, methods):
        self._inner = inner
        self._methods = set(methods)
        self.hits = {m: 0 for m in methods}

    def __getattr__(self, name):
        attr = getattr(self._inner, name)
        if name in self._methods and callable(attr):

            def _wrapped(*a, **k):
                self.hits[name] += 1
                return attr(*a, **k)

            return _wrapped
        return attr


def _all_pass_importance():
    class _S:
        def metrics(self, i):  # noqa: ARG002
            return ResidueMetrics(
                plddt=50.0,
                rsa=0.5,
                contact_number=4,
                is_catalytic=False,
                is_disulfide=False,
                is_positive_phi_gly=False,
                is_structural_pro=False,
            )

    return _S()


def test_every_gate_scorer_is_invoked_on_a_full_pass():
    splice = _Rec(SpliceFake(delta=0.0), ["delta"])
    structure = _Rec(FullStructureFake(pairs={}), ["pair_probs"])
    importance = _Rec(_all_pass_importance(), ["metrics"])
    dock = _Rec(
        DockFake(apo_iptm=0.85, candidate_iptm=0.82), ["interface_iptm", "apo_interface_iptm"]
    )
    maxent = _Rec(MaxEntScanFake(), ["score5", "score3"])
    ot = {"hits": 0}

    def offtarget(protospacers):
        ot["hits"] += 1
        return {p: (0.8, 0, False) for p in protospacers}

    res = run_recoded_tjpe(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo=_CLEAN_CARGO,
        genome="ACGATGCATGGCATGACTGA" * 16,
        site_s=150,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
        splice_scorer=splice,
        full_structure_scorer=structure,
        importance_scorer=importance,
        dock_scorer=dock,
        offtarget_scorer=offtarget,
        maxent_scorer=maxent,
    )
    assert res.pegRNA_sequence
    # Authoritative "every gate ran" guard: the winning result records a verdict for every
    # gate in the cascade (robust to scan-window / caching details of individual scorers).
    recorded = {k for k, _ in res.gate_results}
    for gate in (
        "offtarget",
        "recode",
        "inv0",
        "cryptic",
        "splice",
        "importance",
        "structure",
        "docking",
    ):
        assert gate in recorded, f"{gate} gate produced no verdict (silently skipped?)"
    # Plus: the scorers that are called per candidate must actually have been invoked
    # (guards a gate wired to a no-op scorer). maxent.score5/score3 and the docking apo
    # baseline are intentionally NOT asserted here: the cryptic scan only calls the maxent
    # scorer where a splice-site window exists, and the apo baseline is cached.
    assert splice.hits["delta"] >= 1, "splice (SpliceAI) scorer not invoked"
    assert structure.hits["pair_probs"] >= 1, "RNA structure scorer not invoked"
    assert importance.hits["metrics"] >= 1, "AlphaFold importance scorer not invoked"
    assert dock.hits["interface_iptm"] >= 1, "docking scorer not invoked"
    assert ot["hits"] >= 1, "off-target scorer not invoked"
