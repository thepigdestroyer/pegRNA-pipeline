"""PRIDICT-driven ngRNA wiring (Stage 2b of the 2026-07-10 PRIDICT2 wiring spec).

The ngRNA (pegRNA#2 spacer) comes from PRIDICT2's ``*_nicking_guides.csv`` rows,
while the PIPELINE owns *where* it nicks: it selects among the PRIDICT nicking
candidates by DISTANCE WINDOW + opposite strand + DeepCas9 score
(``tjpe_workflow.select_tjpe_ngrna``), gates the choice on off-target
(the existing off-target scorer), and on a gate failure advances to the NEXT
nicking candidate (the NICK lever). rcPBS2 stays genome-derived via
``tjpe_workflow.derive_rc_pbs2``.

These tests pin the additive path:

* the ngRNA is selected from PRIDICT nicking guides (window + opposite strand + score);
* an off-target failure advances to the next candidate (NICK lever);
* rcPBS2 is derived from the genome 5' of the nick;
* the window helper uses the provided window when it has >= min_flank each side,
  else auto-fetches a wider GRCh38 window and localises the provided window by
  exact substring match (fail-loud when neither yields adequate flanks);
* ``run_recoded_tjpe`` consumes an injected PRIDICT ngRNA (+ rcPBS2) and the
  existing genome-enumerated path is untouched when none is injected;
* ``design_recoded`` uses the PRIDICT ngRNA in the auto path and falls back to
  the genome-derived path when PRIDICT yields no admissible ngRNA.

No real PRIDICT2, no network: bespoke deterministic fake runners + a fake
genome fetcher throughout.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc import recoded_driver
from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanFake
from pegrna_pipeline.hpc.data_types import NgRNADesign, PegRNADesign
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.pridict_ngrna import (
    InjectedNgRNA,
    build_pridict_ngrnas,
    nicking_rows_to_ngrnas,
    resolve_pridict_context,
)
from pegrna_pipeline.hpc.recoded_tjpe_workflow import RecodedTJPEResult, run_recoded_tjpe
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics
from pegrna_pipeline.hpc.tjpe_workflow import derive_rc_pbs2

GENOME = "ACGATGCATGGCATGACTGA" * 16  # 320 bp; same fixture the workflow tests use
SITE_S = 150

# Two 20-mers that pass the ngRNA spacer sanity filter (no TTTT / homopolymer>=5,
# GC in [0.30, 0.75]); distinctive markers so a design's ngRNA can be traced to PRIDICT.
_MARKER_A = "ACGTACGTACGTACGTACGT"  # GC 0.50
_MARKER_B = "GATCGATCGATCGATCGATC"  # GC 0.50


def _nrow(proto: str, dist: str, strand: str, score: str) -> dict:
    """A PRIDICT ``*_nicking_guides.csv``-shaped row (csv.DictReader style)."""
    return {
        "Nicking-Protospacer": proto,
        "Nicking-Position-to-edit": dist,
        "Nicking-PAMdisrupt": "Nicking_PAM_disrupt",
        "Target_Strand": strand,
        "DeepCas9score": score,
    }


def _clean(protospacers):
    return {p: (0.8, 0, False) for p in protospacers}


# ============================================================ mapper


def test_nicking_rows_to_ngrnas_maps_required_fields():
    rows = [_nrow(_MARKER_A, "-30", "Fw", "0.55")]
    (ng,) = nicking_rows_to_ngrnas(rows)
    assert isinstance(ng, NgRNADesign)
    assert ng.spacer == _MARKER_A
    assert ng.distance_to_edit == -30
    assert ng.target_strand == "Fw"
    assert ng.deepcas9_score == 0.55


def test_nicking_rows_to_ngrnas_missing_key_fails_loud():
    bad = {"Nicking-Protospacer": _MARKER_A, "Target_Strand": "Fw"}  # no distance / score
    with pytest.raises(ValueError, match="DeepCas9score"):
        nicking_rows_to_ngrnas([bad])


# ============================================================ selection


def test_build_selects_admissible_and_orders_by_deepcas9():
    rows = [
        _nrow(_MARKER_B, "60", "Fw", "0.5"),  # admissible, lower score
        _nrow(_MARKER_A, "50", "Fw", "0.9"),  # admissible, higher score
    ]
    out = build_pridict_ngrnas(
        rows,
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",  # opposite strand is Fw
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        offtarget_scorer=_clean,
    )
    assert [d.deepcas9_score for d in out] == [0.9, 0.5]  # highest DeepCas9 first
    assert [d.nick_position for d in out] == [SITE_S + 50, SITE_S + 60]
    assert [d.protospacer for d in out] == [_MARKER_A, _MARKER_B]


def test_build_excludes_out_of_window_and_same_strand():
    rows = [
        _nrow(_MARKER_A, "100", "Fw", "0.9"),  # |100| outside [40,70]
        _nrow(_MARKER_B, "55", "Rv", "0.9"),  # same strand as pegRNA (Rv) -> excluded
    ]
    out = build_pridict_ngrnas(
        rows,
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        offtarget_scorer=_clean,
    )
    assert out == []


def test_offtarget_failure_advances_to_next_candidate_nick_lever():
    # The higher-DeepCas9 candidate fails off-target -> the pipeline advances to the
    # next (lower-score, off-target-clean) nicking candidate: the NICK lever.
    rows = [
        _nrow(_MARKER_A, "50", "Fw", "0.9"),  # higher score but DIRTY
        _nrow(_MARKER_B, "60", "Fw", "0.5"),  # lower score but CLEAN
    ]

    def scorer(protospacers):
        return {p: ((0.8, 5, True) if p == _MARKER_A else (0.8, 0, False)) for p in protospacers}

    out = build_pridict_ngrnas(
        rows,
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        offtarget_scorer=scorer,
    )
    assert len(out) == 1
    assert out[0].protospacer == _MARKER_B
    assert out[0].distance_to_edit == 60
    assert out[0].nick_position == SITE_S + 60


def test_build_derives_rcpbs2_from_genome_five_prime_of_nick():
    rows = [_nrow(_MARKER_A, "50", "Fw", "0.7")]
    (inj,) = build_pridict_ngrnas(
        rows,
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        offtarget_scorer=_clean,
        pbs2_length=17,
    )
    assert inj.rcpbs2 == derive_rc_pbs2(GENOME, ngrna_nick_position=SITE_S + 50, pbs2_length=17)


def test_build_empty_when_no_admissible_candidate():
    rows = [_nrow(_MARKER_A, "100", "Fw", "0.9")]  # out of window
    out = build_pridict_ngrnas(
        rows,
        genome=GENOME,
        site_s=SITE_S,
        pegrna_strand="Rv",
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        offtarget_scorer=_clean,
    )
    assert out == []


# ============================================================ window strategy


def test_context_passthrough_when_window_has_enough_flank():
    seq, ss = resolve_pridict_context(genome=GENOME, site_s=SITE_S, min_flank=100)
    assert seq == GENOME and ss == SITE_S  # 150 / 170 bp flanks -> used as-is


def test_context_autofetch_localises_and_trims():
    window = "ACGT" * 10  # 40 bp, site_s within window = 20 (too short: left 20 < 100)
    site_in_window = 20
    wide = "G" * 500 + window + "C" * 500  # window absent from the flanks

    seq, ss = resolve_pridict_context(
        genome=window,
        site_s=site_in_window,
        min_flank=100,
        transcript="ENST_FAKE",
        fetcher=lambda tid: wide,
        target_flank=175,
    )
    assert ss >= 100 and (len(seq) - ss) >= 100  # >= min_flank each side after trim
    # the provided window is recovered at the mapped offset
    assert seq[ss - site_in_window : ss - site_in_window + len(window)] == window
    assert len(seq) < len(wide)  # trimmed to ~target_flank each side


def test_context_fails_loud_without_transcript():
    with pytest.raises(ValueError, match="no transcript"):
        resolve_pridict_context(genome="ACGT" * 10, site_s=20, min_flank=100)


def test_context_fails_loud_when_window_not_localisable():
    with pytest.raises(ValueError, match="localise"):
        resolve_pridict_context(
            genome="ACGT" * 10,
            site_s=20,
            min_flank=100,
            transcript="ENST_FAKE",
            fetcher=lambda tid: "G" * 400,  # does not contain the window
        )


def test_context_fails_loud_when_fetched_still_too_short():
    window = "ACGT" * 10
    wide = "T" * 50 + window + "A" * 400  # only 70 bp 5' of the insert -> < min_flank
    with pytest.raises(ValueError, match="still too short"):
        resolve_pridict_context(
            genome=window,
            site_s=20,
            min_flank=100,
            transcript="ENST_FAKE",
            fetcher=lambda tid: wide,
        )


# ============================================================ run_recoded_tjpe injection


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


def _wf_inputs(**over):
    base = dict(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome=GENOME,
        site_s=SITE_S,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        pegrna_strand="Rv",
        max_run=8,
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_all_pass_importance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=_clean,
        maxent_scorer=MaxEntScanFake(),
    )
    base.update(over)
    return base


def test_run_recoded_tjpe_uses_injected_ngrna_and_its_rcpbs2():
    nick = SITE_S + 50  # 200, downstream (arm_side=3prime default)
    rcpbs2 = derive_rc_pbs2(GENOME, ngrna_nick_position=nick, pbs2_length=17)
    injected = InjectedNgRNA(
        nick_position=nick,
        protospacer=_MARKER_B,
        rcpbs2=rcpbs2,
        distance_to_edit=50,
        target_strand="Fw",
        deepcas9_score=0.7,
    )
    res = run_recoded_tjpe(**_wf_inputs(), injected_ngrnas=[injected])
    assert isinstance(res, RecodedTJPEResult)
    assert res.nick_position == nick  # pipeline mapped the PRIDICT distance
    assert res.ngrna_protospacer == _MARKER_B  # ngRNA came from PRIDICT
    assert res.rcpbs2 == rcpbs2  # rcPBS2 = derive_rc_pbs2, injected verbatim
    assert res.site_s == SITE_S  # S never moves


def test_run_recoded_tjpe_genome_path_unchanged_when_no_injection():
    res = run_recoded_tjpe(**_wf_inputs())
    assert res.pegRNA_sequence
    # genome-enumerated protospacer is a real substring of the window, not a marker
    assert res.ngrna_protospacer not in {_MARKER_A, _MARKER_B}


# ============================================================ design_recoded end-to-end


def _peg(spacer: str) -> PegRNADesign:
    return PegRNADesign(
        sequence_name="peg1_p1",
        spacer=spacer,
        pbs="TCGACGCTTAGGAACAC",
        rt="ACGTACGTACGTACGTACGTACGTACGTAC",
        pbs_length=13,
        rt_length=30,
        target_strand="Fw",
        editing_position=50,
        pe2_hek293t=0.5,
        pe2_k562=0.5,
        pegrna_full="",
    )


class _NickRunner:
    """A PRIDICT stand-in emitting a chosen pegRNA#1 spacer + nicking rows."""

    def __init__(self, spacer: str, nicking: list) -> None:
        self._spacer = spacer
        self._nicking = nicking

    def run_batch(self, rows):
        name = rows[0]["sequence_name"]
        return {name: {"pegrnas": [_peg(self._spacer)], "nicking": self._nicking}}


def _auto_inputs() -> dict:
    return dict(
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome=GENOME,  # 320 bp -> 150 / 170 flanks at site_s=150
        site_s=SITE_S,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
    )


def test_design_recoded_auto_uses_pridict_ngrna():
    runner = _NickRunner(
        spacer="TGCATGCATGCATGCATGCA",
        nicking=[_nrow(_MARKER_B, "50", "Fw", "0.7")],  # in-window (|50| in [40,70]), Fw
    )
    out = recoded_driver.design_recoded(**_auto_inputs(), mode="fake", pridict_runner=runner)
    assert out.result.spacer == "TGCATGCATGCATGCATGCA"  # PRIDICT-designed pegRNA#1 spacer
    assert out.result.ngrna_protospacer == _MARKER_B  # ngRNA from the PRIDICT nicking guide
    assert out.result.nick_position == SITE_S + 50
    assert out.result.rcpbs2 == derive_rc_pbs2(
        GENOME, ngrna_nick_position=SITE_S + 50, pbs2_length=17
    )


def test_design_recoded_auto_falls_back_to_genome_ngrna_when_no_pridict_nick():
    # No admissible PRIDICT nicking guide -> the genome-enumerated ngRNA path is used
    # (the design still succeeds; the ngRNA is not the PRIDICT marker).
    runner = _NickRunner(spacer="TGCATGCATGCATGCATGCA", nicking=[])
    out = recoded_driver.design_recoded(**_auto_inputs(), mode="fake", pridict_runner=runner)
    assert out.result.pegRNA_sequence
    assert out.result.spacer == "TGCATGCATGCATGCATGCA"
    assert out.result.ngrna_protospacer not in {_MARKER_A, _MARKER_B}
