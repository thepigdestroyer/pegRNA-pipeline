"""Task 13: ranked-pool-walk orchestrator (run_recoded_tjpe)."""

import pytest

from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    RecodedTJPEResult,
    _insertion_site_footprint,
    _resolve_protein_length,
    run_recoded_tjpe,
)
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics


def _all_pass_importance():
    class _S:
        def metrics(self, i):  # never blocks
            from pegrna_pipeline.hpc.structure_importance import ResidueMetrics

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


def _clean_offtarget(protospacers):
    return {p: (0.8, 0, False) for p in protospacers}


def _inputs(**over):
    base = dict(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome="ACGATGCATGGCATGACTGA" * 16,  # 320 nt, NGG-rich on Fw, in-band PBS2
        site_s=150,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        pegrna_strand="Rv",
        max_run=8,
    )
    base.update(over)
    return base


def test_pool_walk_returns_first_survivor():
    res = run_recoded_tjpe(
        **_inputs(),
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_all_pass_importance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=_clean_offtarget,
    )
    assert isinstance(res, RecodedTJPEResult)
    assert res.pegRNA_sequence
    assert res.deletion_size == 0
    assert res.candidates_tried >= 1
    assert res.site_s == 150  # S unchanged


def test_pool_walk_advances_past_splice_failure():
    # splice always fails => no survivor => fail-loud (proves the gate is wired)
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(
            **_inputs(),
            splice_scorer=SpliceFake(delta=0.9),
            full_structure_scorer=FullStructureFake(pairs={}),
            importance_scorer=_all_pass_importance(),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
            offtarget_scorer=_clean_offtarget,
        )


def test_pool_walk_fail_loud_on_docking_collapse():
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(
            **_inputs(),
            splice_scorer=SpliceFake(delta=0.0),
            full_structure_scorer=FullStructureFake(pairs={}),
            importance_scorer=_all_pass_importance(),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.20),
            offtarget_scorer=_clean_offtarget,
        )


def test_fail_loud_has_diagnostic():
    try:
        run_recoded_tjpe(
            **_inputs(),
            splice_scorer=SpliceFake(delta=0.9),
            full_structure_scorer=FullStructureFake(pairs={}),
            importance_scorer=_all_pass_importance(),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
            offtarget_scorer=_clean_offtarget,
        )
        raise AssertionError("expected NoValidPegRNA")
    except NoValidPegRNA as exc:
        assert "splice" in str(exc).lower()


class _RaisingSplice:
    def delta(self, ref, alt):
        raise RuntimeError("scorer exploded")


def test_raising_scorer_does_not_crash_the_walk():
    # A gate scorer that throws must be caught -> diagnostic -> NoValidPegRNA,
    # NOT propagate a raw exception that aborts the whole run.
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(
            **_inputs(),
            splice_scorer=_RaisingSplice(),
            full_structure_scorer=FullStructureFake(pairs={}),
            importance_scorer=_all_pass_importance(),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
            offtarget_scorer=_clean_offtarget,
        )


class _ShortContextSplice:
    """SpliceAI adapter behaviour on short context: delta() raises ValueError."""

    def delta(self, ref, alt):
        raise ValueError("SpliceAIScorer.delta requires >= 10000 nt of context")


def test_spliceai_short_context_is_advisory_not_fatal():
    # The CTCF-wall fix: SpliceAI can't run on the short design window, but the
    # MaxEnt cryptic gate is the primary check, so a short-context SpliceAI must be
    # skipped (advisory), NOT abort the run. Previously ValueError failed every cand.
    res = run_recoded_tjpe(
        **_inputs(),
        splice_scorer=_ShortContextSplice(),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_all_pass_importance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=_clean_offtarget,
        maxent_scorer=MaxEntScanFake(),  # primary gate: all motifs below threshold
    )
    assert res.pegRNA_sequence
    assert res.candidates_tried >= 1


def test_maxent_gate_wired_and_happy_path_passes():
    res = run_recoded_tjpe(
        **_inputs(),
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_all_pass_importance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=_clean_offtarget,
        maxent_scorer=MaxEntScanFake(),
    )
    assert isinstance(res, RecodedTJPEResult)
    assert res.pegRNA_sequence


def test_allele_hunt_surfaces_options_without_moving_s():
    het = [{"element": "pegrna_pam", "offset": 1, "haplotype": 1, "phased": True}]
    res = run_recoded_tjpe(
        **_inputs(),
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_all_pass_importance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=_clean_offtarget,
        allele_hunt=True,
        het_snps=het,
    )
    assert res.site_s == 150  # S NEVER moves, even in allele_hunt
    assert res.allele_specific_options  # at least one surfaced


def test_insertion_site_footprint_flanks_s():
    # S=150, all-CDS: residues of 149 (49) and 150 (50), +/-1 context -> 48..51.
    assert _insertion_site_footprint(150, [True] * 320, context=1) == [48, 49, 50, 51]


def test_insertion_site_footprint_context_zero():
    assert _insertion_site_footprint(150, [True] * 320, context=0) == [49, 50]


def test_insertion_site_footprint_noncoding_returns_empty():
    assert _insertion_site_footprint(5, [False] * 10, context=1) == []


def test_insertion_site_footprint_maps_local_to_absolute_with_protein_length():
    # 300 bp all-CDS window; last coding residue local index = 299//3 = 99 anchored to
    # the C-terminus of a 727-aa protein (residue 726). S=150 -> local residues 49/50 ->
    # absolute 49+627 / 50+627 = 676/677 (offset 726-99=627), +/-1 context -> 675..678.
    fp = _insertion_site_footprint(150, [True] * 300, context=1, protein_length=727)
    assert fp == [675, 676, 677, 678]
    # Without protein_length the same window stays window-LOCAL (the pre-fix behaviour
    # that silently indexed local residues into the full-protein structure).
    assert _insertion_site_footprint(150, [True] * 300, context=1) == [48, 49, 50, 51]


class _LenScorer:
    def __init__(self, n):
        self._n = n

    def protein_length(self):
        return self._n

    def metrics(self, i):  # pragma: no cover - not exercised here
        raise AssertionError("metrics not needed")


def test_resolve_protein_length_prefers_explicit():
    assert _resolve_protein_length(727, _LenScorer(999)) == 727


def test_resolve_protein_length_uses_scorer_when_none():
    assert _resolve_protein_length(None, _LenScorer(727)) == 727


def test_resolve_protein_length_none_when_scorer_has_no_length():
    class _NoLen:
        def metrics(self, i):  # pragma: no cover
            raise AssertionError

    assert _resolve_protein_length(None, _NoLen()) is None


def test_resolve_protein_length_none_on_scorer_error():
    class _Boom:
        def protein_length(self):
            raise RuntimeError("no AFDB entry")

    assert _resolve_protein_length(None, _Boom()) is None


class _BlockResidues:
    """Importance scorer that marks the given residue indices as buried + packed."""

    def __init__(self, blocked):
        self._b = set(blocked)

    def metrics(self, i):
        buried = i in self._b
        return ResidueMetrics(
            plddt=90.0 if buried else 50.0,
            rsa=0.05 if buried else 0.5,
            contact_number=30 if buried else 4,
            is_catalytic=False,
            is_disulfide=False,
            is_positive_phi_gly=False,
            is_structural_pro=False,
        )


def test_importance_gate_blocks_structured_insertion_site():
    # S=150 -> insertion footprint [48..51]; a buried residue there blocks all cands.
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(
            **_inputs(),
            splice_scorer=SpliceFake(delta=0.0),
            full_structure_scorer=FullStructureFake(pairs={}),
            importance_scorer=_BlockResidues([50]),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
            offtarget_scorer=_clean_offtarget,
            maxent_scorer=MaxEntScanFake(),
        )


def test_importance_gate_ignores_recode_arm_residues():
    # Residue 55 is in the recoded arm (the old gate would have blocked it) but NOT
    # in the insertion footprint [48..51]. The gate now judges the insertion site
    # only, so the synonymous recode arm no longer blocks the design.
    res = run_recoded_tjpe(
        **_inputs(),
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_BlockResidues([55]),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=_clean_offtarget,
        maxent_scorer=MaxEntScanFake(),
    )
    assert res.pegRNA_sequence
