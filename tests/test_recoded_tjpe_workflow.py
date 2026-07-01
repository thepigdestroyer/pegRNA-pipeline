"""Task 13: ranked-pool-walk orchestrator (run_recoded_tjpe)."""

import pytest

from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    RecodedTJPEResult,
    run_recoded_tjpe,
)
from pegrna_pipeline.hpc.splice_screen import SpliceFake


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
