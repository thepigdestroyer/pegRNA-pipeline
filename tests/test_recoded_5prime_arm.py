"""ADR 0003 / Option A: 5'-of-S recode arm for C-terminal recode-into-CDS.

The cargo always inserts at S; the recode/deletion-replacement arm may sit 5' of
S (into the CDS) so the last codons are recoded synonymously in the SENSE frame.
``_assemble_edit`` is the pure geometry helper; the orchestrator gates ``arm_side``.
"""

from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    _assemble_edit,
    run_recoded_tjpe,
)
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics

GENOME = "AAAAAAAAAA" + "CCCCCCCCCC" + "GGGGGGGGGG"  # 30 nt, deterministic


class _ImpAllPass:
    def metrics(self, i):
        return ResidueMetrics(
            plddt=50.0,
            rsa=0.5,
            contact_number=4,
            is_catalytic=False,
            is_disulfide=False,
            is_positive_phi_gly=False,
            is_structural_pro=False,
        )


def test_assemble_edit_3prime_cargo_before_arm():
    sf, ed, ts, th = _assemble_edit(
        genome=GENOME, site_s=10, nick_position=20, arm_sense="TTTTTTTTTT", cargo="ACGT"
    )
    assert sf == "ACGT" + "TTTTTTTTTT"  # cargo, then 3' arm
    assert ed == GENOME[:10] + "ACGT" + "TTTTTTTTTT" + GENOME[20:]
    assert (ts, th) == (10, 20)  # tprt boundaries == (S, nick)


def test_assemble_edit_5prime_arm_before_cargo():
    # nick 5' of S: arm = genome[10:20] replaced by recoded arm, cargo still at S=20
    sf, ed, ts, th = _assemble_edit(
        genome=GENOME, site_s=20, nick_position=10, arm_sense="TTTTTTTTTT", cargo="ACGT"
    )
    assert sf == "TTTTTTTTTT" + "ACGT"  # recoded arm, then cargo
    assert ed == GENOME[:10] + "TTTTTTTTTT" + "ACGT" + GENOME[20:]
    assert (ts, th) == (10, 20)  # tprt boundaries == (nick, S)
    # cargo lands exactly at S in the edited allele
    assert ed[20:24] == "ACGT"


# Deterministic fixture (verified to yield a surviving 5'-of-S candidate): CDS
# [0,101) with a recodable Leu arm, a Fw PAM giving an ngRNA nick at 51, and a
# PBS2 window that clears the Tm/GC gate.
_G5 = (
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGGTTATCTTCGGAT"
    "CTGCTGGTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG" + "T" * 160 + "AAAAAAAAAA"
)
_S5 = 101


def test_5prime_path_recodes_cds_and_preserves_protein():
    """End-to-end: arm_side='5prime' picks a nick 5' of S, recodes the CDS arm
    synonymously (>0 edits), and survives INV-0 + TPRT (proves the geometry and
    that the recode is protein-silent in the SENSE frame)."""
    cds_mask = [i < _S5 for i in range(len(_G5))]
    res = run_recoded_tjpe(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="L",
        cargo="ATGAAACCCGGG",
        genome=_G5,
        site_s=_S5,
        cds_mask=cds_mask,
        pbs_lengths=(15,),
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        pegrna_strand="Rv",
        arm_side="5prime",
        max_run=8,
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_ImpAllPass(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=lambda ps: {p: (0.8, 0, False) for p in ps},
    )
    assert res.nick_position < res.site_s  # arm sits 5' of S (into the CDS)
    assert res.recode_edits > 0  # silent mutations actually made in the CDS
    assert res.site_s == _S5  # S never moves


def test_default_arm_side_excludes_5prime_nicks():
    # default arm_side='3prime' must reject a purely-5' fixture (no 3' nick) -> fail-loud
    arm = "CTG" * 20
    genome = "ATGGCAACCGGTACACCTAGC" + arm + "AA"  # no 3'-of-S room
    s = len(genome) - 2
    cds_mask = [i < s for i in range(len(genome))]
    import pytest

    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(
            spacer="ACGTACGTACGTACGTACGT",
            rtt_aa="L",
            cargo="ATGAAACCCGGG",
            genome=genome,
            site_s=s,
            cds_mask=cds_mask,
            pbs_lengths=(15,),
            ngrna_distance_min=40,
            ngrna_distance_max=70,
            pegrna_strand="Rv",
            max_run=8,
            splice_scorer=SpliceFake(delta=0.0),
            full_structure_scorer=FullStructureFake(pairs={}),
            importance_scorer=_ImpAllPass(),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
            offtarget_scorer=lambda ps: {p: (0.8, 0, False) for p in ps},
        )
