"""ADR 0003: cds_mask-aware arm partition (byte-exact non-coding) for the
recoded TJ-PE pool-walk. RED-first.

The re-templated arm genome[S:nick) must be partitioned by cds_mask:
  * CDS runs   -> synonymous recode (microhomology <= max_run)
  * non-coding -> BYTE-EXACT (Decision 15; zero edits), microhomology advisory
"""

from Bio.Seq import Seq

from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.recoded_tjpe_workflow import _partition_arm, run_recoded_tjpe
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics


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


def test_all_cds_arm_recodes_and_is_protein_silent():
    cdsarm = "GCTGCTGCTGCTGCTGCTGCTGCT"  # 24 nt, repetitive -> recode will edit
    genome = "AAA" + cdsarm + "TTT"
    cds = [False] * 3 + [True] * 24 + [False] * 3
    b = _partition_arm(genome, lo=3, hi=27, cds_mask=cds, max_run=8)
    assert b is not None
    assert b.recode_span == (3, 27)
    assert b.locked_span is None
    assert b.protein_ok is True
    assert len(b.sense_arm) == 24
    assert str(Seq(b.sense_arm).translate()) == str(Seq(cdsarm).translate())


def test_all_noncoding_arm_is_byte_exact():
    arm = "CGGAGCCTTGTGCGTCGCC" * 2  # 38 nt, repetitive (recode WOULD edit if coding)
    genome = "AAA" + arm + "TTT"
    n = len(arm)
    cds = [False] * (3 + n + 3)
    b = _partition_arm(genome, lo=3, hi=3 + n, cds_mask=cds, max_run=8)
    assert b is not None
    assert b.sense_arm == arm  # BYTE-EXACT -- the fix
    assert b.n_edits == 0
    assert b.recode_span is None
    assert b.locked_span == (3, 3 + n)
    assert b.locked_expected == arm
    assert b.microhomology_flag is True  # advisory, not a gate


def test_mixed_cds_then_noncoding_partitions():
    cdsarm = "GCTGCTGCTGCTGCTGCTGCTGCT"  # 24 nt CDS
    utr = "CGGAGCCTTGTGCGTCGCC"  # 19 nt non-coding
    genome = "AAA" + cdsarm + utr + "TTT"
    lo, hi = 3, 3 + 24 + len(utr)
    cds = [False] * 3 + [True] * 24 + [False] * len(utr) + [False] * 3
    b = _partition_arm(genome, lo=lo, hi=hi, cds_mask=cds, max_run=8)
    assert b is not None
    assert b.recode_span == (3, 27)
    assert b.locked_span == (27, hi)
    assert b.locked_expected == utr
    assert b.sense_arm[24:] == utr  # non-coding tail untouched
    assert str(Seq(b.sense_arm[:24]).translate()) == str(Seq(cdsarm).translate())
    assert b.microhomology_flag is True


def test_non_monotone_arm_is_rejected():
    # CDS / non-coding / CDS -> >2 runs -> unsupported -> None (fail-loud upstream)
    genome = "AAA" + "GCT" * 4 + "TTTTTT" + "GCT" * 4 + "AAA"
    lo, hi = 3, 3 + 12 + 6 + 12
    cds = [False] * 3 + [True] * 12 + [False] * 6 + [True] * 12 + [False] * 3
    assert _partition_arm(genome, lo=lo, hi=hi, cds_mask=cds, max_run=8) is None


def test_pool_walk_sense_layout_keeps_noncoding_arm_byte_exact():
    """End-to-end: a sense-layout window (CDS 5' of S, 3'UTR 3' of S) drives a
    valid pair whose non-coding arm is re-templated byte-exact (zero recode
    edits) with a microhomology advisory -- the sense-window het-edit capability."""
    genome = "ACGATGCATGGCATGACTGA" * 16  # 320 nt; same fixture as the passing Rv test
    site_s = 150
    cds_mask = [i < site_s for i in range(320)]  # CDS upstream of S; 3'UTR downstream
    res = run_recoded_tjpe(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome=genome,
        site_s=site_s,
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
    assert res.pegRNA_sequence
    assert res.recode_edits == 0  # non-coding arm is NOT recoded
    assert any("byte-exact" in w for w in res.warnings)  # Decision 15 advisory
    assert res.nick_position > site_s  # arm sits 3' of S (the non-coding side)
