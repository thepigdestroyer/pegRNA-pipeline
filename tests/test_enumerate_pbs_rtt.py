"""Task 11: enumerate PBS1-length x synonymous-RTT variants (spacer fixed)."""

from Bio.Seq import Seq

from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    PbsRttVariant,
    enumerate_pbs_rtt_variants,
)

_SPACER = "ACGTACGTACGTACGTACGT"


def test_enumeration_fixes_spacer_varies_pbs_and_rtt():
    out = enumerate_pbs_rtt_variants(
        spacer=_SPACER, rtt_aa="MKP",
        pbs_lengths=(11, 13, 15, 17), max_rtt_variants=8,
    )
    assert all(isinstance(v, PbsRttVariant) for v in out)
    assert all(v.spacer == _SPACER for v in out)
    assert {v.pbs_length for v in out} <= {11, 13, 15, 17}
    assert all(str(Seq(v.rtt_dna).translate()) == "MKP" for v in out)
    assert len(out) <= 8 * 4


def test_enumeration_is_deterministic():
    a = enumerate_pbs_rtt_variants(spacer=_SPACER, rtt_aa="MK", pbs_lengths=(13,))
    b = enumerate_pbs_rtt_variants(spacer=_SPACER, rtt_aa="MK", pbs_lengths=(13,))
    assert [(v.rtt_dna, v.pbs_length) for v in a] == [(v.rtt_dna, v.pbs_length) for v in b]


def test_enumeration_cartesian_product_sizes():
    out = enumerate_pbs_rtt_variants(
        spacer=_SPACER, rtt_aa="K", pbs_lengths=(11, 13), max_rtt_variants=2,
    )
    # K has 2 codons (AAA/AAG); 2 rtt x 2 lengths = up to 4
    assert len(out) <= 4
    assert {v.pbs_length for v in out} == {11, 13}
