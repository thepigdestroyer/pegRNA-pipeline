"""Task 11b: authoritative tjPE pegRNA assembly order (decision 0 / B1).

5'->3':  Spacer . Scaffold . PBS2 . RTT(insert+deletion-replacement-arm) . PBS1 . PolyU

Implemented as a NEW function ``assemble_tj_pegrna_v2`` so the legacy
``assemble_tj_pegrna`` (consumed by the old scarless ``run_tjpe_workflow`` and
its passing tests) is left byte-for-byte unchanged. The recoded workflow uses
the v2 layout. See the deviation note in the Task 11b commit message.
"""

import pytest

from pegrna_pipeline.hpc.tjpe_workflow import (
    STANDARD_SGRNA_SCAFFOLD,
    assemble_tj_pegrna_v2,
)


def test_authoritative_order_spacer_scaffold_pbs2_rtt_pbs1_polyu():
    out = assemble_tj_pegrna_v2(
        spacer="ACGTACGTACGTACGTACGT",
        pbs2="GCAGTCGATCG",
        rtt="ATGAAACCCGGGTTTAAATTT",
        pbs1="CCCAGTAGCTAGC",
        polyu="UUU",
        scaffold=STANDARD_SGRNA_SCAFFOLD,
    )
    expected = (
        "ACGTACGTACGTACGTACGT" + STANDARD_SGRNA_SCAFFOLD
        + "GCAGTCGATCG" + "ATGAAACCCGGGTTTAAATTT" + "CCCAGTAGCTAGC" + "UUU"
    )
    assert out == expected
    assert out.index("GCAGTCGATCG") < out.index("ATGAAACCCGGGTTTAAATTT") < out.index("CCCAGTAGCTAGC")
    assert out.endswith("UUU")


def test_polyu_must_be_only_uridines():
    with pytest.raises(ValueError, match="uridine"):
        assemble_tj_pegrna_v2(spacer="ACGT" * 5, pbs2="GCAGTCGATCG",
                              rtt="ATGAAA", pbs1="CCCAGTAGCTAGC", polyu="UUA")


def test_polyu_length_2_to_3():
    # 2 U ok
    out = assemble_tj_pegrna_v2(spacer="ACGT" * 5, pbs2="GCAGTCGATCG",
                                rtt="ATGAAA", pbs1="CCCAGTAGCTAGC", polyu="UU")
    assert out.endswith("UU")
    with pytest.raises(ValueError, match="2-3"):
        assemble_tj_pegrna_v2(spacer="ACGT" * 5, pbs2="GCAGTCGATCG",
                              rtt="ATGAAA", pbs1="CCCAGTAGCTAGC", polyu="UUUU")


def test_dna_components_reject_non_acgt():
    with pytest.raises(ValueError):
        assemble_tj_pegrna_v2(spacer="ACGT" * 5, pbs2="GCAGTCGATCN",  # N
                              rtt="ATGAAA", pbs1="CCCAGTAGCTAGC", polyu="UUU")
