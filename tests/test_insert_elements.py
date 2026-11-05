"""T7-G / rcT7-G insertable promoter elements (two G's after TATA) + cargo-insert assembly.

T7 RNA polymerase initiates at a guanine (consensus TAATACGACTCACTATAGGGAGA, transcription begins at the first G
after TATA; Rong et al. 1998 PNAS 95:515). The element is the coding strand, so the initiator is G on the sequence
as written. The element is the promoter plus the +1 G and one more: ``TAATACGACTCACTATAGG``, 19 nt. An earlier
version carried ``GGGG`` (21 nt, a multiple of 3); at the 3' end of the CTCF cassette that GGGG is the reverse
complement of CCCC and paired with the rcPBS2's GGGG, which failed the folding gate on every candidate (ADR 0028
addendum). The element is therefore NOT a multiple of 3 and no longer keeps a reading frame on its own.
"""

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.insert_elements import (
    RC_T7G_PROMOTER,
    T7_CORE,
    T7_MIN_PROMOTER,
    T7G_PROMOTER,
    InsertElement,
    assemble_insert,
    frame_shift_after,
    get_promoter,
    pad_t7_to_frame,
    transcribed_segments,
)


def test_t7_element_is_the_promoter_with_exactly_two_guanines_after_tata():
    assert T7_CORE == "TAATACGACTCACTATAGG"
    assert T7_CORE.startswith(T7_MIN_PROMOTER) and T7_MIN_PROMOTER == "TAATACGACTCACTATAG"
    assert T7_CORE.split("TATA")[-1] == "GG"  # two G's, the +1 and one more, and no more


def test_t7g_is_the_two_g_element_and_ends_in_gg_not_gggg():
    assert T7G_PROMOTER == "TAATACGACTCACTATAGG" and len(T7G_PROMOTER) == 19
    assert T7G_PROMOTER.endswith("TAGG") and not T7G_PROMOTER.endswith("GGG")


def test_the_element_is_not_frame_neutral_and_says_so():
    assert len(T7G_PROMOTER) % 3 == 1
    assert frame_shift_after(T7G_PROMOTER) == 1 and frame_shift_after(RC_T7G_PROMOTER) == 1
    assert frame_shift_after("ACG") == 0 and frame_shift_after("ACGT") == 1


def test_pad_t7_is_still_available_and_is_guanine_and_in_frame():
    # Kept for callers that want a frame-neutral promoter: the minimum G's to reach a multiple of 3.
    for core in ["A", "AC", "ACG", "ACGT", "ACGTA", T7_CORE]:
        padded = pad_t7_to_frame(core)
        added = padded[len(core) :]
        assert padded.startswith(core)
        assert set(added) <= {"G"}  # only ever G's (never C)
        assert len(padded) % 3 == 0
    assert pad_t7_to_frame("ACG") == "ACG"  # already in-frame -> unchanged
    assert pad_t7_to_frame(T7_CORE) == T7_CORE + "GG"  # the old, longer element is one call away


def test_rct7g_is_reverse_complement_of_t7g_not_reverse():
    # The reverse element MUST be the reverse complement, so the functional promoter reads 5'->3' on the
    # non-targeted strand.
    assert RC_T7G_PROMOTER == str(Seq(T7G_PROMOTER).reverse_complement())
    assert RC_T7G_PROMOTER == "CCTATAGTGAGTCGTATTA"
    # NOT merely the reversed string (that would break directionality):
    assert RC_T7G_PROMOTER != T7G_PROMOTER[::-1]
    # Round-trips, and starts with the two C's that pair the GG start on the far strand:
    assert str(Seq(RC_T7G_PROMOTER).reverse_complement()) == T7G_PROMOTER
    assert RC_T7G_PROMOTER.startswith("CC") and not RC_T7G_PROMOTER.startswith("CCC")
    assert len(RC_T7G_PROMOTER) == 19


def test_old_and_cytosine_names_are_gone():
    import pegrna_pipeline.hpc.insert_elements as ie

    for gone in ("T7_PROMOTER", "T7_PROMOTER_REV", "T7C_PROMOTER", "RC_T7C_PROMOTER"):
        assert not hasattr(ie, gone), gone
    assert set(ie.PROMOTERS) == {"T7-G", "rcT7-G"}  # no T7 / T7_rev / *-C keys


def test_get_promoter_directions():
    assert get_promoter("T7-G").direction == "forward"
    assert get_promoter("rcT7-G").direction == "reverse"
    assert get_promoter("T7-G").sequence == T7G_PROMOTER
    assert get_promoter("rcT7-G").sequence == RC_T7G_PROMOTER


def test_get_promoter_unknown_raises():
    with pytest.raises(KeyError):
        get_promoter("SP6")


def test_assemble_insert_concatenates_named_segments():
    segs = [
        ("FLAG", "GACTACAAAGACGATGACGACAAG"),
        ("GGGGS", "GGTGGAGGCGGTTCA"),
        ("T7-G", T7G_PROMOTER),
        ("P2A", "GGAAGCGGAGCT"),
        ("ZeoR", "ATGGCCAAGTTGACC"),
    ]
    out = assemble_insert(segs)
    assert out == "".join(s for _, s in segs)


def test_assemble_insert_rejects_non_acgt():
    with pytest.raises(ValueError, match="ACGT"):
        assemble_insert([("bad", "ACGTN")])


def test_forward_t7_transcribes_downstream():
    # 5'-FLAG-GGGGS-T7-G-P2A-ZeoR-3' : forward T7-G transcribes P2A-ZeoR (downstream)
    segs = [
        ("FLAG", "GACTACAAAGACGATGACGACAAG"),
        ("GGGGS", "GGTGGAGGCGGTTCA"),
        ("T7-G", T7G_PROMOTER),
        ("P2A", "GGAAGCGGAGCT"),
        ("ZeoR", "ATGGCCAAGTTGACC"),
    ]
    assert transcribed_segments(segs, promoter_index=2) == ["P2A", "ZeoR"]


def test_reverse_t7_transcribes_upstream():
    # 5'-FLAG-GGGGS-rcT7-G-P2A-ZeoR-3' : rcT7-G transcribes GGGGS-FLAG (upstream)
    segs = [
        ("FLAG", "GACTACAAAGACGATGACGACAAG"),
        ("GGGGS", "GGTGGAGGCGGTTCA"),
        ("rcT7-G", RC_T7G_PROMOTER),
        ("P2A", "GGAAGCGGAGCT"),
        ("ZeoR", "ATGGCCAAGTTGACC"),
    ]
    # upstream segments, reported 5'->3' from the promoter outward
    assert transcribed_segments(segs, promoter_index=2) == ["GGGGS", "FLAG"]


def test_insert_element_is_frozen():
    el = InsertElement(name="T7-G", sequence=T7G_PROMOTER, direction="forward")
    with pytest.raises(Exception):
        el.sequence = "X"  # frozen dataclass
