"""T7-G / rcT7-G frame-safe insertable promoter elements + cargo-insert assembly.

T7 RNA polymerase initiates at a guanine (consensus TAATACGACTCACTATAGGGAGA,
transcription begins at the first G of GGG; Rong et al. 1998 PNAS 95:515). The
element is the coding strand, so the initiator is G on the sequence as written.
The 20-nt core is padded with the minimum G's to a multiple of 3 so it never
frameshifts a downstream CDS while keeping the start G-rich (C would be wrong).
"""

from Bio.Seq import Seq

import pytest

from pegrna_pipeline.hpc.insert_elements import (
    RC_T7G_PROMOTER,
    T7G_PROMOTER,
    T7_CORE,
    InsertElement,
    assemble_insert,
    get_promoter,
    pad_t7_to_frame,
    transcribed_segments,
)


def test_t7_core_canonical_sequence():
    assert T7_CORE == "TAATACGACTCACTATAGGG"  # canonical phi10 +1 GGG guanine start


def test_t7g_is_core_plus_frame_guanines():
    # 20 nt core -> +1 G -> 21 nt, divisible by 3 (no frameshift), start stays G-rich.
    assert T7G_PROMOTER == "TAATACGACTCACTATAGGGG"
    assert T7G_PROMOTER.startswith(T7_CORE)
    added = T7G_PROMOTER[len(T7_CORE) :]
    assert set(added) == {"G"} and 1 <= len(added) <= 3  # guanines, not cytosines
    assert len(T7G_PROMOTER) % 3 == 0  # frame-safe


def test_pad_t7_is_guanine_and_in_frame():
    # Padding is the minimum G's to reach a multiple of 3 (0 if already in-frame).
    for core in ["A", "AC", "ACG", "ACGT", "ACGTA", T7_CORE]:
        padded = pad_t7_to_frame(core)
        added = padded[len(core) :]
        assert padded.startswith(core)
        assert set(added) <= {"G"}  # only ever G's (never C)
        assert len(padded) % 3 == 0
    assert pad_t7_to_frame("ACG") == "ACG"  # already in-frame -> unchanged


def test_rct7g_is_reverse_complement_of_t7g_not_reverse():
    # The reverse element MUST be the reverse complement of the padded promoter,
    # so the functional promoter reads 5'->3' on the non-targeted strand.
    assert RC_T7G_PROMOTER == str(Seq(T7G_PROMOTER).reverse_complement())
    assert RC_T7G_PROMOTER == "CCCCTATAGTGAGTCGTATTA"
    # NOT merely the reversed string (that would break directionality):
    assert RC_T7G_PROMOTER != T7G_PROMOTER[::-1]
    # Round-trips, and starts with the C's that pair the GGGG start on the far strand:
    assert str(Seq(RC_T7G_PROMOTER).reverse_complement()) == T7G_PROMOTER
    assert RC_T7G_PROMOTER.startswith("CCCC")
    assert len(RC_T7G_PROMOTER) % 3 == 0  # also frame-safe


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
