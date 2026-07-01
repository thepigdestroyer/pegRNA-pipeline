"""T7 / T7-rev insertable promoter elements + cargo-insert assembly."""

from Bio.Seq import Seq

import pytest

from pegrna_pipeline.hpc.insert_elements import (
    T7_PROMOTER,
    T7_PROMOTER_REV,
    InsertElement,
    assemble_insert,
    get_promoter,
    transcribed_segments,
)


def test_t7_canonical_sequence():
    assert T7_PROMOTER == "TAATACGACTCACTATAGGG"  # canonical phi10 +1 GGG start


def test_t7_rev_is_reverse_complement_of_t7():
    assert T7_PROMOTER_REV == str(Seq(T7_PROMOTER).reverse_complement())
    assert T7_PROMOTER_REV == "CCCTATAGTGAGTCGTATTA"


def test_get_promoter_directions():
    assert get_promoter("T7").direction == "forward"
    assert get_promoter("T7_rev").direction == "reverse"
    assert get_promoter("T7").sequence == T7_PROMOTER


def test_get_promoter_unknown_raises():
    with pytest.raises(KeyError):
        get_promoter("SP6")


def test_assemble_insert_concatenates_named_segments():
    segs = [
        ("FLAG", "GACTACAAAGACGATGACGACAAG"),
        ("GGGGS", "GGTGGAGGCGGTTCA"),
        ("T7", T7_PROMOTER),
        ("P2A", "GGAAGCGGAGCT"),
        ("ZeoR", "ATGGCCAAGTTGACC"),
    ]
    out = assemble_insert(segs)
    assert out == "".join(s for _, s in segs)


def test_assemble_insert_rejects_non_acgt():
    with pytest.raises(ValueError, match="ACGT"):
        assemble_insert([("bad", "ACGTN")])


def test_forward_t7_transcribes_downstream():
    # 5'-FLAG-GGGGS-T7-P2A-ZeoR-3' : forward T7 transcribes P2A-ZeoR (downstream)
    segs = [
        ("FLAG", "GACTACAAAGACGATGACGACAAG"),
        ("GGGGS", "GGTGGAGGCGGTTCA"),
        ("T7", T7_PROMOTER),
        ("P2A", "GGAAGCGGAGCT"),
        ("ZeoR", "ATGGCCAAGTTGACC"),
    ]
    result = transcribed_segments(segs, promoter_index=2)
    assert result == ["P2A", "ZeoR"]


def test_reverse_t7_transcribes_upstream():
    # 5'-FLAG-GGGGS-T7rev-P2A-ZeoR-3' : T7-rev transcribes GGGGS-FLAG (upstream)
    segs = [
        ("FLAG", "GACTACAAAGACGATGACGACAAG"),
        ("GGGGS", "GGTGGAGGCGGTTCA"),
        ("T7_rev", T7_PROMOTER_REV),
        ("P2A", "GGAAGCGGAGCT"),
        ("ZeoR", "ATGGCCAAGTTGACC"),
    ]
    result = transcribed_segments(segs, promoter_index=2)
    # upstream segments, reported 5'->3' from the promoter outward
    assert result == ["GGGGS", "FLAG"]


def test_insert_element_is_frozen():
    el = InsertElement(name="T7", sequence=T7_PROMOTER, direction="forward")
    with pytest.raises(Exception):
        el.sequence = "X"  # frozen dataclass
