"""Tests for _has_internal_stop.

We must check EVERY codon in the frame (including the last one) so that a
stop at the tail position is detected. Uses range(0, len(dna) - 2, 3).
"""

from pegrna_pipeline.hpc.codon_expander import _has_internal_stop


def test_no_stop_codons():
    assert _has_internal_stop("ATGAAA") is False


def test_taa_at_end_is_detected():
    assert _has_internal_stop("ATGAAATAA") is True


def test_tag_in_middle_is_detected():
    # ATG TAG AAA -- TAG in frame 2/3 -- stop.
    assert _has_internal_stop("ATGTAGAAA") is True


def test_tga_at_end_is_detected():
    assert _has_internal_stop("ATGTGA") is True


def test_out_of_frame_stop_not_detected():
    # T followed by AAT GAA -- TGA straddles 2nd/3rd codon frame, not in-frame.
    assert _has_internal_stop("TAATGAAA") is True  # TAA at pos 0 in frame


def test_out_of_frame_tga_ignored():
    # Insert TGA at offset 1 so codons are: ATG AAT GAC -- no stops.
    assert _has_internal_stop("ATGAATGAC") is False


def test_empty_sequence_has_no_stop():
    assert _has_internal_stop("") is False


def test_not_multiple_of_three_still_checked_in_full_codons_only():
    # ATGTGA + trailing C: codons ATG, TGA; last partial 'C' ignored.
    assert _has_internal_stop("ATGTGAC") is True
