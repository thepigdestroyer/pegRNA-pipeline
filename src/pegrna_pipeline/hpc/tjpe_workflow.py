"""Tethered-junction prime editing (TJ-PE) workflow.

Entry point: :func:`run_tjpe_workflow`. See Plan 2 v1.1 for the pipeline
contract. Tasks here:

* :func:`rank_variants_by_pe2` — T6
* :func:`select_tjpe_ngrna`   — T7
* :func:`derive_rc_pbs2`      — T8
* :func:`assemble_tj_pegrna`  — T9
* :class:`TJPEResult`         — T10
* :func:`run_tjpe_workflow`   — T10
"""

from __future__ import annotations

import re
from typing import List, Optional

from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import (
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)


OPPOSITE_STRAND = {"Fw": "Rv", "Rv": "Fw"}


#: Standard 76-nt sgRNA scaffold sequence (tracrRNA fusion used in
#: prime-editing constructs).
STANDARD_SGRNA_SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
)


_ACGT_ONLY = re.compile(r"^[ACGT]+$")


# ------------------------------------------------------------------ T6


def rank_variants_by_pe2(variants: List[VariantScore]) -> List[VariantScore]:
    """Return a new list, sorted by ``top_score`` descending, stably.

    Stable sort ensures callers can rely on insertion order for ties.
    """
    return sorted(variants, key=lambda v: v.top_score, reverse=True)


# ------------------------------------------------------------------ T7


def _spacer_passes_sanity(spacer: str) -> bool:
    """Reject unusable spacers on their own sequence content.

    * ``TTTT`` is a RNA-pol III terminator — guide will truncate.
    * A homopolymer run of >= 5 bp kills expression.
    * GC content outside [0.30, 0.75] gives poor DNA binding.
    """
    s = spacer.upper()
    if "TTTT" in s:
        return False
    # Homopolymer >= 5
    run = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            run += 1
            if run >= 5:
                return False
        else:
            run = 1
    if not s:
        return False
    gc = (s.count("G") + s.count("C")) / len(s)
    if gc < 0.30 or gc > 0.75:
        return False
    return True


def select_tjpe_ngrna(
    candidates: List[NgRNADesign],
    *,
    min_distance: int = 40,
    max_distance: int = 100,
    pegrna_strand: Optional[str] = None,
) -> Optional[NgRNADesign]:
    """Pick the best (highest ``deepcas9_score``) surviving ngRNA.

    Parameters
    ----------
    candidates:
        List of :class:`NgRNADesign` rows, typically parsed from a PRIDICT
        nicking-guides CSV.
    min_distance, max_distance:
        Inclusive bounds on ``|distance_to_edit|``.
    pegrna_strand:
        If given (``"Fw"`` or ``"Rv"``), only keep opposite-strand
        candidates: ``{"Fw": "Rv", "Rv": "Fw"}``. Raises :class:`KeyError`
        for any other value — unknown strands are always a bug, not a
        silent filter.

    Returns
    -------
    The candidate with the maximum ``deepcas9_score`` that passes every
    filter, or ``None`` if no candidate survives.
    """
    required_strand: Optional[str] = None
    if pegrna_strand is not None:
        required_strand = OPPOSITE_STRAND[pegrna_strand]  # raises KeyError if unknown

    survivors: List[NgRNADesign] = []
    for c in candidates:
        d = abs(c.distance_to_edit)
        if d < min_distance or d > max_distance:
            continue
        if required_strand is not None and c.target_strand != required_strand:
            continue
        if not _spacer_passes_sanity(c.spacer):
            continue
        survivors.append(c)

    if not survivors:
        return None
    return max(survivors, key=lambda c: c.deepcas9_score)


# ------------------------------------------------------------------ T8


def derive_rc_pbs2(
    top_strand_sequence: str,
    *,
    ngrna_nick_position: int,
    pbs2_length: int = 17,
) -> str:
    """Return the reverse complement of the PBS2 window.

    PBS2 is the primer-binding site used by the **downstream** pegRNA in
    a TJ-PE construct. It's extracted from the top strand immediately
    5' of the ngRNA nick, then reverse-complemented so it anneals to the
    opposite (nicked) strand.

    Parameters
    ----------
    top_strand_sequence:
        mRNA-sense sequence (the same strand as ``GeneRecord.mrna_sequence``).
    ngrna_nick_position:
        0-indexed position on ``top_strand_sequence`` where the ngRNA cuts.
    pbs2_length:
        Number of nt upstream of the nick to use as PBS2 (default 17).

    Raises
    ------
    ValueError
        If the window ``[nick - pbs2_length, nick)`` falls outside the
        sequence.
    """
    start = ngrna_nick_position - pbs2_length
    if start < 0:
        raise ValueError(
            f"PBS2 window start {start} < 0 "
            f"(nick={ngrna_nick_position}, pbs2_length={pbs2_length})"
        )
    if ngrna_nick_position > len(top_strand_sequence):
        raise ValueError(
            f"ngrna_nick_position {ngrna_nick_position} exceeds sequence length "
            f"{len(top_strand_sequence)}"
        )
    window = top_strand_sequence[start:ngrna_nick_position]
    return str(Seq(window).reverse_complement())


# ------------------------------------------------------------------ T9


def _check_acgt(label: str, value: str) -> None:
    if not value or not _ACGT_ONLY.match(value):
        raise ValueError(f"{label} must be non-empty uppercase ACGT; got {value!r}")


def assemble_tj_pegrna(
    *,
    spacer: str,
    insert_dna: str,
    pbs1: str,
    rc_pbs2: str,
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    """Return the full TJ-PE pegRNA DNA string.

    Layout (5' -> 3'):

        spacer + scaffold + insert_dna + pbs1 + rc_pbs2

    All four sequence components (plus the scaffold) must be non-empty
    uppercase ACGT. Any lowercase letter, ``U``, or ambiguity code raises
    :class:`ValueError`.
    """
    for label, value in (
        ("spacer", spacer),
        ("insert_dna", insert_dna),
        ("pbs1", pbs1),
        ("rc_pbs2", rc_pbs2),
        ("scaffold", scaffold),
    ):
        _check_acgt(label, value)
    return spacer + scaffold + insert_dna + pbs1 + rc_pbs2
