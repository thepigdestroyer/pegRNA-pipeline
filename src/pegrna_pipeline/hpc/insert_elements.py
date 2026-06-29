"""Insertable cargo elements, including the T7 / T7-rev promoters.

The design pipeline can splice named elements into the cargo. The two added here
are directional T7 RNA-polymerase promoters that enable *linear amplification*
of RNA off the integrated cassette:

* ``T7``      — canonical phi10 promoter (5'->3'), transcribes the segments
                **downstream** (3') of it.
* ``T7_rev``  — the reverse complement, transcribes the segments **upstream**
                (5') of it (the promoter is directional, so flipping it flips the
                direction of transcription).

Example cargo ``5'-FLAG-GGGGS-[T7|T7_rev]-P2A-ZeoR-3'``:
  * ``T7``     -> transcribes ``P2A-ZeoR`` (the 3' end)
  * ``T7_rev`` -> transcribes ``GGGGS-FLAG`` (the 5' end)

All element sequences are stored 5'->3' DNA (uppercase ACGT).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

from Bio.Seq import Seq

__all__ = [
    "T7_PROMOTER",
    "T7_PROMOTER_REV",
    "InsertElement",
    "PROMOTERS",
    "get_promoter",
    "assemble_insert",
    "transcribed_segments",
]

#: Canonical T7 class III phi10 promoter with the +1 GGG transcription start.
T7_PROMOTER = "TAATACGACTCACTATAGGG"
#: T7 pointing the opposite way (directional) — reverse complement of T7.
T7_PROMOTER_REV = str(Seq(T7_PROMOTER).reverse_complement())


@dataclass(frozen=True)
class InsertElement:
    name: str
    sequence: str
    #: "forward" transcribes downstream (3'); "reverse" transcribes upstream (5').
    direction: str


PROMOTERS = {
    "T7": InsertElement(name="T7", sequence=T7_PROMOTER, direction="forward"),
    "T7_rev": InsertElement(name="T7_rev", sequence=T7_PROMOTER_REV, direction="reverse"),
}


def get_promoter(name: str) -> InsertElement:
    """Return the registered promoter element, or raise ``KeyError``."""
    if name not in PROMOTERS:
        raise KeyError(f"unknown promoter {name!r}; known: {sorted(PROMOTERS)}")
    return PROMOTERS[name]


def _check_acgt(label: str, seq: str) -> None:
    if not seq or set(seq.upper()) - set("ACGT"):
        raise ValueError(f"{label} must be non-empty uppercase ACGT; got {seq!r}")


def assemble_insert(segments: Sequence[Tuple[str, str]]) -> str:
    """Concatenate ordered ``(name, dna)`` segments into the full cargo (5'->3').

    Every segment's DNA must be uppercase ACGT (promoter elements included).
    """
    parts: List[str] = []
    for name, dna in segments:
        _check_acgt(name, dna)
        parts.append(dna.upper())
    return "".join(parts)


def transcribed_segments(segments: Sequence[Tuple[str, str]], *, promoter_index: int) -> List[str]:
    """Return the names of segments transcribed by the promoter at
    ``promoter_index``, in 5'->3' transcription order from the promoter outward.

    A ``forward`` promoter transcribes the downstream segments (indices >
    promoter_index) in their natural order. A ``reverse`` promoter transcribes
    the upstream segments (indices < promoter_index) in reverse order (nearest
    first), reflecting transcription proceeding 5'-ward.
    """
    name = segments[promoter_index][0]
    element = get_promoter(name)
    if element.direction == "forward":
        return [n for n, _ in segments[promoter_index + 1 :]]
    upstream = [n for n, _ in segments[:promoter_index]]
    return list(reversed(upstream))
