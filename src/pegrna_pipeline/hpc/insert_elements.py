"""Insertable cargo elements: the frame-safe T7 promoters ``T7-G`` / ``rcT7-G``.

The design pipeline can splice named elements into the cargo. The two here are
directional T7 RNA-polymerase promoters that enable *linear amplification* of RNA
off the integrated cassette. Two correctness requirements shape them:

1. **Reverse orientation = reverse COMPLEMENT, not reverse.** A promoter is
   directional and single-stranded-meaningful. To transcribe the *other* way, you
   insert the **reverse complement** of the promoter into the (sense) strand: the
   functional promoter then reads 5'->3' on the **non-targeted** strand, pointing
   the opposite way, with its directionality preserved. Merely reversing the
   letters would put the promoter 3'->5' and break it.

2. **Frame-safe + transcription-competent (guanine start).** T7 RNA polymerase
   initiates at a **guanine** and the start region is G-rich: the consensus is
   ``TAATACGACTCACTATAGGGAGA`` and transcription begins at the first ``G`` of the
   ``GGG`` (Rong, He, McAllister & Durbin 1998, PNAS 95:515). Because the element
   is written as the NON-TEMPLATE (coding) strand and the RNA copies that strand,
   the initiating bases are **G's on the element as written** (they are C's only
   on the opposite/template strand). The canonical core here already carries the
   ``GGG`` start. The core is 20 nt, which is **not** a multiple of 3, so dropped
   mid-cargo it would frameshift the downstream coding sequence (e.g. ZeoR). We
   append the minimum ``G`` bases (1 for the 20-mer) that make the element length
   divisible by 3 — this fixes the frame while keeping the start G-rich. (Padding
   with C would be wrong: C's on the coding strand are not the initiator and do
   not aid transcription.)

Elements (both stored 5'->3' DNA, uppercase ACGT):

* ``T7-G``   — forward: ``TAATACGACTCACTATAGGGG`` (21 nt). Transcribes the
              segments **downstream** (3') of it.
* ``rcT7-G`` — reverse: ``reverse_complement(T7-G)`` = ``CCCCTATAGTGAGTCGTATTA``
              (21 nt). On the non-targeted strand this reads 5'->3' as ``T7-G``
              (GGGG start intact), transcribing the segments **upstream** (5').

Example cargo ``5'-FLAG-GGGGS-[T7-G|rcT7-G]-P2A-ZeoR-3'``:
  * ``T7-G``   -> transcribes ``P2A-ZeoR`` (the 3' end)
  * ``rcT7-G`` -> transcribes ``GGGGS-FLAG`` (the 5' end)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

from Bio.Seq import Seq

__all__ = [
    "T7_CORE",
    "T7G_PROMOTER",
    "RC_T7G_PROMOTER",
    "InsertElement",
    "PROMOTERS",
    "pad_t7_to_frame",
    "get_promoter",
    "assemble_insert",
    "transcribed_segments",
]

#: Canonical T7 class III phi10 promoter with the +1 GGG guanine start.
T7_CORE = "TAATACGACTCACTATAGGG"


def pad_t7_to_frame(core: str = T7_CORE) -> str:
    """Append the minimum ``G`` bases so ``core`` is length-divisible-by-3.

    The G's keep the element from frameshifting a downstream coding insert while
    keeping the start G-rich (T7 initiates at guanine). For the 20-nt core this is
    a single G -> 21 nt. Returns ``core`` unchanged if it is already in-frame.
    """
    n = (3 - len(core) % 3) % 3
    return core + "G" * n


#: Forward, frame-safe T7 promoter: T7 core + frame G's (21 nt, GGGG start).
T7G_PROMOTER = pad_t7_to_frame(T7_CORE)
#: Reverse-orientation T7: the REVERSE COMPLEMENT of T7-G (so the functional
#: promoter reads 5'->3' on the non-targeted strand, GGGG start intact). Starts
#: with the C's that complement the GGGG. 21 nt -> also frame-safe.
RC_T7G_PROMOTER = str(Seq(T7G_PROMOTER).reverse_complement())


@dataclass(frozen=True)
class InsertElement:
    name: str
    sequence: str
    #: "forward" transcribes downstream (3'); "reverse" transcribes upstream (5').
    direction: str


PROMOTERS = {
    "T7-G": InsertElement(name="T7-G", sequence=T7G_PROMOTER, direction="forward"),
    "rcT7-G": InsertElement(name="rcT7-G", sequence=RC_T7G_PROMOTER, direction="reverse"),
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
