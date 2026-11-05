"""Insertable cargo elements: the T7 promoters ``T7-G`` / ``rcT7-G`` (two G's after TATA).

The design pipeline can splice named elements into the cargo. The two here are
directional T7 RNA-polymerase promoters that enable *linear amplification* of RNA
off the integrated cassette. Two correctness requirements shape them:

1. **Reverse orientation = reverse COMPLEMENT, not reverse.** A promoter is
   directional and single-stranded-meaningful. To transcribe the *other* way, you
   insert the **reverse complement** of the promoter into the (sense) strand: the
   functional promoter then reads 5'->3' on the **non-targeted** strand, pointing
   the opposite way, with its directionality preserved. Merely reversing the
   letters would put the promoter 3'->5' and break it.

2. **Transcription-competent (guanine start), and no longer than it has to be.** T7 RNA polymerase
   initiates at a **guanine** and the start region is G-rich: the consensus is
   ``TAATACGACTCACTATAGGGAGA`` and transcription begins at the first ``G`` after ``TATA``
   (Rong, He, McAllister & Durbin 1998, PNAS 95:515). The element is written as the NON-TEMPLATE
   (coding) strand and the RNA copies that strand, so the initiating bases are **G's on the element as
   written**. The element here is ``TAATACGACTCACTATA`` + ``GG``: the promoter, the initiating G and the
   next G. (An earlier version carried ``GGGG``: the ``GGG`` of the consensus plus a padding ``G`` that
   made it 21 nt, a multiple of 3. At the 3' end of a cargo that ``GGGG`` is the reverse complement of
   ``CCCC``, which paired at p >= 0.996 with the rcPBS2's ``GGGG`` in the CTCF cassette and failed the
   folding gate on every candidate: ADR 0028 addendum, ADR 0029.)

   **Consequence: the element is 19 nt, NOT a multiple of 3.** It no longer keeps a reading frame by
   itself. That is harmless at the 3' end of a cargo or after a stop codon (where the CTCF cassette has it)
   and wrong in the middle of an ORF: anything coding after it is read one base off. The builders say so
   (:func:`frame_shift_after`) instead of the length hiding it.

3. **It is never recoded.** A promoter is not protein, and its consensus cannot change without breaking it.
   The recoder treats these elements as locked wherever they occur (:mod:`cargo_spec`).

Elements (both stored 5'->3' DNA, uppercase ACGT):

* ``T7-G``   — forward: ``TAATACGACTCACTATAGG`` (19 nt). Transcribes the
              segments **downstream** (3') of it.
* ``rcT7-G`` — reverse: ``reverse_complement(T7-G)`` = ``CCTATAGTGAGTCGTATTA``
              (19 nt). On the non-targeted strand this reads 5'->3' as ``T7-G``
              (GG start intact), transcribing the segments **upstream** (5').

Example cargo ``5'-FLAG-GGGGS-[T7-G|rcT7-G]-P2A-ZeoR-3'``:
  * ``T7-G``   -> transcribes ``P2A-ZeoR`` (the 3' end)
  * ``rcT7-G`` -> transcribes ``GGGGS-FLAG`` (the 5' end)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

from Bio.Seq import Seq

__all__ = [
    "T7_MIN_PROMOTER",
    "T7_CORE",
    "T7G_PROMOTER",
    "RC_T7G_PROMOTER",
    "InsertElement",
    "PROMOTERS",
    "pad_t7_to_frame",
    "frame_shift_after",
    "get_promoter",
    "assemble_insert",
    "transcribed_segments",
]

#: The minimal T7 class III (phi10) promoter, ending at the +1 guanine: what the recoder searches for.
T7_MIN_PROMOTER = "TAATACGACTCACTATAG"
#: The promoter with its first two transcribed bases (+1 G, +2 G). This is the element.
T7_CORE = "TAATACGACTCACTATAGG"


def pad_t7_to_frame(core: str = T7_CORE) -> str:
    """Append the minimum ``G`` bases so ``core`` is length-divisible-by-3.

    No longer applied to the shipped elements (see the module docstring): kept for callers that want a
    frame-neutral promoter and accept the longer G run. For the 19-nt element this is 2 G's -> 21 nt.
    """
    n = (3 - len(core) % 3) % 3
    return core + "G" * n


#: Forward T7 promoter: the promoter + the +1 G and one more (19 nt, GG start). Not a multiple of 3.
T7G_PROMOTER = T7_CORE
#: Reverse-orientation T7: the REVERSE COMPLEMENT of T7-G (so the functional
#: promoter reads 5'->3' on the non-targeted strand, GG start intact). Starts
#: with the C's that complement the GG. 19 nt.
RC_T7G_PROMOTER = str(Seq(T7G_PROMOTER).reverse_complement())


def frame_shift_after(sequence: str) -> int:
    """How many bases (0, 1 or 2) this element moves the reading frame of anything coding after it."""
    return len(sequence) % 3


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
