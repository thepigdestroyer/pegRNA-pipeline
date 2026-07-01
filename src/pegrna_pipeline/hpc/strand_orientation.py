"""Strand orientation for the manually-designed pegRNA extension elements.

CONVENTION (LOCKED 2026-06-19): every sequence is stored **5'->3'**. Antiparallel
/ template-directed relationships are therefore ``reverse_complement`` in code.
The "just complement, not reverse" description is the same physical molecule
viewed in the 3'->5' synthesis frame; serialised back to a 5'->3' oligo it is a
reverse complement.

Two manually-designed elements need explicit orientation as a function of the
pegRNA target strand (PRIDICT-designed PBS1 + RTT-seed are already correct and
are never re-complemented here):

* **RTT tail** (insert-after-FLAG + deletion-replacement-arm): the 3' extension
  is reverse-transcribed into the new flap. For a reverse-targeting pegRNA the
  oligo carries the **reverse complement** of the desired sense flap; for a
  forward-targeting pegRNA it carries the sense flap literally.
* **PBS2**: it is NOT bound by the pegRNA — it ends up at the 3' end of the new
  flap and must anneal to the strand the nicking-sgRNA exposes. The net result
  (complement-of-complement) is that PBS2-on-the-pegRNA equals the captured
  strand segment. For a reverse-targeting pegRNA the captured strand is the
  sense window taken literally; for a forward-targeting pegRNA it is the reverse
  complement.

``simulate_tprt`` (tprt_sim.py) is the source of truth: it models the full
reaction and asserts the integrated **sense** allele equals reference+cargo, so
any misuse of these helpers fails loudly.
"""

from __future__ import annotations

from Bio.Seq import Seq

__all__ = ["orient_rtt_tail", "orient_pbs2"]

_VALID_STRANDS = ("Fw", "Rv")


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _check_strand(pegrna_strand: str) -> None:
    if pegrna_strand not in _VALID_STRANDS:
        raise ValueError(f"pegrna_strand must be one of {_VALID_STRANDS}; got {pegrna_strand!r}")


def orient_rtt_tail(sense_flap: str, *, pegrna_strand: str) -> str:
    """Return the RTT-tail as it must appear in the 5'->3' pegRNA oligo.

    ``sense_flap`` is the desired product on the sense (top) strand. For a
    reverse-targeting pegRNA the oligo carries its reverse complement; for a
    forward-targeting pegRNA it carries it literally.
    """
    _check_strand(pegrna_strand)
    return _rc(sense_flap) if pegrna_strand == "Rv" else sense_flap


def orient_pbs2(captured_strand_seq: str, *, pegrna_strand: str) -> str:
    """Return PBS2 as it must appear in the 5'->3' pegRNA oligo.

    ``captured_strand_seq`` is the segment (5'->3') of the strand the nicking
    sgRNA exposes that the new flap must capture. For a reverse-targeting pegRNA
    PBS2 equals that segment literally (complement-of-complement = identity); for
    a forward-targeting pegRNA it is the reverse complement.
    """
    _check_strand(pegrna_strand)
    return captured_strand_seq if pegrna_strand == "Rv" else _rc(captured_strand_seq)
