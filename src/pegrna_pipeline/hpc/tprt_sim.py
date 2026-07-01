"""Template-jumping prime editing (TPRT) simulator.

Reconstructs the integrated **sense** allele from the designed pegRNA RTT tail
and the genomic geometry, using only the one bedrock physical fact: a template
and the strand synthesized from it are **reverse complements when both are
written 5'->3'**. The reconstructed sense allele is then compared to the intended
product; any wrong strand operation (plain ``complement`` or plain ``reverse``
instead of ``reverse_complement``) produces a garbled allele and raises
:class:`TPRTMismatch`.

This is the source of truth for strand correctness — orientation helpers
(strand_orientation.py) are validated *through* this simulator, not by
hand-asserted reverse-complement values.

Convention (LOCKED): everything 5'->3'. The simulator shares the strand
convention validated by ``test_scarless_golden.py`` (Fw: the RTT tail carries
the sense flap literally; Rv: it carries the reverse complement). The binding
guarantee is: reconstructed sense allele == intended product.
"""

from __future__ import annotations

from typing import Optional

from Bio.Seq import Seq

__all__ = ["simulate_tprt", "TPRTMismatch"]


class TPRTMismatch(RuntimeError):
    """Raised when the reconstructed sense allele does not equal the intended
    product — i.e. a strand-orientation error in the designed pegRNA."""


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def simulate_tprt(
    *,
    genome: str,
    site_s: int,
    hi: int,
    rtt_tail_on_peg: str,
    pegrna_strand: str,
    expected_allele: Optional[str] = None,
) -> str:
    """Reconstruct the integrated sense allele produced by the designed RTT tail.

    Physics: the new DNA flap synthesized from the RTT-tail template is the
    reverse complement of the RTT tail (both 5'->3'). For a reverse-targeting
    (``Rv``) pegRNA the convention stores the RTT tail as the reverse complement
    of the desired sense flap, so applying the physical reverse-complement
    recovers the sense flap directly. For a forward (``Fw``) pegRNA the RTT tail
    already carries the sense flap, and the physical synthesis frame is undone by
    the same single reverse-complement bookkeeping used to place it.

    The recovered sense flap is spliced into the genome at ``site_s`` (replacing
    nothing endogenous up to ``hi`` — the arm is part of the flap), giving the
    integrated sense allele:

        genome[:site_s] + sense_flap + genome[hi:]

    ``hi`` is the 3' boundary of the re-templated arm on the sense genome (so the
    flap = cargo + restored/recoded arm). If ``expected_allele`` is supplied and
    the reconstruction differs, :class:`TPRTMismatch` is raised; otherwise the
    reconstructed allele is returned.
    """
    if pegrna_strand not in ("Fw", "Rv"):
        raise ValueError(f"pegrna_strand must be 'Fw' or 'Rv'; got {pegrna_strand!r}")

    # Recover the sense flap from the pegRNA RTT tail.
    # Rv: tail was stored as reverse_complement(sense_flap) -> RC recovers sense.
    # Fw: tail is the sense flap literally.
    sense_flap = _rc(rtt_tail_on_peg) if pegrna_strand == "Rv" else rtt_tail_on_peg

    allele = genome[:site_s] + sense_flap + genome[hi:]

    if expected_allele is not None and allele != expected_allele:
        raise TPRTMismatch(
            "reconstructed sense allele != intended product — strand-orientation "
            f"error (got len {len(allele)} / expected len {len(expected_allele)}; "
            "check complement vs reverse vs reverse_complement on the RTT tail)"
        )
    return allele
