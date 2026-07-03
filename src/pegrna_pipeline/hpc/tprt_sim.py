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

__all__ = ["simulate_tprt", "simulate_rcpbs2_capture", "TPRTMismatch"]


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


def simulate_rcpbs2_capture(
    *,
    rcpbs2: str,
    ts_capture_window: str,
    expected_pbs2: Optional[str] = None,
) -> str:
    """Reconstruct the flap-borne PBS2 from the tj-pegRNA rcPBS2 and assert capture.

    Mechanism (Zheng 2023; ADR 0014): the tj-pegRNA carries ``rcPBS2 = RC(TS)``. The
    RT copies it into the 3' end of the newly-synthesised flap, where it appears as
    its reverse complement, ``PBS2 = RC(rcPBS2)``. That flap-borne PBS2 must equal the
    **target-strand** window (5'->3') it captures at the ngRNA nick — i.e. it recovers
    the TS window when ``rcpbs2`` was correctly oriented as ``RC(TS)``.

    Physics: the flap is synthesised as the reverse complement of the pegRNA 3'
    extension template, so a template segment ``rcPBS2 = RC(TS)`` becomes
    ``RC(RC(TS)) = TS`` on the flap. If a non-RC'd rcPBS2 (the TS literal — the old
    ``orient_pbs2``=identity-for-Rv bug) is supplied, ``RC(rcpbs2) != TS`` and
    :class:`TPRTMismatch` is raised.     ``expected_pbs2`` optionally pins the
    reconstructed PBS2 for an extra caller-supplied cross-check.

    Locus-dependent caveat: on an RC-palindromic capture window (``RC(TS) == TS``)
    a TS-literal rcPBS2 is byte-identical to ``RC(TS)``, so the orientation error is
    undetectable from sequence alone — this guard cannot flag it there.
    """
    pbs2_on_flap = _rc(rcpbs2)
    if pbs2_on_flap != ts_capture_window:
        raise TPRTMismatch(
            "rcPBS2 does not capture the target-strand window — orientation error: "
            f"RC(rcPBS2)={pbs2_on_flap!r} != TS window {ts_capture_window!r}. The "
            "tj-pegRNA must carry rcPBS2 = RC(TS) (ADR 0014); a TS-literal rcPBS2 fails "
            "here (check reverse_complement vs identity on the rcPBS2)."
        )
    if expected_pbs2 is not None and pbs2_on_flap != expected_pbs2:
        raise TPRTMismatch(f"reconstructed PBS2 {pbs2_on_flap!r} != expected {expected_pbs2!r}")
    return pbs2_on_flap
