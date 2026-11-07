"""Executable self-check of a Prime Assembly design's strand arithmetic (spec 4.3, ``docs/STRUCTURE.md``).

The two pegRNAs and the donor are derived with reverse-complement arithmetic that is easy to get silently wrong
(ADR 0004 records how it went wrong in tj-PE). So, like ``simulate_tprt``, this works from what would be ORDERED
and reproduces the molecular mechanism instead of re-running the formulas that built it:

* A guide nicks the PAM strand. The 5' side of the nicked strand ends in a free 3'-OH, the **primer**.
* The PBS must anneal to that primer, antiparallel: the PBS is the reverse complement of the primer.
* The reverse transcriptase extends the primer by copying the RTT from its 3' end toward its 5' end (away from the
  PBS), so the new DNA is the **reverse complement of the RTT**.
* Guide A is on the top strand, so its new DNA already reads in the sense orientation. Guide B is on the bottom
  strand, which reads right to left in sense coordinates, so its new DNA is the reverse complement of what it
  writes once it is read on the sense strand.
* Each flap must equal the end of the donor it is meant to anneal to; then the donor stands where the span between
  the nicks was.

It is a self-consistency check, not an experimental proof of the strand logic (ADR 0004 says the same of
``simulate_tprt``): it proves that what is ordered does what the design says, given the mechanism above.
"""

from __future__ import annotations

from typing import Optional, Tuple

from Bio.Seq import Seq

from pegrna_pipeline.hpc.pa_geometry import SPACER_LEN, Guide
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

__all__ = ["PAMismatch", "split_pegrna", "simulate_side", "simulate_pa"]


class PAMismatch(RuntimeError):
    """What is ordered does not do what the Prime Assembly design says it does."""


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def split_pegrna(
    peg: str,
    *,
    side: str,
    flap_len: int,
    pbs_len: int,
    tail_len: int = 3,
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> Tuple[str, str, str, str]:
    """Cut an ordered pegRNA into ``(spacer, rtt, pbs, tail)`` by the stated lengths, checking the scaffold."""
    expected = SPACER_LEN + len(scaffold) + flap_len + pbs_len + tail_len
    if len(peg) != expected:
        raise PAMismatch(
            f"pegRNA {side} has length {len(peg)}, expected {expected} "
            f"(spacer {SPACER_LEN} + scaffold {len(scaffold)} + RTT {flap_len} + PBS {pbs_len} "
            f"+ tail {tail_len})"
        )
    i = SPACER_LEN
    spacer, seen_scaffold = peg[:i], peg[i : i + len(scaffold)]
    if seen_scaffold != scaffold:
        raise PAMismatch(f"pegRNA {side}: the scaffold is not the expected one")
    i += len(scaffold)
    rtt, pbs, tail = (
        peg[i : i + flap_len],
        peg[i + flap_len : i + flap_len + pbs_len],
        peg[-tail_len:],
    )
    if set(tail) != {"U"}:
        raise PAMismatch(f"pegRNA {side}: its 3' tail must be {tail_len} uridines; it is {tail!r}")
    return spacer, rtt, pbs, tail


def _check_guide(genome: str, guide: Guide, side: str) -> None:
    """The guide object must be what the genome really holds, on the strand its side requires."""
    want = "top" if side == "A" else "bottom"
    if guide.strand != want:
        raise PAMismatch(
            f"guide {side} must be on the {want} strand; it is on the {guide.strand} strand"
        )
    p0, p1 = guide.protospacer_span
    q0, q1 = guide.pam_span
    if p1 - p0 != SPACER_LEN or q1 - q0 != 3:
        raise PAMismatch(f"guide {side}: the protospacer and PAM spans have the wrong length")
    if want == "top":
        # protospacer then PAM on the top strand; the nick is 3 bp 5' of the PAM
        adjacent, spacer, pam, cut = q0 == p1, genome[p0:p1], genome[q0:q1], p1 - 3
    else:
        # the guide reads on the bottom strand: both appear reverse-complemented in sense coordinates,
        # with the PAM just left of the protospacer
        adjacent, spacer, pam, cut = q1 == p0, _rc(genome[p0:p1]), _rc(genome[q0:q1]), p0 + 3
    if not adjacent:
        raise PAMismatch(f"guide {side}: the PAM is not next to the protospacer")
    if spacer != guide.spacer:
        raise PAMismatch(
            f"guide {side}: its spacer is not what the genome holds at its protospacer"
        )
    if pam != guide.pam or pam[1:] != "GG":
        raise PAMismatch(
            f"guide {side}: its PAM is {guide.pam!r} but the genome holds {pam!r} (needs NGG)"
        )
    if guide.cut != cut:
        raise PAMismatch(
            f"guide {side}: its cut ({guide.cut}) is not the first base right of the nick, 3 bp from "
            f"its PAM ({cut})"
        )


def simulate_side(
    *,
    side: str,
    genome: str,
    guide: Guide,
    donor: str,
    peg: str,
    flap_len: int,
    pbs_len: int,
    tail_len: int = 3,
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    """Run one pegRNA through the mechanism and return the flap it writes, in SENSE orientation.

    Raises :class:`PAMismatch` unless the guide matches the genome, the spacer is the guide's, the PBS anneals
    to the primer at the nick, and the flap the RT writes equals the donor end it must anneal to.
    """
    if side not in ("A", "B"):
        raise ValueError(f"side must be 'A' or 'B'; got {side!r}")
    _check_guide(genome, guide, side)
    spacer, rtt, pbs, _tail = split_pegrna(
        peg, side=side, flap_len=flap_len, pbs_len=pbs_len, tail_len=tail_len, scaffold=scaffold
    )
    if spacer != guide.spacer:
        raise PAMismatch(f"pegRNA {side}: its spacer is not guide {side}'s spacer")
    if flap_len > len(donor):
        raise PAMismatch(f"flap {side} ({flap_len} nt) is longer than the donor ({len(donor)} nt)")

    cut = guide.cut
    if side == "A":
        # the top strand 5' of the nick; its 3' end is the primer
        nicked = genome[:cut]
    else:
        # the bottom strand 5' of the nick, written 5'->3'; its 3' end is the primer
        nicked = _rc(genome[cut:])
    if len(nicked) < pbs_len:
        raise PAMismatch(f"PBS {side} ({pbs_len} nt) is longer than the DNA 5' of the nick")
    primer = nicked[-pbs_len:]
    if _rc(pbs) != primer:
        raise PAMismatch(
            f"PBS {side} does not anneal to the 3' end of the nicked strand at the cut "
            f"(it should be the reverse complement of {primer})"
        )
    # the RT copies the RTT from its 3' end toward its 5' end, so the new DNA is the reverse complement of the RTT
    written = _rc(rtt)
    # guide B extends the bottom strand: read on the sense strand, its new DNA is the reverse complement
    flap_sense = written if side == "A" else _rc(written)
    end = donor[:flap_len] if side == "A" else donor[-flap_len:]
    if flap_sense != end:
        raise PAMismatch(
            f"flap {side} does not match the {'5' if side == 'A' else '3'}' end of the donor: the RTT "
            f"writes {flap_sense} where the donor has {end}"
        )
    return flap_sense


def simulate_pa(
    *,
    genome: str,
    a: Guide,
    b: Guide,
    donor: str,
    peg_a: str,
    peg_b: str,
    flap_a: int,
    pbs_a: int,
    flap_b: int,
    pbs_b: int,
    tail_len: int = 3,
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
    expected_product: Optional[str] = None,
) -> str:
    """Rebuild the product from the two ORDERED pegRNAs and the donor, and return it.

    The product is the top strand left of nick A, the flap A writes, the donor's middle, the flap B writes
    (read on the sense strand), and the bottom strand right of nick B. It must equal the donor standing where
    the span between the nicks was, and ``expected_product`` when given. Raises :class:`PAMismatch` on any
    inconsistency.
    """
    _check_guide(genome, a, "A")
    _check_guide(genome, b, "B")
    n1, n2 = a.cut, b.cut
    if not n1 < n2:
        raise PAMismatch(
            f"guide A (cut {n1}) must lie left of guide B (cut {n2}): the flaps grow toward each other"
        )
    if flap_a + flap_b > len(donor):
        raise PAMismatch(
            f"the two flaps ({flap_a} + {flap_b} nt) are longer than the donor ({len(donor)} nt): "
            "they would overlap each other"
        )
    written_a = simulate_side(
        side="A",
        genome=genome,
        guide=a,
        donor=donor,
        peg=peg_a,
        flap_len=flap_a,
        pbs_len=pbs_a,
        tail_len=tail_len,
        scaffold=scaffold,
    )
    written_b = simulate_side(
        side="B",
        genome=genome,
        guide=b,
        donor=donor,
        peg=peg_b,
        flap_len=flap_b,
        pbs_len=pbs_b,
        tail_len=tail_len,
        scaffold=scaffold,
    )
    # Everything left of nick A and right of nick B is untouched; the new DNA is the two flaps around the
    # donor's middle, which the donor's own strands supply.
    assembled = (
        genome[:n1] + written_a + donor[flap_a : len(donor) - flap_b] + written_b + genome[n2:]
    )
    intended = genome[:n1] + donor + genome[n2:]
    if (
        assembled != intended
    ):  # unreachable if the flap checks passed; kept as a self-check of this function
        raise PAMismatch(
            "the assembled product differs from the donor standing in place of the span"
        )
    if expected_product is not None and assembled != expected_product:
        raise PAMismatch("the assembled product differs from the expected product")
    return assembled
