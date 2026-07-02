"""Strand orientation for the manually-designed pegRNA extension elements.

CONVENTION (LOCKED 2026-06-19): every sequence is stored **5'->3'**. Antiparallel
/ template-directed relationships are therefore ``reverse_complement`` in code.
The "just complement, not reverse" description is the same physical molecule
viewed in the 3'->5' synthesis frame; serialised back to a 5'->3' oligo it is a
reverse complement.

TARGET-STRAND (TS) CONVENTION (ADR 0014): the whole pegRNA 3'-extension is the
reverse complement of the target strand (the strand the tj-pegRNA targets). Two
manually-designed elements need explicit orientation (PRIDICT-designed PBS1 +
RTT-seed are already correct and are never re-complemented here):

* **RTT tail** (insert-after-FLAG + deletion-replacement-arm): the 3' extension
  is reverse-transcribed into the new flap. For a reverse-targeting pegRNA the
  oligo carries the **reverse complement** of the desired sense flap; for a
  forward-targeting pegRNA it carries the sense flap literally.
* **rcPBS2**: PBS2 does NOT sit on the pegRNA — it ends up at the 3' end of the
  newly-synthesised flap and captures the flap the ngRNA generates on the
  non-target strand. The tj-pegRNA therefore carries the reverse complement,
  ``rcPBS2 = RC(PBS2) = RC(TS)`` (Zheng 2023). :func:`orient_rcpbs2` takes the
  target-strand window at the capture locus and returns its reverse complement.

The old strand-dependent :func:`orient_pbs2` emitted the TS-literal window for a
reverse-targeting pegRNA (i.e. PBS2, not rcPBS2) — a bug. It is now a deprecated
alias of :func:`orient_rcpbs2`.

``simulate_tprt`` (tprt_sim.py) is the source of truth: it models the full
reaction and asserts the integrated **sense** allele equals reference+cargo, so
any misuse of these helpers fails loudly.
"""

from __future__ import annotations

import warnings
from typing import Optional

from Bio.Seq import Seq

__all__ = ["orient_rtt_tail", "orient_rcpbs2", "orient_pbs2"]

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


def orient_rcpbs2(ts_window: str) -> str:
    """Return rcPBS2 as carried by the tj-pegRNA (ADR 0014).

    ``ts_window`` is the target-strand (TS) segment, 5'->3', at the capture locus
    immediately downstream of the ngRNA nick. The tj-pegRNA carries the reverse
    complement of that window: ``rcPBS2 = RC(TS)`` — because PBS2 physically lives
    on the newly-synthesised flap (it captures the flap the ngRNA exposes on the
    non-target strand), and the pegRNA-borne element is its reverse complement
    (Zheng 2023). TS-centric: there is no strand branch here; the caller supplies
    the target-strand window.
    """
    return _rc(ts_window)


def orient_pbs2(ts_window: str, *, pegrna_strand: Optional[str] = None) -> str:
    """Deprecated alias for :func:`orient_rcpbs2` (ADR 0014).

    The old two-branch, strand-dependent ``orient_pbs2`` returned the TS-literal
    window for a reverse-targeting pegRNA — i.e. PBS2, not rcPBS2 — which was a
    bug (the tj-pegRNA carries ``rcPBS2 = RC(TS)``). This shim emits a
    :class:`DeprecationWarning` and returns ``orient_rcpbs2(ts_window)``. The
    ``pegrna_strand`` argument is accepted for backward compatibility and ignored.
    """
    warnings.warn(
        "orient_pbs2 is deprecated; use orient_rcpbs2(ts_window). BEHAVIOUR CHANGED: "
        "it now returns rcPBS2 = RC(TS) for every strand (the old Rv branch returned "
        "the TS-literal window, which was a bug — ADR 0014); the pegrna_strand argument "
        "is accepted for backward compatibility and ignored.",
        DeprecationWarning,
        stacklevel=2,
    )
    return orient_rcpbs2(ts_window)
