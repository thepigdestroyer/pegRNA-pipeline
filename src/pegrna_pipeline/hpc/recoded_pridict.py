"""Two-pegRNA PRIDICT adapter, orientation, and ngRNA lifecycle (Stage C; ADR 0014/0016).

The recoded tj-pegRNA is assembled from **two** PRIDICT2 designs:

* **pegRNA#1** is designed on the target strand (TS) for the S-end of the insert. PRIDICT
  already emits its PBS and RT as ``PBSrevcomp``/``RTrevcomp`` (i.e. RC(TS)), so its ``PBS1``
  and first insert end are used **untouched**.
* **pegRNA#2** is designed on the **opposite** strand (OS) for the PBS2-end. Because PRIDICT
  emits ``PBSrevcomp = RC(OS) = TS`` (a TS window, literally), reverse-complementing it gives
  ``rcPBS2 = RC(peg2.PBSrevcomp) = RC(TS)``; likewise ``last_end = RC(peg2.RTrevcomp)``. Its
  **spacer + scaffold becomes the ngRNA**, bound to this design.

Every 3'-extension component therefore ends up RC(TS), consistent with the target-strand
convention (ADR 0014). PRIDICT2 exposes no per-spacer off-target-risk field and no separate
"RTgDNA" column (verified on CX3; ADR 0016), so ``offtarget_risk`` defaults to ``None``
(Cas-OFFinder handles spacer QC) and the RT genomic-homology tail lives inside ``RTrevcomp``.

NOTE (C5/C2 pending): ``first_end``/``last_end`` here are the **full** ``RTrevcomp`` (insert end
plus its genomic-homology tail). Excluding the tail's overlap with the recoded between-nicks arm
is C5; the real PRIDICT adapter is C2. ``TwoPegDesign`` is a design representation only and is not
yet wired into real construct assembly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, Protocol, Tuple

from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import PegRNADesign
from pegrna_pipeline.hpc.strand_orientation import orient_rcpbs2

__all__ = ["TwoPegDesign", "TwoPegPridict", "TwoPegPridictFake", "combine_two_peg"]


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


@dataclass(frozen=True)
class TwoPegDesign:
    """A combined tj-pegRNA + bound ngRNA design (components, all RC(TS)).

    ``pbs1`` + ``first_end`` come untouched from pegRNA#1; ``rcpbs2`` + ``last_end`` are the
    reverse-complemented pegRNA#2 outputs; ``ngrna`` = pegRNA#2 spacer+scaffold (bound -- it is
    discarded together with this design). ``pbs1_len``/``rt_len`` and ``editing_position``/
    ``target_strand`` are taken from pegRNA#1 (authoritative for the edit); ``rt_len`` is the
    length of ``first_end`` (pegRNA#1's RT). ``first_end``/``last_end`` are the full RTrevcomp
    (incl. genomic-homology tail) until C5 excludes the tail overlapping the recode arm.
    ``offtarget_risk`` is ``None`` unless supplied (PRIDICT2 has no such field; ADR 0016).
    """

    pbs1: str
    first_end: str
    rcpbs2: str
    last_end: str
    ngrna: str
    pbs1_len: int
    rt_len: int
    editing_position: int
    target_strand: Literal["Fw", "Rv"]
    offtarget_risk: Optional[float] = None


class TwoPegPridict(Protocol):
    """Interface for a two-pegRNA design source (fake in CI, real PRIDICT on CX3).

    Returns ``(pegRNA#1, pegRNA#2)`` as :class:`PegRNADesign` rows shaped exactly like
    :func:`pegrna_pipeline.hpc.pridict_parser.parse_pegrna_csv` output (``.pbs``/``.rt`` are
    ``PBSrevcomp``/``RTrevcomp``, already reverse-complemented by PRIDICT). Orientation into a
    :class:`TwoPegDesign` is applied by :func:`combine_two_peg`, not here.

    Parameters
    ----------
    first_end, last_end, pbs1_window, pbs2_window:
        All **target-strand** (TS) sequences, 5'->3': the desired first/last insert ends and the
        PBS windows on the target strand. An implementation must emit pegRNA#1 designed on TS and
        pegRNA#2 designed on the opposite strand (so ``peg2.PBSrevcomp == RC(OS) == pbs2_window``).
    target_strand:
        ``"Fw"`` or ``"Rv"`` (as in :class:`PegRNADesign`).
    """

    def design_pair(
        self,
        *,
        first_end: str,
        last_end: str,
        pbs1_window: str,
        pbs2_window: str,
        spacer1: str,
        spacer2: str,
        target_strand: Literal["Fw", "Rv"],
        editing_position: int = 0,
    ) -> Tuple[PegRNADesign, PegRNADesign]: ...


class TwoPegPridictFake:
    """Deterministic PRIDICT stand-in for CI: emits PBSrevcomp/RTrevcomp shaped like
    ``parse_pegrna_csv``. All window params are **target-strand** (see :class:`TwoPegPridict`).
    pegRNA#1's fields are RC of the given TS windows (real PRIDICT designs on TS). pegRNA#2's
    fields equal the TS windows **literally** -- because PRIDICT designs peg2 on the opposite
    strand, so ``peg2.PBSrevcomp = RC(OS) = TS`` -- so :func:`combine_two_peg`'s RC then yields
    RC(TS). No PRIDICT run required.
    """

    def design_pair(
        self,
        *,
        first_end: str,
        last_end: str,
        pbs1_window: str,
        pbs2_window: str,
        spacer1: str,
        spacer2: str,
        target_strand: Literal["Fw", "Rv"],
        editing_position: int = 0,
    ) -> Tuple[PegRNADesign, PegRNADesign]:
        peg1 = PegRNADesign(
            sequence_name="peg1",
            spacer=spacer1,
            pbs=_rc(pbs1_window),
            rt=_rc(first_end),
            pbs_length=len(pbs1_window),
            rt_length=len(first_end),
            target_strand=target_strand,
            editing_position=editing_position,
            pe2_hek293t=50.0,
            pe2_k562=50.0,
            pegrna_full="",
        )
        peg2 = PegRNADesign(
            sequence_name="peg2",
            spacer=spacer2,
            pbs=pbs2_window,
            rt=last_end,
            pbs_length=len(pbs2_window),
            rt_length=len(last_end),
            target_strand=target_strand,
            editing_position=editing_position,
            pe2_hek293t=50.0,
            pe2_k562=50.0,
            pegrna_full="",
        )
        return peg1, peg2


def combine_two_peg(
    peg1: PegRNADesign,
    peg2: PegRNADesign,
    *,
    scaffold: str,
    offtarget_risk: Optional[float] = None,
) -> TwoPegDesign:
    """Combine two PRIDICT pegRNA designs into a tj-pegRNA + bound ngRNA (C3; ADR 0014/0016).

    pegRNA#1 (target strand) contributes ``PBS1`` + first insert end **untouched** (already
    RC(TS)). pegRNA#2 (opposite strand) contributes ``rcpbs2 = orient_rcpbs2(peg2.PBSrevcomp)``
    (= RC(TS), the ADR-0014 rcPBS2 helper) and ``last_end = RC(peg2.RTrevcomp)``; its
    ``spacer + scaffold`` becomes the ngRNA, bound to this design (dropped together on discard).
    Edit metadata (editing_position, target_strand, lengths) is taken from pegRNA#1.
    """
    return TwoPegDesign(
        pbs1=peg1.pbs,
        first_end=peg1.rt,
        rcpbs2=orient_rcpbs2(peg2.pbs),
        last_end=_rc(peg2.rt),
        ngrna=peg2.spacer + scaffold,
        pbs1_len=peg1.pbs_length,
        rt_len=peg1.rt_length,
        editing_position=peg1.editing_position,
        target_strand=peg1.target_strand,
        offtarget_risk=offtarget_risk,
    )
