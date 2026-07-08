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

``split_rt_genomic_tail`` and ``exclude_rtgdna_tail`` handle the C5 rule: remove only the
genomic tail bases that overlap the separately recoded between-nicks arm, preserving any tail
suffix outside that arm. ``TwoPegPridictReal`` handles C2 submission and parsing. ``TwoPegDesign``
is still a design representation only and is not yet wired into real construct assembly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Literal, Optional, Protocol, Tuple

from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import DEFAULT_SCORE_CELL_LINE, PegRNADesign, _score_attr
from pegrna_pipeline.hpc.strand_orientation import orient_rcpbs2

__all__ = [
    "TwoPegDesign",
    "TwoPegPridict",
    "TwoPegPridictFake",
    "TwoPegPridictReal",
    "combine_two_peg",
    "split_rt_genomic_tail",
    "exclude_rtgdna_tail",
]


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def split_rt_genomic_tail(rt: str, insert_part: str) -> Tuple[str, str]:
    """Split a PRIDICT RT template into ``(insert_part, genomic_tail)`` (C5; ADR 0016).

    PRIDICT2 has no separate "RTgDNA" column (CX3-verified); the RT template ``RTrevcomp``
    is the desired insert end followed by a **genomic-homology tail** (the RT-template portion
    past the insert that primes the flap onto the genome). ``insert_part`` is the intended
    insert end (same frame/orientation as ``rt``). Returns ``(insert_part, rt[len(insert_part):])``.
    Raises :class:`ValueError` if ``rt`` does not start with ``insert_part`` (RT/insert are
    misaligned -- fail loud rather than silently mis-trim).
    """
    if not insert_part:
        raise ValueError("insert_part must be non-empty when splitting a PRIDICT RT template")
    if not rt.startswith(insert_part):
        raise ValueError(
            f"RT template {rt!r} does not start with the intended insert part {insert_part!r} "
            "-- RT/insert misaligned"
        )
    return insert_part, rt[len(insert_part) :]


def exclude_rtgdna_tail(rt: str, insert_part: str, *, overlap_len: Optional[int] = None) -> str:
    """Return the RT end with its genomic-homology tail EXCLUDED (C5; ADR 0016).

    In the recoded tj-PE workflow the between-nicks arm is separately re-templated
    (synonymously recoded, breaking microhomology), so it authoritatively covers that genomic
    span to the nucleotide. Any genomic tail PRIDICT appended to the RT would duplicate -- and
    un-recode -- that span, re-introducing the microhomology the recode exists to break, so it
    is excluded. When ``overlap_len`` is supplied, only that many bases of the tail are removed;
    any non-overlapping tail suffix is preserved. Without ``overlap_len`` the full tail is treated
    as overlapping, preserving the initial fake/CI call shape until integration supplies geometry.
    """
    insert, tail = split_rt_genomic_tail(rt, insert_part)
    if overlap_len is None:
        overlap_len = len(tail)
    if overlap_len < 0 or overlap_len > len(tail):
        raise ValueError(
            f"overlap_len must be between 0 and the RT genomic tail length ({len(tail)}); "
            f"got {overlap_len}"
        )
    return insert + tail[overlap_len:]


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


@dataclass
class TwoPegPridictReal:
    """Real two-pegRNA design source (C2; ADR 0016): submit two editseqs to PRIDICT2 on CX3,
    parse the outputs, and select the top-PE2 pegRNA per end.

    Wraps :func:`pegrna_pipeline.hpc.pridict_submitter.submit_pridict_batches` (which qsubs
    ``pridict2_pegRNA_design.py batch`` and waits) and :func:`...pridict_parser.parse_pegrna_csv`.
    ``design_pair`` takes the two **editseqs** (bracket notation, built upstream via
    :func:`...input_formatter.format_bracket_insertion`) rather than pre-split windows -- that is
    the real design input. Only the submit/parse/select orchestration is unit-tested (with a mock
    ssh); the live PRIDICT run is validated on CX3 in a real job (deferred). ``score_cell_line``
    picks which PRIDICT PE2 column ranks the designs (``"k562"`` by default, or ``"hek293t"``).
    """

    ssh: object
    remote_run_dir: str
    pridict_install_path: str
    conda_env: str = "pridict2"
    score_cell_line: str = DEFAULT_SCORE_CELL_LINE

    #: A PRIDICT ``sequence_name`` is interpolated into a remote path
    #: (``<output_dir>/<name>_pegRNA_Pridict_full.csv``) and a batch-CSV cell, so it is
    #: restricted to a filesystem-safe alphabet: no ``/``, ``..``, whitespace, or shell metachars.
    _NAME_RE = re.compile(r"^[A-Za-z0-9_-]+$")

    def __post_init__(self) -> None:
        # Fail fast on an invalid score column: otherwise the ValueError only surfaces in
        # ``_best_for`` -- AFTER a multi-hour CX3 batch has already run (reviewers flagged this).
        _score_attr(self.score_cell_line)

    def design_pair(
        self,
        *,
        editseq1: str,
        editseq2: str,
        name1: str = "peg1_end",
        name2: str = "peg2_end",
        **submit_kwargs: object,
    ) -> Tuple[PegRNADesign, PegRNADesign]:
        # Deferred import to avoid a module-load cycle and keep the fake path import-light.
        from pegrna_pipeline.hpc.pridict_submitter import submit_pridict_batches

        # Validate BEFORE submitting the (expensive, real) CX3 batch. Re-check score_cell_line
        # here too: the dataclass is mutable, so a reassignment could slip past __post_init__.
        _score_attr(self.score_cell_line)
        self._validate_names(name1, name2)

        fixed_submit_fields = {
            "ssh",
            "variants",
            "remote_run_dir",
            "pridict_install_path",
            "conda_env",
        }
        clobbered = sorted(fixed_submit_fields.intersection(submit_kwargs))
        if clobbered:
            raise ValueError(
                "submit_kwargs cannot override fixed PRIDICT submit fields: " + ", ".join(clobbered)
            )

        variants = [
            {"sequence_name": name1, "editseq": editseq1},
            {"sequence_name": name2, "editseq": editseq2},
        ]
        result = submit_pridict_batches(
            ssh=self.ssh,
            variants=variants,
            remote_run_dir=self.remote_run_dir,
            pridict_install_path=self.pridict_install_path,
            conda_env=self.conda_env,
            **submit_kwargs,
        )
        if result.failed:
            raise RuntimeError(f"PRIDICT batch failed for {len(result.failed)} chunk(s)")
        return (
            self._best_for(result.output_dirs, name1),
            self._best_for(result.output_dirs, name2),
        )

    def _validate_names(self, name1: str, name2: str) -> None:
        for nm in (name1, name2):
            if not self._NAME_RE.match(nm):
                raise ValueError(
                    f"pegRNA name {nm!r} must match [A-Za-z0-9_-]+ (it is interpolated into a "
                    "remote path and a batch CSV cell)"
                )
        if name1 == name2:
            # Equal names collide on one <name>_pegRNA_Pridict_full.csv, so both ends would parse
            # the SAME file and return identical designs -- a silently wrong pair. Fail loud.
            raise ValueError(
                f"name1 and name2 must differ; both are {name1!r} (they would collide on a single "
                "PRIDICT output CSV and yield an identical pegRNA pair)"
            )

    def _best_for(self, output_dirs: List[str], name: str) -> PegRNADesign:
        designs: List[PegRNADesign] = []
        for out_dir in output_dirs:
            designs.extend(self._load_designs(out_dir, name))
        if not designs:
            raise RuntimeError(f"no PRIDICT pegRNA designs found for {name!r}")
        attr = _score_attr(self.score_cell_line)
        return max(designs, key=lambda d: getattr(d, attr))

    def _load_designs(self, output_dir: str, name: str) -> List[PegRNADesign]:
        """Download ``<output_dir>/<name>_pegRNA_Pridict_full.csv`` and parse it.

        Thin I/O layer over ``ssh.download`` + :func:`parse_pegrna_csv`; validated on CX3 (a
        missing CSV yields an empty list so :meth:`_best_for` fails loud). Unit tests monkeypatch
        this method to exercise the submit/select orchestration without real I/O.
        """
        import os
        import tempfile

        from pegrna_pipeline.hpc.pridict_parser import parse_pegrna_csv

        remote_csv = f"{output_dir}/{name}_pegRNA_Pridict_full.csv"
        with tempfile.TemporaryDirectory() as tmp:
            local = os.path.join(tmp, f"{name}_pegRNA_Pridict_full.csv")
            try:
                self.ssh.download(remote_csv, local)
            except FileNotFoundError:
                return []
            return parse_pegrna_csv(local)
