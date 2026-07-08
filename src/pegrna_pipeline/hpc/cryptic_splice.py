"""Cryptic splice-site gate via MaxEntScan motif strength (ADR 0011).

Why this exists (and is separate from :mod:`splice_screen`): SpliceAI is a 10 kb
context-aware model and rightly refuses to run on the short design window the
recoded workflow carries, and it is out-of-distribution for the large exogenous
cassette insertion. The dominant splice risk for a cassette knock-in is a *cryptic*
splice site introduced by the inserted/recoded sequence itself. MaxEntScan (Yeo &
Burge 2004) scores the intrinsic strength of a candidate 5' (donor) or 3' (acceptor)
splice site from a short motif window -- no genomic context needed -- so it is the
right instrument for "did the edit create a strong new splice site?".

This module provides:
  * :class:`MaxEntScorer` -- Protocol (``score5``/``score3``).
  * :class:`MaxEntScanScorer` -- real adapter over ``maxentpy`` (CX3-only; staged
    in the ``[splice]`` extra). Lazy import; raises a helpful error if absent.
  * :class:`MaxEntScanFake` -- deterministic CI/test fake.
  * :func:`passes_cryptic_splice_gate` -- the hard gate: scans the *variant region*
    of the edited construct (the inserted/recoded bases plus the junctions that span
    them) for a donor/acceptor motif at or above ``strong_threshold`` that the edit
    created or strengthened relative to the reference.

MaxEntScan conventions (Yeo & Burge 2004):
  * 5'ss (donor): a 9-mer = 3 exonic + 6 intronic nt, with the ``GT`` at positions
    4-5 (0-indexed 3-4). Strong donor MaxEnt > 8, weak < 5.
  * 3'ss (acceptor): a 23-mer = 20 intronic + 3 exonic nt, with the ``AG`` at the
    end of the intron (0-indexed 18-19). Strong acceptor MaxEnt > 8, weak < 5.
"""

from __future__ import annotations

from typing import List, Optional, Protocol, Tuple

__all__ = [
    "MaxEntScorer",
    "MaxEntScanScorer",
    "MaxEntScanFake",
    "passes_cryptic_splice_gate",
    "STRONG_MAXENT",
]

#: MaxEnt bits at/above which a splice-site motif is "strong" (Yeo & Burge 2004).
STRONG_MAXENT = 8.0

_DONOR_LEN = 9  # 3 exon + 6 intron; GT at 0-indexed 3-4
_ACCEPTOR_LEN = 23  # 20 intron + 3 exon; AG at 0-indexed 18-19


class MaxEntScorer(Protocol):
    def score5(self, nine_mer: str) -> float:
        """MaxEnt score of a 9-nt 5' splice-site (donor) window."""
        ...

    def score3(self, twentythree_mer: str) -> float:
        """MaxEnt score of a 23-nt 3' splice-site (acceptor) window."""
        ...


class MaxEntScanScorer:
    """Real MaxEntScan adapter over ``maxentpy`` (CX3-only; lazy heavy import)."""

    def __init__(self) -> None:
        self._score5 = None
        self._score3 = None
        self._m5 = None
        self._m3 = None

    def _load(self):  # pragma: no cover - exercised only on CX3 with maxentpy staged
        if self._score5 is not None:
            return
        try:
            from maxentpy.maxent import score3 as _s3  # noqa: PLC0415
            from maxentpy.maxent import score5 as _s5  # noqa: PLC0415
        except Exception as exc:  # pragma: no cover
            raise ImportError(
                "MaxEntScanScorer requires the '[splice]' extra (maxentpy) staged on "
                "CX3; run real cryptic-splice scoring there. "
                f"Original error: {exc}"
            ) from exc
        self._score5, self._score3 = _s5, _s3
        # Newer maxentpy loads matrices internally; older needs them passed. Try to
        # load them so both APIs work.
        try:
            from maxentpy.maxent import load_matrix3, load_matrix5  # noqa: PLC0415

            self._m5, self._m3 = load_matrix5(), load_matrix3()
        except Exception:  # pragma: no cover - newer maxentpy has no load_matrix*
            self._m5 = self._m3 = None

    def score5(self, nine_mer: str) -> float:  # pragma: no cover - CX3 only
        self._load()
        return float(self._score5(nine_mer, self._m5) if self._m5 else self._score5(nine_mer))

    def score3(self, twentythree_mer: str) -> float:  # pragma: no cover - CX3 only
        self._load()
        return float(
            self._score3(twentythree_mer, self._m3) if self._m3 else self._score3(twentythree_mer)
        )


class MaxEntScanFake:
    """Deterministic fake for CI/tests.

    ``score5_map`` / ``score3_map`` map an exact motif window to a score; any window
    not present returns ``default`` (0.0 -> below-threshold, i.e. no cryptic site).
    """

    def __init__(
        self,
        *,
        score5_map: dict | None = None,
        score3_map: dict | None = None,
        default: float = 0.0,
    ) -> None:
        self._s5 = {k.upper(): v for k, v in (score5_map or {}).items()}
        self._s3 = {k.upper(): v for k, v in (score3_map or {}).items()}
        self._default = default

    def score5(self, nine_mer: str) -> float:
        return self._s5.get(nine_mer.upper(), self._default)

    def score3(self, twentythree_mer: str) -> float:
        return self._s3.get(twentythree_mer.upper(), self._default)


def _common_prefix(a: str, b: str) -> int:
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _common_suffix(a: str, b: str, cap: int) -> int:
    n = min(len(a), len(b), cap)
    i = 0
    while i < n and a[len(a) - 1 - i] == b[len(b) - 1 - i]:
        i += 1
    return i


def _scan_sites(seq: str, scorer: MaxEntScorer) -> List[Tuple[str, int, int, int, float]]:
    """Return ``(kind, motif_start, motif_end, anchor, score)`` for every ACGT
    donor (``GT``) and acceptor (``AG``) motif window in ``seq``.

    ``anchor`` is the position of the G of GT (donor) / A of AG (acceptor); the
    motif window is ``[motif_start, motif_end)``.
    """
    u = seq.upper()
    n = len(u)
    out: List[Tuple[str, int, int, int, float]] = []
    # Donor: GT at (i, i+1); window seq[i-3 : i+6].
    for i in range(3, n - 5):
        if u[i] == "G" and u[i + 1] == "T":
            start, end = i - 3, i + 6
            nine = u[start:end]
            if set(nine) <= {"A", "C", "G", "T"} and len(nine) == _DONOR_LEN:
                out.append(("donor", start, end, i, scorer.score5(nine)))
    # Acceptor: AG at (k, k+1) as the intron's last 2 nt; window seq[k-18 : k+5].
    for k in range(18, n - 4):
        if u[k] == "A" and u[k + 1] == "G":
            start, end = k - 18, k + 5
            tt = u[start:end]
            if set(tt) <= {"A", "C", "G", "T"} and len(tt) == _ACCEPTOR_LEN:
                out.append(("acceptor", start, end, k, scorer.score3(tt)))
    return out


def passes_cryptic_splice_gate(
    ref_seq: str,
    edited_seq: str,
    *,
    scorer: MaxEntScorer,
    strong_threshold: float = STRONG_MAXENT,
    metrics: Optional[dict] = None,
) -> Tuple[bool, str]:
    """HARD GATE: reject if the edit creates/strengthens a strong cryptic splice site.

    Aligns ``ref_seq`` and ``edited_seq`` on their common prefix/suffix to locate the
    variant region (inserted + recoded bases). Scans the edited construct for donor/
    acceptor motifs whose window overlaps the variant region and scores them with
    MaxEntScan. A motif at or above ``strong_threshold`` is a candidate cryptic site;
    it fails the gate unless an equally strong motif already existed in the reference
    (i.e. the edit did not create or strengthen it). Distant motifs entirely within
    the unchanged flanks are ignored (they are native, not introduced by the edit).
    """
    lr, le = len(ref_seq), len(edited_seq)
    p = _common_prefix(ref_seq, edited_seq)
    s = _common_suffix(ref_seq, edited_seq, cap=min(lr, le) - p)
    var_lo, var_hi = p, le - s  # edited variant region [var_lo, var_hi)

    def _overlaps_variant(start: int, end: int) -> bool:
        return start < var_hi and end > var_lo

    edited_sites = _scan_sites(edited_seq, scorer)
    ref_best = 0.0
    for _kind, _start, _end, _anchor, sc in _scan_sites(ref_seq, scorer):
        if sc > ref_best:
            ref_best = sc

    # Strongest MaxEnt motif the edit introduces in/over the variant region
    # (regardless of threshold) — the raw value surfaced to the progress display.
    var_scores = [sc for (_k, s, e, _a, sc) in edited_sites if _overlaps_variant(s, e)]
    observed_max = max(var_scores) if var_scores else 0.0
    if metrics is not None:
        metrics.update({
            "metric": "MaxEnt", "value": round(observed_max, 1),
            "cutoff": round(strong_threshold, 1), "cmp": "<",
            "extra": f"ref best {ref_best:.1f}",
        })

    created = [
        (kind, start, end, anchor, sc)
        for (kind, start, end, anchor, sc) in edited_sites
        if _overlaps_variant(start, end) and sc >= strong_threshold and sc > ref_best
    ]
    if not created:
        return True, f"no cryptic splice site created (MaxEnt < {strong_threshold:g})"
    kind, _start, _end, anchor, sc = max(created, key=lambda t: t[4])
    return (
        False,
        f"cryptic {kind} splice site MaxEnt {sc:.1f} >= {strong_threshold:g} "
        f"created by edit near pos {anchor}",
    )
