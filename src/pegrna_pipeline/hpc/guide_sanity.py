"""Sequence-content rules for a guide spacer: one definition, for every path that picks a guide.

Two kinds of finding, kept apart on purpose:

HARD rules. A spacer that fails one is unusable and is never chosen.

* ``TTTT``: RNA polymerase III terminates transcription at a run of four T, so a guide
  expressed from a U6 promoter is truncated and does not work.
* GC content outside 30-75%: too low binds weakly, too high mis-folds.

ADVISORY findings. Reported on the design that is chosen, never a reason to reject it.

* a homopolymer run of 5 or more of any base.

Why a homopolymer is advisory and not a rejection. TTTT is special because Pol III reads it
as a stop signal; that is a mechanism, and it is why it is a hard rule. A run of five C, A or G
has no such mechanism, and no primary evidence was found that it prevents a guide from being
made or from working (general design guides advise avoiding runs of C and G, but cite nothing
specific). Treating it as a hard rule had a concrete cost: the same real guide
(``GCTGCCACAGATGCCCCCAA``) was rejected on one candidate path and accepted on the other,
because only one of them applied the check. So it is reported, with its position, and the
chooser decides.

Every site that selects a guide calls :func:`assess_spacer` / :func:`spacer_passes_sanity`.
Nothing re-implements these rules.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

__all__ = [
    "GC_MIN",
    "GC_MAX",
    "HOMOPOLYMER_ADVISORY_RUN",
    "SpacerAssessment",
    "assess_spacer",
    "spacer_passes_sanity",
]

GC_MIN = 0.30
GC_MAX = 0.75
HOMOPOLYMER_ADVISORY_RUN = 5


@dataclass(frozen=True)
class SpacerAssessment:
    spacer: str
    #: HARD-rule violations. Non-empty means the spacer must not be used.
    failures: Tuple[str, ...] = ()
    #: Findings worth showing to the user that do not disqualify the spacer.
    advisories: Tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return not self.failures


def _homopolymer_runs(s: str, min_run: int) -> list:
    """(base, start, length) of every run of >= ``min_run`` identical bases; start is 0-based."""
    runs = []
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j - i + 1 >= min_run:
            runs.append((s[i], i, j - i + 1))
        i = j + 1
    return runs


def assess_spacer(spacer: str) -> SpacerAssessment:
    """Check ``spacer`` against the hard rules and collect advisories."""
    s = (spacer or "").upper()
    failures = []
    advisories = []

    if not s:
        return SpacerAssessment(spacer=s, failures=("empty spacer",))

    if "TTTT" in s:
        failures.append(
            "contains TTTT: RNA polymerase III terminates there, so a U6-expressed guide "
            "is truncated"
        )

    gc = (s.count("G") + s.count("C")) / len(s)
    if gc < GC_MIN or gc > GC_MAX:
        failures.append(
            f"GC content {100 * gc:.0f}% is outside {100 * GC_MIN:.0f}-{100 * GC_MAX:.0f}%"
        )

    for base, start, length in _homopolymer_runs(s, HOMOPOLYMER_ADVISORY_RUN):
        advisories.append(
            f"homopolymer run {base * length} ({length} nt, positions {start + 1}-"
            f"{start + length}): not known to block expression; check activity if the "
            "guide underperforms"
        )

    return SpacerAssessment(spacer=s, failures=tuple(failures), advisories=tuple(advisories))


def spacer_passes_sanity(spacer: str) -> bool:
    """True when ``spacer`` violates no HARD rule (advisories do not count)."""
    return assess_spacer(spacer).ok
