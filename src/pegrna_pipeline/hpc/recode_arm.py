"""Minimum-edit synonymous recode of the deletion-replacement arm.

Goal: change the FEWEST bases needed to drop the longest contiguous
perfect-identity run with the displaced genomic strand below ``max_run``, while
keeping every codon synonymous (protein-identical). Pure greedy implementation;
:func:`recode_min_edits_dnachisel` (Task 4) swaps in a constrained optimiser
behind the same signature, and :func:`recode_arm` dispatches to whichever is
available.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from Bio.Seq import Seq

from pegrna_pipeline.hpc.codon_expander import codons_for_aa
from pegrna_pipeline.hpc.microhomology import longest_common_run

__all__ = ["recode_min_edits", "RecodeResult"]


@dataclass(frozen=True)
class RecodeResult:
    recoded: str
    n_edits: int
    max_run_after: int


def _edit_distance(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def recode_min_edits(
    cds_span: str, *, genome: str, frame: int, max_run: int
) -> Optional[RecodeResult]:
    """Return a synonymous recode of ``cds_span`` whose longest common run with
    ``genome`` is ``<= max_run``, using a greedy minimum-edit heuristic, or
    ``None`` if no synonymous choice achieves it.

    ``frame`` is the reading-frame offset of ``cds_span``; only ``frame == 0``
    (codon-aligned) spans are accepted — the caller is responsible for passing a
    span that starts on a codon boundary.
    """
    if frame != 0:
        raise ValueError("recode_min_edits requires a codon-aligned span (frame=0)")

    n_codons = len(cds_span) // 3
    codons = [cds_span[i * 3 : i * 3 + 3] for i in range(n_codons)]
    tail = cds_span[n_codons * 3 :]
    best: List[str] = list(codons)

    def assembled(cs: List[str]) -> str:
        return "".join(cs) + tail

    # Greedy: while the longest run is too long, pick the single codon swap that
    # most reduces the run (ties broken by fewest base changes), repeat. Bounded
    # by the number of codons to guarantee termination.
    for _ in range(n_codons + 1):
        cur = assembled(best)
        if longest_common_run(cur, genome) <= max_run:
            return RecodeResult(cur, _edit_distance(cur, cds_span),
                                longest_common_run(cur, genome))
        cur_run = longest_common_run(cur, genome)
        improved = False
        for ci, codon in enumerate(best):
            aa = str(Seq(codon).translate())
            alts = [c for c in codons_for_aa(aa) if c != codon]
            if not alts:
                continue
            scored = []
            for alt in alts:
                trial = best.copy()
                trial[ci] = alt
                run = longest_common_run(assembled(trial), genome)
                scored.append((run, _edit_distance(alt, codon), ci, alt))
            scored.sort()
            run0, _, _, alt0 = scored[0]
            if run0 < cur_run:
                best[ci] = alt0
                improved = True
                break
        if not improved:
            break

    final = assembled(best)
    run = longest_common_run(final, genome)
    if run <= max_run:
        return RecodeResult(final, _edit_distance(final, cds_span), run)
    return None
