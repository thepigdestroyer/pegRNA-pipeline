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

__all__ = ["recode_min_edits", "recode_min_edits_dnachisel", "recode_arm", "RecodeResult"]


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


def recode_min_edits_dnachisel(
    cds_span: str, *, genome: str, frame: int, max_run: int
) -> Optional[RecodeResult]:
    """Constrained-optimisation recode using DnaChisel (optional dep ``[recode]``).

    Constraints: ``EnforceTranslation`` (protein identity) + an ``AvoidPattern``
    for every ``(max_run + 1)``-mer of ``genome`` (so no contiguous identity run
    can exceed ``max_run``). Objective: ``AvoidChanges`` (minimise edits).
    Returns a :class:`RecodeResult` or ``None`` if infeasible / not below run.
    """
    if frame != 0:
        raise ValueError("recode_min_edits_dnachisel requires frame=0 (codon-aligned)")
    try:
        from dnachisel import (  # noqa: PLC0415
            AvoidChanges,
            AvoidPattern,
            DnaOptimizationProblem,
            EnforceTranslation,
        )
    except Exception as exc:  # pragma: no cover - exercised only without dnachisel
        raise RuntimeError("DnaChisel not available; install '.[recode]'") from exc

    genome = genome.upper()
    k = max_run + 1
    kmers = {genome[i : i + k] for i in range(len(genome) - max_run)}
    kmers = {km for km in kmers if len(km) == k}
    constraints = [EnforceTranslation()]
    constraints.extend(AvoidPattern(km) for km in kmers)
    try:
        problem = DnaOptimizationProblem(
            sequence=cds_span.upper(),
            constraints=constraints,
            objectives=[AvoidChanges()],
            logger=None,
        )
        problem.resolve_constraints()
        problem.optimize()
    except Exception:
        return None
    recoded = problem.sequence
    run = longest_common_run(recoded, genome)
    if run > max_run:
        return None
    return RecodeResult(recoded, _edit_distance(recoded, cds_span), run)


def recode_arm(
    cds_span: str, *, genome: str, frame: int, max_run: int
) -> Optional[RecodeResult]:
    """Dispatch to the DnaChisel optimiser when importable, else the pure greedy.

    Both backends honour the same contract: protein-identical, longest common run
    with ``genome`` ``<= max_run``, minimised edits, or ``None`` if infeasible.
    """
    try:
        import dnachisel  # noqa: F401, PLC0415
    except Exception:
        return recode_min_edits(cds_span, genome=genome, frame=frame, max_run=max_run)
    result = recode_min_edits_dnachisel(cds_span, genome=genome, frame=frame, max_run=max_run)
    if result is not None:
        return result
    # Fall back to greedy if the optimiser failed to find a solution.
    return recode_min_edits(cds_span, genome=genome, frame=frame, max_run=max_run)
