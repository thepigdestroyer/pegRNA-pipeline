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

from pegrna_pipeline.hpc.codon_expander import codon_frequency, codons_for_aa
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
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0/1/2; got {frame}")
    # ``frame`` is the reading-frame offset of cds_span[0]. Only whole codons are
    # eligible to change; the partial-codon prefix (pad bases) and any trailing
    # partial codon stay byte-exact so we NEVER introduce a non-synonymous change.
    pad = (3 - frame) % 3
    prefix = cds_span[:pad]
    body = cds_span[pad:]
    n_codons = len(body) // 3
    codons = [body[i * 3 : i * 3 + 3] for i in range(n_codons)]
    trailing = body[n_codons * 3 :]
    best: List[str] = list(codons)

    def assembled(cs: List[str]) -> str:
        return prefix + "".join(cs) + trailing

    # Greedy: while the longest run is too long, pick the single codon swap that
    # most reduces the run (ties broken by fewest base changes), repeat. Bounded
    # by the number of codons to guarantee termination.
    for _ in range(n_codons + 1):
        cur = assembled(best)
        if longest_common_run(cur, genome) <= max_run:
            return RecodeResult(cur, _edit_distance(cur, cds_span), longest_common_run(cur, genome))
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
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0/1/2; got {frame}")
    # Recode only the codon-aligned interior; keep partial-codon edges byte-exact.
    pad = (3 - frame) % 3
    if pad or (len(cds_span) - pad) % 3:
        interior_len = ((len(cds_span) - pad) // 3) * 3
        prefix = cds_span[:pad]
        interior = cds_span[pad : pad + interior_len]
        suffix = cds_span[pad + interior_len :]
        inner = recode_min_edits_dnachisel(interior, genome=genome, frame=0, max_run=max_run)
        if inner is None:
            return None
        recoded = prefix + inner.recoded + suffix
        run = longest_common_run(recoded, genome)
        if run > max_run:
            return None
        return RecodeResult(recoded, _edit_distance(recoded, cds_span), run)
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


def _prefer_human_codons(
    recoded: str, *, original: str, genome: str, frame: int, max_run: int
) -> str:
    """Upgrade CHANGED codons to the highest-frequency human synonym (ADR 0009).

    A pure within-min-edits tiebreak: for every codon position where ``recoded``
    already differs from ``original`` (a base the recode was forced to change to
    break microhomology), swap in the most frequent human synonymous codon that
    still keeps ``longest_common_run(., genome) <= max_run``. Positions equal to
    ``original`` are left byte-exact, so this never adds edits and a no-op recode
    stays a no-op. Iterated to a fixpoint (each swap strictly raises that codon's
    frequency, so it terminates) so the result is locally frequency-optimal at the
    changed positions given the run constraint. ``frame`` is the reading-frame
    offset of ``recoded[0]``; the partial-codon prefix/suffix is never touched.
    """
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0/1/2; got {frame}")
    pad = (3 - frame) % 3
    n_codons = (len(recoded) - pad) // 3
    seq = recoded.upper()
    orig = original.upper()
    gen = genome.upper()
    for _ in range(n_codons + 1):
        changed_any = False
        for i in range(n_codons):
            s = pad + 3 * i
            codon = seq[s : s + 3]
            if codon == orig[s : s + 3]:
                continue  # unchanged position: keep byte-exact (min-edits / no-op)
            aa = str(Seq(codon).translate())
            best, best_f = codon, codon_frequency(aa, codon)
            for alt in codons_for_aa(aa):
                af = codon_frequency(aa, alt)
                if af <= best_f:
                    continue
                trial = seq[:s] + alt + seq[s + 3 :]
                if longest_common_run(trial, gen) <= max_run:
                    best, best_f = alt, af
            if best != codon:
                seq = seq[:s] + best + seq[s + 3 :]
                changed_any = True
        if not changed_any:
            break
    return seq


def recode_arm(cds_span: str, *, genome: str, frame: int, max_run: int) -> Optional[RecodeResult]:
    """Dispatch to the DnaChisel optimiser when importable, else the pure greedy,
    then apply the human-frequency codon tiebreak (:func:`_prefer_human_codons`).

    Both backends honour the same contract: protein-identical, longest common run
    with ``genome`` ``<= max_run``, minimised edits, or ``None`` if infeasible. On
    top of that, ``recode_arm`` prefers the most frequent human synonymous codon at
    each position the recode changed (ADR 0009) -- so the forced microhomology-break
    edits land on human-preferred codons -- without adding edits or disturbing the
    byte-exact unchanged positions.
    """
    try:
        import dnachisel  # noqa: F401, PLC0415
    except Exception:
        result = recode_min_edits(cds_span, genome=genome, frame=frame, max_run=max_run)
    else:
        result = recode_min_edits_dnachisel(cds_span, genome=genome, frame=frame, max_run=max_run)
        if result is None:
            # Fall back to greedy if the optimiser failed to find a solution.
            result = recode_min_edits(cds_span, genome=genome, frame=frame, max_run=max_run)
    if result is None:
        return None
    tuned = _prefer_human_codons(
        result.recoded, original=cds_span, genome=genome, frame=frame, max_run=max_run
    )
    if tuned == result.recoded:
        return result
    return RecodeResult(tuned, _edit_distance(tuned, cds_span), longest_common_run(tuned, genome))
