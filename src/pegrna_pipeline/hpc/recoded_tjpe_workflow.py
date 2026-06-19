"""Recoded TJ-pegRNA workflow: enumeration, gates, and the ranked-pool-walk.

This module wires together the pure building blocks (microhomology, recode_arm,
nick_search, INV-0) and the injected gate scorers (splice, importance, full
structure, complex docking) into the candidate enumeration and the
ranked-pool-walk orchestrator. Heavy scorers are passed in (fakes in CI, real
adapters on CX3); this module imports none of them at module load.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

from pegrna_pipeline.hpc.codon_expander import enumerate_variants

__all__ = ["PbsRttVariant", "enumerate_pbs_rtt_variants"]


@dataclass(frozen=True)
class PbsRttVariant:
    spacer: str
    rtt_dna: str
    pbs_length: int


def enumerate_pbs_rtt_variants(
    *, spacer: str, rtt_aa: str, pbs_lengths: Tuple[int, ...],
    max_rtt_variants: int = 100,
) -> List[PbsRttVariant]:
    """Enumerate (synonymous RTT) × (PBS1 length); the spacer is FIXED by site S.

    These are the candidates handed to PRIDICT2 (reverse-targeting) for ranking.
    The spacer is never varied — there is one (or a few) protospacer(s) that nick
    at S, and that placement is fixed (decision 1). Deterministic order:
    RTT-variant rank, then PBS length.
    """
    rtts = enumerate_variants(rtt_aa, max_variants=max_rtt_variants)
    out: List[PbsRttVariant] = []
    for rv in rtts:
        for length in pbs_lengths:
            out.append(PbsRttVariant(spacer=spacer, rtt_dna=rv.dna, pbs_length=length))
    return out
