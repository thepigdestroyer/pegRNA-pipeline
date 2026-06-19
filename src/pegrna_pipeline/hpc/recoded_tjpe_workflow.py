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

from pegrna_pipeline.hpc.codon_expander import _gc_content, enumerate_variants

__all__ = ["PbsRttVariant", "enumerate_pbs_rtt_variants", "passes_pbs2_gate"]


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


def passes_pbs2_gate(
    pbs2: str, *, repetitive: bool,
    min_len: int = 11, max_len: int = 17,
    gc_lo: float = 0.40, gc_hi: float = 0.70,
    tm_lo: float = 37.0, tm_hi: float = 42.0,
) -> Tuple[bool, str]:
    """Joint PBS2 gate (decision 8): length 11-17, GC 40-70 %, Tm 37-42 C, not
    repetitive.

    Length and Tm are coupled — a short PBS2 rarely reaches 37 C unless GC-rich —
    so the feasible set tends to collapse to ~14-17 nt; that is expected. Tm uses
    Biopython's nearest-neighbour model (imported lazily so the core import stays
    light).
    """
    if not (min_len <= len(pbs2) <= max_len):
        return False, f"length {len(pbs2)} outside [{min_len},{max_len}]"
    gc = _gc_content(pbs2)
    if not (gc_lo <= gc <= gc_hi):
        return False, f"GC {gc:.2f} outside [{gc_lo},{gc_hi}]"
    from Bio.SeqUtils import MeltingTemp as mt  # noqa: PLC0415
    tm = float(mt.Tm_NN(pbs2))
    if not (tm_lo <= tm <= tm_hi):
        return False, f"Tm {tm:.1f} outside [{tm_lo},{tm_hi}]"
    if repetitive:
        return False, "repetitive region"
    return True, "ok"
