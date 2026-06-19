"""Diploid genomic context for a transcript: exon map + coordinate conversions.

The recoded TJ-PE pipeline (decision 9, B3-C) designs against the **T2T diploid
KOLF2.1J assembly** as the source of truth for all emitted/avoided sequence, on
**both haplotypes**, and maps each locus to GRCh38 only to feed the
GRCh38-trained gates (SpliceAI/Pangolin, AlphaFold-DB).

The real fetch (T2T-KOLF diploid FASTA + per-haplotype annotations + a
KOLF<->GRCh38 chain) is a separate follow-on plan; this module defines the
diploid-aware type, the mRNA<->genomic and KOLF<->GRCh38 conversions, and the
fail-loud guard that forbids any recoding/splice work when context is absent.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple


class MissingGenomicContextError(RuntimeError):
    """Raised when a genome-dependent operation is attempted without context."""


@dataclass(frozen=True)
class ExonBlock:
    g_start: int
    g_end: int
    m_start: int
    m_end: int


@dataclass(frozen=True)
class Haplotype:
    name: str
    genomic_sequence: Optional[str]
    exons: Tuple[ExonBlock, ...]


@dataclass(frozen=True)
class GenomicContext:
    genome_build: str
    chrom: str
    strand: int  # +1 / -1
    haplotypes: Tuple[Haplotype, ...]
    het_snps: Tuple[Dict[str, object], ...] = ()
    #: Optional KOLF->GRCh38 offset map placeholder. A real chain file is wired
    #: in the follow-on plan; until then ``to_grch38`` is the identity.
    grch38_chain: Optional[object] = None

    def to_grch38(self, pos: int) -> int:
        """Map a KOLF-assembly coordinate to GRCh38.

        Placeholder identity map until a real chain file is loaded; the follow-on
        ingest plan replaces this with a liftover. Callers feeding SpliceAI/AFDB
        must use this rather than raw KOLF coordinates.
        """
        if self.grch38_chain is None:
            return pos
        raise NotImplementedError("real KOLF<->GRCh38 chain mapping is a follow-on plan")


def mrna_to_genomic(ctx: Optional[GenomicContext], m_pos: int, *, haplotype: int = 0) -> int:
    """Map an mRNA position to a genomic coordinate on the given haplotype.

    Raises :class:`MissingGenomicContextError` if ``ctx`` is ``None`` (fail-loud:
    recoding/splice work must never silently proceed on mRNA-only data), and
    ``ValueError`` if ``m_pos`` lies outside every exon block.
    """
    if ctx is None:
        raise MissingGenomicContextError(
            "genomic context required (recoding/splice forbidden on mRNA-only data)"
        )
    hap = ctx.haplotypes[haplotype]
    for ex in hap.exons:
        if ex.m_start <= m_pos < ex.m_end:
            return ex.g_start + (m_pos - ex.m_start)  # + strand
    raise ValueError(f"mRNA position {m_pos} not within any exon block")
