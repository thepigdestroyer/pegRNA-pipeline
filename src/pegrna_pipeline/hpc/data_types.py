"""Shared frozen dataclasses used across Plan-2 workflows.

This module intentionally has no runtime dependencies beyond the standard
library so it can be imported by any layer of the pipeline (tests, HPC
jobs, the Click CLI) without dragging in Biopython or pandas.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GeneRecord:
    """Immutable snapshot of a gene's mRNA context.

    Attributes
    ----------
    gene_name:
        The HGNC-style symbol the caller asked for (case-preserved).
    refseq_accession:
        The RefSeq accession that was resolved (e.g. ``NM_006565.4``).
    mrna_sequence:
        Full mRNA sequence as an uppercase ACGT string (no Us).
    cds_start, cds_end:
        Half-open indices of the CDS within ``mrna_sequence``
        (``[cds_start, cds_end)``), 0-indexed.
    stop_codon_start:
        Index of the first nucleotide of the stop codon, i.e.
        ``cds_end - 3``.
    stop_codon:
        3-letter stop codon (``TAA`` | ``TAG`` | ``TGA``).
    upstream_flank:
        ``flanking_bp`` nucleotides immediately preceding the stop codon
        (ends at ``stop_codon_start - 1``).
    downstream_flank:
        ``flanking_bp`` nucleotides immediately after the stop codon
        (starts at ``stop_codon_start + 3``).
    flanking_bp:
        Length of the flanks actually returned.
    """

    gene_name: str
    refseq_accession: str
    mrna_sequence: str
    cds_start: int
    cds_end: int
    stop_codon_start: int
    stop_codon: str
    upstream_flank: str
    downstream_flank: str
    flanking_bp: int
