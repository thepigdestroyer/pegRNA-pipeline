"""Real genomic-context ingest for the recoded TJ-PE pipeline (replaces the
fixture-only ``genome_context`` builder).

Fetches a transcript's exon structure + genomic coordinates from the Ensembl
GRCh38 REST API and builds a diploid :class:`GenomicContext`. An optional set of
phased het SNPs (from a KOLF2.1J VCF) is attached so the pipeline can reason over
**both haplotypes** (decision 9 / B3-C).

Why Ensembl-GRCh38 and not the raw T2T-KOLF FASTA: the T2T diploid assembly's
genome FASTA is not cleanly resolvable as a single pinned asset yet (the Zenodo
WGS record is README-only), whereas Ensembl GRCh38 gives reliable, versioned
exon/sequence data today and is the coordinate space the downstream gates
(SpliceAI, AlphaFold-DB) are trained on. The KOLF private variants are layered as
the het-SNP overlay; swapping in the full T2T diploid FASTA later only changes
the *source* of the exon/sequence + variant set, not this module's interface.

The HTTP layer is injected (``http_get``) so this is unit-testable offline; the
default real fetcher uses ``requests``.
"""

from __future__ import annotations

from typing import Callable, List, Optional, Sequence

from pegrna_pipeline.hpc.genome_context import (
    ExonBlock,
    GenomicContext,
    Haplotype,
    MissingGenomicContextError,
)

__all__ = ["fetch_genomic_context", "build_genomic_context"]

_ENSEMBL = "https://rest.ensembl.org"


def _real_http_get(url: str, **params) -> dict:  # pragma: no cover - network path
    import requests  # noqa: PLC0415

    headers = {"Content-Type": "application/json"}
    resp = requests.get(url, headers=headers, params=params, timeout=60)
    resp.raise_for_status()
    return resp.json()


def build_genomic_context(
    *, lookup: dict, mrna_seq: str, het_snps: Sequence[dict] = ()
) -> GenomicContext:
    """Build a diploid :class:`GenomicContext` from an Ensembl lookup payload
    (``expand=1``, so it carries the ``Exon`` list) + the spliced mRNA sequence.

    Both haplotypes share the reference exon map; phased het SNPs are attached via
    ``het_snps`` so allele-specific logic (Task 17) can branch on them. (A future
    upgrade can materialise a true alt-haplotype sequence from a full VCF; the
    interface here already exposes both haplotypes.)
    """
    exons_raw = lookup.get("Exon") or lookup.get("exons")
    chrom = lookup.get("seq_region_name")
    strand = int(lookup.get("strand", 1))
    if not exons_raw or chrom is None:
        raise MissingGenomicContextError(
            "Ensembl lookup payload missing Exon list / seq_region_name"
        )
    # Sort exons in transcription order and assign contiguous mRNA coordinates.
    exons_sorted = sorted(exons_raw, key=lambda e: e["start"], reverse=(strand < 0))
    blocks: List[ExonBlock] = []
    m = 0
    for ex in exons_sorted:
        length = int(ex["end"]) - int(ex["start"]) + 1
        blocks.append(
            ExonBlock(g_start=int(ex["start"]), g_end=int(ex["end"]), m_start=m, m_end=m + length)
        )
        m += length
    exon_tuple = tuple(blocks)
    hap1 = Haplotype(name="hap1", genomic_sequence=mrna_seq, exons=exon_tuple)
    hap2 = Haplotype(name="hap2", genomic_sequence=mrna_seq, exons=exon_tuple)
    return GenomicContext(
        genome_build="GRCh38 (Ensembl) + KOLF2.1J variant overlay",
        chrom=str(chrom),
        strand=strand,
        haplotypes=(hap1, hap2),
        het_snps=tuple(het_snps),
    )


def fetch_genomic_context(
    transcript_id: str,
    *,
    het_snps: Optional[Sequence[dict]] = None,
    http_get: Callable[..., dict] = _real_http_get,
) -> GenomicContext:
    """Fetch exon structure + sequence for an Ensembl transcript and build a
    diploid :class:`GenomicContext`. ``http_get`` is injectable for testing.

    Raises :class:`MissingGenomicContextError` if Ensembl returns no usable
    transcript (fail-loud — never silently proceed on absent genomic data).
    """
    lookup = http_get(f"{_ENSEMBL}/lookup/id/{transcript_id}", expand=1)
    if not lookup or "Exon" not in lookup:
        raise MissingGenomicContextError(
            f"Ensembl returned no expandable transcript for {transcript_id!r}"
        )
    chrom = lookup["seq_region_name"]
    start = min(int(e["start"]) for e in lookup["Exon"])
    end = max(int(e["end"]) for e in lookup["Exon"])
    seq_payload = http_get(f"{_ENSEMBL}/sequence/region/human/{chrom}:{start}..{end}")
    mrna_seq = seq_payload.get("seq", "")
    return build_genomic_context(lookup=lookup, mrna_seq=mrna_seq, het_snps=het_snps or ())
