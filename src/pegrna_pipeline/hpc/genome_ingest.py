"""Real genomic-context ingest for the recoded TJ-PE pipeline (replaces the
fixture-only ``genome_context`` builder).

Fetches a transcript's exon structure + genomic coordinates from the Ensembl
GRCh38 REST API and builds a diploid :class:`GenomicContext`. An optional set of
phased het SNPs (from a KOLF2.1J VCF) is attached so the pipeline can reason over
**both haplotypes** (decision 9 / B3-C).

Why Ensembl-GRCh38 coordinates (and not the raw KOLF FASTA as truth): the
downstream gates (SpliceAI, AlphaFold-DB) are GRCh38-trained, so all coordinates
must stay in GRCh38 or their scores are silently wrong. A clean diploid KOLF2.1J
assembly now exists (NCBI GCA_058280895.1 ``pri`` + GCA_058280805.1 ``alt``,
NHGRI, 2026-06-17), but making it the truth genome would require a validated
KOLF<->GRCh38 liftover chain and a GENCODE->KOLF annotation projection (follow-on
"scope B"). Instead this module stays in GRCh38 coordinates and materialises the
real diploid **coordinate-safely** (scope A): KOLF2.1J phased substitutions are
applied onto the GRCh38 region sequence per haplotype (:func:`apply_phased_variants`).
Same-length substitutions preserve the coordinate frame, so ``to_grch38`` stays a
valid identity and the GRCh38-trained gates are unaffected. Indels are advisory
only here (indel-accurate per-haplotype truth is scope B).

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
from pegrna_pipeline.hpc.variant_apply import apply_phased_variants

__all__ = ["fetch_genomic_context", "build_genomic_context"]

_ENSEMBL = "https://rest.ensembl.org"


def _real_http_get(url: str, **params) -> dict:  # pragma: no cover - network path
    import requests  # noqa: PLC0415

    headers = {"Content-Type": "application/json"}
    resp = requests.get(url, headers=headers, params=params, timeout=60)
    resp.raise_for_status()
    return resp.json()


def build_genomic_context(
    *,
    lookup: dict,
    mrna_seq: str,
    het_snps: Sequence[dict] = (),
    phased_variants: Sequence[dict] = (),
) -> GenomicContext:
    """Build a diploid :class:`GenomicContext` from an Ensembl lookup payload
    (``expand=1``, so it carries the ``Exon`` list) + the spliced mRNA sequence.

    Both haplotypes share the reference GRCh38 exon map. ``phased_variants`` (KOLF2.1J
    phased calls: ``pos``/``ref``/``alt`` + ``haplotype`` or ``homozygous``) are
    materialised into each haplotype's sequence coordinate-safely via
    :func:`apply_phased_variants` (substitutions only; scope A), so the two haplotypes
    can genuinely differ while staying in GRCh38 coordinates. When ``phased_variants``
    is empty, both haplotypes equal the reference sequence (legacy behaviour).

    ``het_snps`` remains the advisory allele-specificity overlay (element/offset/phased
    metadata consumed by :mod:`allele_specificity`); it is orthogonal to the sequence
    materialisation done from ``phased_variants``.
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
    # ``mrna_seq`` is the forward-strand sequence of the region min..max exon; its first
    # base sits at the minimum genomic exon start (1-based). That anchors variant coords.
    region_start = min(int(ex["start"]) for ex in exons_raw)
    hap1_seq = apply_phased_variants(
        mrna_seq, phased_variants, region_start=region_start, haplotype=0
    )
    hap2_seq = apply_phased_variants(
        mrna_seq, phased_variants, region_start=region_start, haplotype=1
    )
    hap1 = Haplotype(name="hap1", genomic_sequence=hap1_seq, exons=exon_tuple)
    hap2 = Haplotype(name="hap2", genomic_sequence=hap2_seq, exons=exon_tuple)
    build = "GRCh38 (Ensembl) coordinates"
    build += (
        " + materialised KOLF2.1J diploid (phased substitutions)"
        if phased_variants
        else " + KOLF2.1J variant overlay"
    )
    return GenomicContext(
        genome_build=build,
        chrom=str(chrom),
        strand=strand,
        haplotypes=(hap1, hap2),
        het_snps=tuple(het_snps),
    )


def fetch_genomic_context(
    transcript_id: str,
    *,
    het_snps: Optional[Sequence[dict]] = None,
    phased_variants: Optional[Sequence[dict]] = None,
    http_get: Callable[..., dict] = _real_http_get,
) -> GenomicContext:
    """Fetch exon structure + sequence for an Ensembl transcript and build a
    diploid :class:`GenomicContext`. ``http_get`` is injectable for testing.

    ``phased_variants`` (KOLF2.1J phased calls) are materialised into each haplotype
    coordinate-safely (see :func:`build_genomic_context`). Raises
    :class:`MissingGenomicContextError` if Ensembl returns no usable transcript
    (fail-loud — never silently proceed on absent genomic data).
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
    return build_genomic_context(
        lookup=lookup,
        mrna_seq=mrna_seq,
        het_snps=het_snps or (),
        phased_variants=phased_variants or (),
    )
