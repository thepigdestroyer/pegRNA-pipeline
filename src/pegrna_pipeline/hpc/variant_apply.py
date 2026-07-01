"""Coordinate-safe diploid materialization (scope A).

Applies KOLF2.1J phased **substitutions** onto a GRCh38 region sequence to build a
per-haplotype sequence, staying entirely in GRCh38 coordinates. Because only
same-length substitutions are applied, the coordinate frame and sequence length are
preserved for both haplotypes, so :meth:`GenomicContext.to_grch38` stays a valid
identity and the GRCh38-trained gates (SpliceAI, AlphaFold-DB) are unaffected.

Indels are intentionally NOT spliced into the coordinate-anchored sequence: a single
fixed GRCh38 coordinate frame cannot represent an insertion/deletion on one haplotype
without drifting downstream coordinates. Indel-accurate per-haplotype truth is scope B
(ingest the KOLF assembly FASTAs + a validated KOLF<->GRCh38 chain). Indels are skipped
here (advisory only).

A variant is a dict with at least ``pos`` (1-based genomic), ``ref``, ``alt``. It
applies to a haplotype when ``homozygous`` is truthy, or its ``haplotype`` field equals
the target haplotype index. A ref-base mismatch is fail-loud (it means the wrong build
or coordinates), never silently ignored.
"""

from __future__ import annotations

from typing import Sequence

__all__ = ["apply_phased_variants"]


def _applies_to(variant: dict, haplotype: int) -> bool:
    if variant.get("homozygous"):
        return True
    hap = variant.get("haplotype")
    return hap is not None and int(hap) == haplotype


def apply_phased_variants(
    reference: str,
    variants: Sequence[dict],
    *,
    region_start: int,
    haplotype: int,
) -> str:
    """Return ``reference`` with the phased substitutions for ``haplotype`` applied.

    ``region_start`` is the 1-based genomic coordinate of ``reference[0]``. Only
    substitutions (``len(ref) == len(alt)``) inside the window are applied; indels and
    out-of-window variants are skipped. Raises :class:`ValueError` if a variant's
    declared ``ref`` does not match ``reference`` at its position (wrong build/coords).
    """
    if not variants:
        return reference
    seq = list(reference)
    n = len(seq)
    for v in variants:
        if not _applies_to(v, haplotype):
            continue
        ref = str(v["ref"]).upper()
        alt = str(v["alt"]).upper()
        if len(ref) != len(alt):
            # Indel: out of scope for coordinate-anchored materialization (scope B).
            continue
        idx = int(v["pos"]) - region_start
        if idx < 0 or idx + len(ref) > n:
            # Variant lies outside the fetched window.
            continue
        observed = "".join(seq[idx : idx + len(ref)]).upper()
        if observed != ref:
            raise ValueError(
                f"ref-base mismatch at pos {v['pos']}: variant claims {ref!r} but "
                f"reference has {observed!r} (wrong genome build or coordinates)"
            )
        seq[idx : idx + len(alt)] = list(alt)
    return "".join(seq)
