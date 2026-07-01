# Plan: coordinate-safe diploid materialization (scope A)

Date: 2026-07-01
Status: In progress

## Context

A clean diploid KOLF2.1J assembly now exists on NCBI (GCA_058280895.1 `pri`,
GCA_058280805.1 `alt`, NHGRI, 2026-06-17), so the availability blocker recorded in
`genome_ingest.py` is stale. But making that FASTA the *truth* genome would break the
GRCh38-trained gates (SpliceAI, AlphaFold) unless a validated KOLF<->GRCh38 chain and a
GENCODE->KOLF annotation projection exist. `genome_context.to_grch38` is still an
identity placeholder. Doing the full swap naively silently corrupts the ML gates.

Chosen scope: **A — coordinate-safe diploid.** Materialize the real second haplotype by
applying KOLF2.1J phased variants onto the GRCh38 region sequence, staying entirely in
GRCh38 coordinates. Zero ML-gate risk (`to_grch38` identity stays correct). Full
KOLF-as-truth (chain + annotation projection) remains scoped follow-on B.

## Design

New pure module `variant_apply.py`:

- `apply_phased_variants(reference, variants, *, region_start, haplotype) -> str`
  - `region_start` = 1-based genomic coordinate of `reference[0]`.
  - Applies only **substitutions** (`len(ref) == len(alt)`) so the coordinate frame and
    sequence length are preserved for BOTH haplotypes (the point of scope A).
  - A variant applies to a haplotype when `homozygous` is truthy, or its `haplotype`
    field equals the target haplotype index.
  - **Fail-loud**: if `reference[idx:idx+len(ref)] != ref`, raise `ValueError` (wrong
    build/coordinates). Variants outside the window are ignored. Indels are skipped for
    the coordinate-anchored sequence (advisory-only; indel-as-truth needs scope B).

Wire into `genome_ingest.build_genomic_context` / `fetch_genomic_context`:

- New optional `phased_variants: Sequence[dict] = ()`. When empty, behaviour is
  unchanged (hap1 == hap2 == reference: backward compatible).
- `region_start` derived from the min exon start (same value already computed in
  `fetch_genomic_context`).
- hap1 = apply(..., haplotype=0); hap2 = apply(..., haplotype=1).
- `genome_build` string updated to state GRCh38 coordinates with a materialized
  KOLF2.1J diploid overlay.

## Tests (RED first) — `tests/test_variant_apply.py`

1. het SNV on hap2 only: hap0 == ref, hap1 changed at the site.
2. homozygous SNV: both haplotypes changed.
3. ref-base mismatch raises `ValueError` (build/coordinate safety).
4. variant outside window ignored.
5. insertion/deletion skipped (sequence unchanged) in coordinate-safe mode.
6. empty variants returns the reference unchanged (backward compat).
7. multiple SNVs on one haplotype all applied (order-independent).
8. `build_genomic_context` materializes a real diploid (hap1 != hap2) while keeping
   `to_grch38` identity and a GRCh38 `genome_build`.

## Docs

- Update `genome_ingest.py` + `genome_context.py` headers (assembly now exists; diploid
  now materialized coordinate-safe; full-truth = follow-on B).
- ADR 0007 (coordinate-safe diploid materialization); note on ADR 0006.
- `docs/DESIGN_HISTORY.md` §6.

## Out of scope (follow-on B)

Ingesting the assembly FASTAs as truth, building/validating a KOLF<->GRCh38 liftover
chain, GENCODE->KOLF annotation projection, and indel-accurate per-haplotype coordinates.
