# 0007. Coordinate-safe diploid materialization (KOLF2.1J on GRCh38)

Date: 2026-07-01
Status: Accepted (partially supersedes the "data-change-only" framing in 0006)

## Context

ADR 0006 chose a T2T-KOLF2.1J truth model with GRCh38 coordinates for the ML gates, and
recorded that the raw diploid FASTA was not a resolvable pinned asset, so the pipeline
used an Ensembl-GRCh38 + het-SNP overlay whose two haplotypes shared the reference
sequence (a structural, not sequence-materialised, diploid).

Two facts changed the decision:

1. A clean diploid KOLF2.1J assembly now exists on NCBI: GCA_058280895.1 (`pri`) and
   GCA_058280805.1 (`alt`), submitted by NHGRI on 2026-06-17, chromosome-level.
2. Making that FASTA the *truth* genome is NOT a drop-in data change. SpliceAI and
   AlphaFold-DB are GRCh38-trained, and `genome_context.to_grch38` is an identity
   placeholder. A truth swap needs a validated KOLF<->GRCh38 liftover chain plus a
   GENCODE->KOLF annotation projection; doing it naively feeds KOLF coordinates to
   GRCh38-trained gates and silently corrupts their scores.

## Decision

Implement **scope A** now: stay in GRCh38 coordinates and materialise the real diploid by
applying KOLF2.1J phased **substitutions** onto the GRCh38 region sequence, per
haplotype (`variant_apply.apply_phased_variants`, wired into
`genome_ingest.build_genomic_context` / `fetch_genomic_context` via a new optional
`phased_variants` argument).

- Same-length substitutions preserve the coordinate frame and sequence length for both
  haplotypes, so `to_grch38` stays a valid identity and the GRCh38-trained gates are
  unaffected. Zero ML-gate risk.
- Ref-base mismatch is fail-loud (`ValueError`): it means the wrong build/coordinates.
- Indels are advisory only (skipped for the coordinate-anchored sequence); a single
  fixed GRCh38 frame cannot represent an indel on one haplotype without downstream drift.
- Empty `phased_variants` reproduces the legacy behaviour (both haplotypes == reference),
  so the change is backward compatible.

Defer **scope B** (ingest the KOLF assembly FASTAs as truth, build + validate a
KOLF<->GRCh38 chain, GENCODE->KOLF annotation projection, indel-accurate per-haplotype
coordinates) as a separate, CX3-run follow-on now that the assembly is available.

## Alternatives considered

| Alternative | Why rejected (for now) |
|---|---|
| Swap the KOLF FASTA in as truth with the identity `to_grch38` | Silently feeds KOLF coords to GRCh38-trained SpliceAI/AlphaFold; corrupts scores. |
| Full scope B immediately | Multi-day liftover + annotation-projection + validation; higher risk before the chain is validated; not required to get a real diploid. |
| Splice indels into the GRCh38-anchored sequence | Downstream coordinate drift; a fixed frame cannot hold an indel on one haplotype. |
| Keep hap1 == hap2 (status quo) | Not a real diploid; wastes the now-available alt-haplotype information. |

## Consequences

- Positive: a genuine diploid (hap1 != hap2) with zero risk to the GRCh38-trained gates.
- Positive: fail-loud ref check catches build/coordinate errors immediately.
- Positive: backward compatible; existing callers unaffected.
- Negative: indels are not represented in the per-haplotype sequence (advisory only)
  until scope B; substitution-only diploid.
- Neutral: the KOLF phased-variant *data source* (a GRCh38 phased VCF, e.g. from dipcall
  of the assemblies) is a plug-in input; this ADR delivers the application capability.

## References

- `src/pegrna_pipeline/hpc/variant_apply.py`, `genome_ingest.py`; `tests/test_variant_apply.py`.
- `docs/decisions/0006-genome-model-t2t-kolf-grch38.md`.
- `docs/superpowers/plans/2026-07-01-coordinate-safe-diploid.md`.
- NCBI GCA_058280895.1 (`KOLF2.1Jv1.1.pri`), GCA_058280805.1 (`KOLF2.1Jv1.1.alt`), NHGRI 2026-06-17.
- `docs/DESIGN_HISTORY.md` §6.
