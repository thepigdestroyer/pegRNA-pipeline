# 0006. Genome model: T2T-KOLF2.1J truth, GRCh38 for ML gates

Date: 2026-06-18
Status: Accepted

## Context

The lab edits KOLF2.1J iPSCs. The design must reason about the actual edited genome
(including KOLF2.1J variants and both alleles), but the ML gates (SpliceAI, AlphaFold)
were trained on GRCh38 and only produce calibrated scores in GRCh38 context. A single
genome choice could not satisfy both truth and tool-validity.

## Decision

Use the T2T diploid KOLF2.1J genome as the truth model and GRCh38 coordinates for the
GRCh38-trained ML gates (decision B3-C). In the current implementation this is wired as
**Ensembl GRCh38 ingest plus a KOLF2.1J het-SNP overlay** (`genome_ingest.py`,
`genome_context.py`), because a clean T2T-KOLF FASTA was not resolvable at build time.
The interface is identical, so substituting the T2T FASTA later is a data change, not a
code change. KOLF2.1J is male/XY.

The genome context supplies the local sequence, the CDS mask (drives arm partitioning),
and splice context. Het-SNPs are surfaced as advisory allele-specific options; in
`allele_hunt` mode they may move the ngRNA nick but never the insert site `S`.
Non-coding bases are re-templated byte-exact to the targeted allele so the other allele
is not silently overwritten (Decision 15).

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| GRCh38 only | Ignores KOLF2.1J variants and diploidy; can design against bases the cell does not have. |
| T2T-KOLF only | ML gates are GRCh38-trained; their scores are not calibrated in a different coordinate/sequence frame. |
| Collapse to a haploid consensus | Loses allele-specificity and the ability to target one allele. |

## Consequences

- Positive: design reasons about the real edited line while ML gates stay in their
  trained frame.
- Positive: allele-aware, with explicit byte-exact handling of the untargeted allele.
- Negative: until the true T2T-KOLF FASTA is wired in, the truth model is an
  approximation (Ensembl-GRCh38 + SNP overlay).

## References

- `src/pegrna_pipeline/hpc/genome_ingest.py`, `genome_context.py`, `allele_specificity.py`.
- `docs/decisions/0003-recoded-arm-cds-mask-partition-and-het-edit-pair.md`.
- `docs/DESIGN_HISTORY.md` §6.
