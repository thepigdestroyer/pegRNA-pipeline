# 0001. TJ-PE ngRNA safety-gated selection

Date: 2026-06-08
Status: Accepted

## Context

The TJ-PE pipeline produced a CTCF design whose nicking sgRNA was 19 nt and had four exact genomic hits. Zheng et al. 2023 define TJ-PE as one TJ-pegRNA plus one nicking sgRNA; the nicking sgRNA must expose the PBS2-compatible second nick, so selecting a multi-hit ngRNA is a biological safety and validity failure. The implementation must remain auditable on CX3 and must not silently choose an unsafe ngRNA.

## Decision

TJ-PE will rank pegRNA/variant candidates by PRIDICT2 PE2 score, then try them in PE2 order and accept the first candidate whose compatible ngRNA survives strict TJ-PE filters: opposite strand, 19-21 nt with 20 nt preferred in ranking, default absolute nick distance 70-90 bp, sequence sanity filters, and exact genomic uniqueness (`mm0_count == 1`, where the one exact hit is presumed to be the intended locus). ngRNA candidates are ranked only after uniqueness filtering. Reports will expose the ngRNA selection window, spacer length, exact-hit count, and all TJ-pegRNA components including PBS1 and RC-PBS2.

## Alternatives considered

| Alternative | Pros | Cons | Why rejected |
|---|---|---|---|
| Keep current selection (`abs(distance)` 40-100, DeepCas9 only, post-hoc off-target scoring) | Minimal code change | Can select multi-hit ngRNAs; ignores paper-backed TJ-PE large-insert window; off-target safety comes too late | Already produced an unsafe CTCF ngRNA with `mm0_count = 4` |
| Pick the top PE2 pegRNA and fail if it lacks a safe ngRNA | Preserves a literal two-stage interpretation | Brittle: a near-top pegRNA may have a safe ngRNA while rank 1 does not | Rejected after consensus review; top-N PE2 fallback is safer and more useful |
| Jointly score every pegRNA/ngRNA pair globally | Theoretically optimal | Requires scanning thousands of ngRNAs genome-wide per run; high CX3 cost and harder debugging | Too expensive for current R/Biostrings scorer and not required by the user's stated workflow |

## Consequences

- Positive: TJ-PE cannot silently emit a multi-hit ngRNA as the final design.
- Positive: If the best PE2 pegRNA lacks a safe nicking partner, the pipeline can still find the best safe design by PE2-ranked fallback.
- Negative: Runs can fail explicitly when no safe ngRNA exists in the default 70-90 bp TJ-PE window.
- Neutral: DeepSpCas9/DeepHF/CFD remain reported when available, but the current scorer still needs follow-up work to compute non-NA values from genomic context.

## References

- Zheng et al. 2023, Nature Communications 14:3369, Fig. 1a and Results: TJ-pegRNA + nicking sgRNA, 90 bp examples at AAVS1/HEK3.
- Supplementary Table 1, Zheng et al. 2023: sgRNA scaffold and nicking sgRNA sequences.
- CTCF run `CTCF_tjpe_ecf7a6b8`: selected ngRNA `GCAGAGGGAGAGGAAGAGG`, `mm0_count = 4`.
