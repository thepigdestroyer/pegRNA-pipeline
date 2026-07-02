# 0012. Colour-annotated verification DOCX for manual cross-checking

Date: 2026-07-02
Status: Accepted

## Context

A winning recoded TJ-PE design is emitted as a text report and a component CSV, but
neither makes it easy to verify *by eye* that the pipeline assembled a functional pegRNA
with nothing missed or misplaced. The user asked for a dedicated verification file that
shows, colour-coded and annotated: (1) the genomic context annotated with every region the
pegRNA/tjPE-pegRNA binds or uses and every region that will change; (2) the full pegRNA
broken into every component; (3) the predicted fully-edited genomic DNA annotated with the
insert and every changed base; (4) all component sequences used to assemble the insert and
the winning pegRNA. Output as a DOCX with text colour, highlight, annotations, and legends.

## Decision

Add `verification.py` with a **pure, renderer-agnostic annotation model** separated from a
**lazy python-docx renderer**:

- Model: `Annotation(start, end, label, category)`, `AnnotatedSequence`, `ComponentRow`,
  `VerificationReport`, and `build_verification_report(result, *, genome, cargo)`. The model
  is where verification correctness lives (span coordinates, categories) and is unit-tested
  exhaustively **without** python-docx.
- Renderer: `render_verification_docx(report, out_path)` — lazily imports python-docx and
  colour/highlights each base by the highest-priority annotation category, with a legend,
  per-section annotation tables, and a component table. Raises `ImportError` if python-docx
  is absent (so the model stays dependency-free).
- Section 3 (edited allele) computes component spans from the **known construction geometry**
  (`edited = genome[:S] + cargo + arm_sense + genome[nick:]`, deletion_size == 0), not a
  common-prefix/suffix diff — a diff would absorb the matching tail of the recoded arm and
  mis-size it, so the changed-base highlighting would silently drop.
- Surfaced `edited_allele`, `recode_span`, `locked_span` on `RecodedTJPEResult` (loop-computed
  facts needed for sections 1 and 3).
- `design_recoded` emits `recoded_tjpe_verification.docx` beside the text/CSV, **guarded**: if
  python-docx is not installed the run still succeeds (text + CSV written, DOCX skipped).
  python-docx is declared in a new `[report]` extra and staged in `environment.yml`.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Annotate the existing text report only | No colour/highlight; hard to scan; the user asked for a visual DOCX |
| Emit HTML instead of DOCX | User specified DOCX (Word colour/highlight/legend) |
| Put python-docx calls inside the model | Couples the model to a heavy optional dep; not CI-testable without it. The model/renderer split keeps the correctness-critical part dependency-free |
| Diff-infer section 3 spans (prefix/suffix) | Absorbs the arm's matching tail -> wrong arm size -> changed-base highlights silently dropped; geometry is exact |

## Consequences

- Positive: a single DOCX lets a user verify the whole design by eye -- binding sites,
  changed bases, every component -- with a legend and per-section annotation tables.
- Positive: the annotation model is fully CI-tested (spans, categories, geometry) without
  python-docx; the renderer is a thin layer tested where python-docx is installed.
- Positive: fully backward compatible -- text + CSV unchanged; the DOCX is additive and
  skipped gracefully where python-docx is absent.
- Neutral: the DOCX is emitted only on a winning design (as with the text/CSV). python-docx
  is a normal PyPI dependency (unlike maxentpy), staged via the `[report]` extra.

## References

- `src/pegrna_pipeline/hpc/verification.py`, `recoded_driver.py` (emission),
  `recoded_tjpe_workflow.py` (`edited_allele`/`recode_span`/`locked_span`).
- `tests/test_verification.py`; `docs/FEATURES.md`.
