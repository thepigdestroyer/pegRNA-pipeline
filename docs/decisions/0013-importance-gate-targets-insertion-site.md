# 0013. Structural-importance gate judges the insertion site, not the recode arm

Date: 2026-07-02
Status: Accepted

## Context

After the splice gate was fixed (ADR 0011), the real CTCF run advanced and then failed all
18 candidates at the AlphaFold structural-importance gate ("importance blocked a residue").
The gate computed the residue footprint of the **synonymous recode arm** and blocked the
candidate if any of those residues was buried / packed / special.

This is the wrong target. The recode is **synonymous** -- INV-0 already asserts the protein
over the recode span is byte-for-byte identical in amino-acid sequence -- so the recoded
residues' structure is unchanged and their burial is irrelevant to a synonymous edit.
Because a folded C-terminus is mostly buried/packed, the mis-targeted gate rejected every
candidate for no protein-level reason, and it defeated the core "recode to break
microhomology" strategy in any structured region.

The gate's real purpose (per the design intent) is **insertion-site tolerance**: the cargo
(epitope tag + linker + T7-promoter-leftover junk) is physically inserted at ``S``, and that
insert must not land in a structured / interaction-prone region where it would disrupt the
fold or a functional interface. So the gate must judge the **insertion site**, and the
ability to forbid inserting into a structured region must be retained.

## Decision

Re-target the structural-importance gate from the recode-arm footprint to the
**insertion-site footprint**:

- New `_insertion_site_footprint(site_s, cds_mask, *, context=1)` returns the CDS residue
  indices flanking the insertion point (the residue at ``S-1`` (5') and ``S`` (3', if
  coding), plus ``context`` residues each side).
- `run_recoded_tjpe` gates on `span_blocked(importance_scorer, insertion_footprint)` instead
  of the recode-span footprint. `insertion_importance_context` (default 1) is a parameter.
- The special-residue hard block (catalytic / disulfide / structural-Pro / positive-phi-Gly)
  and the buried/packed criteria (pLDDT >= 70 & RSA < 0.25; contacts > 18) are unchanged --
  they now correctly apply to the insertion locus.
- The diagnostic names the offending insertion-site residue(s) and advises choosing a
  surface / loop / terminal insertion site.

Because ``S`` is fixed (the insert site never moves, Decision 1), this is effectively a
GO/NO-GO on the user's chosen insertion site, constant across candidates -- exactly the
"don't insert a tag into a structured region" check the gate was meant to be.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Keep judging the recode-arm residues | Biologically wrong: synonymous recode leaves the protein unchanged, so arm-residue burial is irrelevant; blocks every structured-region design for no reason |
| Make the gate fully advisory (never block) | User wants to RETAIN the ability to forbid inserting a tag into a structured / interaction-prone region -- that protection is the point of the gate |
| Relax the buried/packed thresholds | Treats a targeting problem as a threshold problem; the real fix is judging the right residues |

## Consequences

- Positive: the gate now protects against the real risk -- a tag/linker/junk insertion
  disrupting the protein -- and no longer spuriously rejects synonymous recodes.
- Positive: a poorly-chosen insertion site ``S`` (buried / in an interface) is still
  correctly rejected with an actionable diagnostic; the special-residue block is retained.
- Positive: unblocks structured-region and C-terminal designs where the insertion locus
  itself is tolerant (surface / terminal), which is the common tagging case.
- Neutral: the human-codon-preference recode (ADR 0009) already mitigates the only residual
  concern about the recode arm (co-translational folding kinetics), so dropping the arm
  check loses no real protection.

## References

- `src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py` (`_insertion_site_footprint`, gate),
  `structure_importance.py` (`residue_blocked`, unchanged criteria).
- `tests/test_recoded_tjpe_workflow.py`; `docs/FEATURES.md`.
- Diagnosis: CX3 job 3164183 ("cand1..18: AlphaFold importance blocked a residue").
