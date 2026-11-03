# 0026. Report when the extension pairs with the scaffold; do not reject on it

Date: 2026-11-03
Status: Accepted. Settles issue #27.

## Context

ADR 0024's sweep found a gap: nothing checks whether the pegRNA's 3' extension pairs with its own
scaffold.

* The folding gate counts cross-region base pairs but **excludes every pair that touches the
  scaffold**, by design (the scaffold folds on itself).
* The docking gate catches extensive sequestration (every complement of 40-60 nt over the
  scaffold's 5' two-thirds) but only 4 of 12 24-nt complements, and it is blind to the 3' hairpin.

A realistic way to hit this: an insert that carries a guide-RNA scaffold, so that its reverse
complement sits in the extension and pairs with the pegRNA's own scaffold.

## Evidence

The measure is the number of **scaffold nucleotides paired with the extension at p >= 0.99** in
ViennaRNA's partition function, the folding gate's own calibrated cutoff.

| designs | n | scaffold nt paired |
|---|---|---|
| synthetic neutral extensions, 126-336 nt | 56 | 0-11, median 0 |
| built to pair with the scaffold (the sweep) | 12 | **18-29**, median 21 |
| the pipeline's own two designed pegRNAs | 2 | 0-1 |
| **real candidate pool** (CTCF: 8 cargos, 2 spacers, 4 nicks; 194-317 nt) | 130 | 0-9, median 0 |

* The measure **does not grow with length**: the longest neutral designs (256-336 nt) pair at most
  5. Under ViennaRNA's MFE structure, by contrast, a random 296 nt extension pairs with up to 24
  scaffold nucleotides by chance, which would have made a length-dependent threshold.
* The probability cutoff matters. The same 130 real pegRNAs at p >= 0.5 have a median of 11 and 16
  of them would exceed 15; at 0.99 none do. 0.5 is noise, 0.99 is a signal.
* The two classes are separated with margin: 11 against 18.

## Decision

An **advisory** (`THRESHOLDS.scaffold_pairing_advisory_nt = 15`, between 11 and 18) computed with
`scaffold_extension_pairing` from the probability map the folding gate already obtains (the
sequence is folded once and shared, so a GPU folder does not run twice). It uses the folding gate's
own resolved cutoff. When it reaches 15 the design carries a warning, shown in the "Check these
before you order" callout; the structure row states the count for every design. **Nothing is
rejected on it.**

Why not reject. Whether pairing at this level lowers editing for a real pegRNA is not established
here. The sweep constructs are perfect complements built to pair; real designs sit somewhere
between, and none of the 130 real ones does. It is widely reported that a 3' extension folding back
on the guide can lower prime-editing efficiency, but nothing in this repository measures it for the
scaffold, and the numbers above are a model's structure, not an experiment. A rejecting rule on an
unvalidated signal repeats the mistake ADR 0024 corrected. Making it reject should wait for
literature or a wet-lab control.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Reject at the threshold | Unvalidated biology; no real design is near it, so there is also nothing to gain yet |
| MFE count instead of a probability | Grows with extension length (24 nt for random 296 nt), so no single threshold works |
| p >= 0.5 | Median 11 on real designs: it would fire constantly |
| Count scaffold-scaffold pairs too | They are the scaffold's own fold and always present |

## Consequences

* It stays quiet on ordinary designs (0 of 130 flagged) and fires on constructs that tie up the
  scaffold (12 of 12). A synthetic cargo containing the scaffold is flagged (tested).
* The threshold was calibrated on ViennaRNA. For another folder (RiNALMo) the same count is
  unvalidated, which is acceptable for a warning and not for a gate.
* One locus, a synthetic sweep, and no wet-lab data. The three-way separation is a property of a
  model's structure prediction.

## References

* `src/pegrna_pipeline/hpc/full_structure.py` (`scaffold_extension_pairing`),
  `gate_thresholds.py`, `recoded_tjpe_workflow.py`.
* Tests: `tests/test_scaffold_pairing_advisory.py` (pins the calibration to the committed data and
  to real pool candidates). ADR 0024 (the sweep).
