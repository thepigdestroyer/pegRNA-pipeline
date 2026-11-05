# 0028. The SpliceAI cross-check binds only on real genomic context; without it the number is advisory

Date: 2026-11-04
Status: Accepted. Amends ADR 0011.

## Context

A real run of the CTCF example with the 501-bp GGGGS-FLAG-2A-resistance cassette (`api` mode, every gate
a real model) was rejected on all five candidates by the SpliceAI cross-check:

    cand1..cand5: splice (SpliceAI cross-check): splice delta 0.336 >= 0.2 (predicted disruption)

The same example with the 24-bp FLAG passes the same gate with 0.150. The candidates all passed the
primary splice gate, MaxEntScan (best created site 5.0 bits against the strong-site cutoff of 8.0, best
native site 3.6).

ADR 0011 settled what SpliceAI is for here: an **advisory** cross-check that binds only "if a caller supplies
>= 10 kb context", because a 500-bp non-genomic insertion is out of its distribution and MaxEntScan, which needs
no context, is the primary gate. On a short window it raises, and the walk records an "advisory skip".
`docs/FEATURES.md` still says "skipped on short context".

The `api` bundle (added later) wraps the scorer in `_PaddedSpliceAIScorer`, which pads the 297-nt window with N
up to 10 kb. A number then comes back where the original code raised, and the walk applied it as a hard gate.
The documented behaviour and the actual one had drifted apart.

## Evidence

`scripts/analysis/probe_splice_gate.py`, run on the DSI box with the real SpliceAI weights, on the failing
config (ENST00000264010, S = 147, 297-nt window, 501-nt cassette). It first reproduces the
failure through the real pipeline path (5 calls, all 0.336), then takes the number apart.

**Where 0.336 comes from.** `SpliceAIScorer.delta` is the largest of several terms:

| term | value | where |
|---|---|---|
| a native donor site in the right flank, its score changed by the edit | **0.336** | 14 nt 3' of the insert, not in it |
| strongest donor anywhere in the inserted span ("gain") | 0.321 | cargo +272 of 501 |
| strongest acceptor in the inserted span | 0.133 | cargo +75 |
| strongest native site lost from the replaced span | 0.008 | ref, 13 nt 5' of S |

The 24-bp FLAG example, run through the same decomposition in the same arm, gives 0.150, exactly what the shipped
gate reported, so the probe measures what the gate measures.

**The "gain" term is a maximum over every inserted base, so it grows with the insert.** Prefixes of the cassette:

| insert length (nt) | 24 | 60 | 120 | 240 | 360 | 501 |
|---|---|---|---|---|---|---|
| gain | 0.025 | 0.024 | 0.066 | 0.175 | 0.160 | 0.321 |

**Padded with N, an insert of this size cannot pass.** The same alleles with the cassette replaced by something
else of the same length:

| 501-nt insert | n | min | median | max | >= 0.2 |
|---|---|---|---|---|---|
| random synonymous recodings of the same protein | 10 | 0.286 | 0.511 | 0.712 | **10 / 10** |
| random sequences | 10 | 0.288 | 0.639 | 0.774 | **10 / 10** |
| this cassette | 1 | | 0.336 | | 1 / 1 |

A gate that rejects 20 of 20 unrelated sequences, and this cassette sits at the bottom of that range, is not
measuring the cassette.

**With the real flanks it is informative.** Ensembl is reachable from the DSI box. The design window sits at exactly
one position of ENST00000264010's region (chr16:67562526-67639177, + strand), and the real 5 kb either side
replaces the N:

| | N-padded | real +/-5 kb |
|---|---|---|
| this cassette, overall | 0.336 | **0.159** |
| the native donor 14 nt downstream | 0.336 | 0.159 |
| strongest site inside the cassette | 0.321 | **0.007** |
| random recodings (n = 6), min / median / max | 0.29 / 0.51 / 0.71 | 0.165 / 0.448 / 0.838 |
| random recodings that reach 0.2 | 10 / 10 | 5 / 6 |

With real context the real cassette passes and most random recodings of it do not, so the check separates a
designed sequence from an arbitrary one. That is the case ADR 0011 said should bind.

**Not measured.** One locus, one cassette. Ten and six random sequences are a small null, and random codon usage
is harsher than a real expression cassette, so "5 of 6 fail" says the check discriminates, not what fraction of
real cassettes it would pass. Nothing here validates 0.2 (the SNV high-recall cutoff) for a 500-bp insert, and ADR
0011's point that such an insert is out of distribution still stands; RT-PCR across the junctions is the test.

## Decision

1. **A splice scorer says whether its verdict may reject.** `binding` (default `True`, so every existing scorer and
   test is unchanged) and `context_note`. The `api` scorer is binding only when it holds complete real flanks (or
   was handed >= 10 kb by its caller). Otherwise the walk **reports and does not reject**: the gate result reads
   "advisory only" with the context it was computed in; a number at or above the cutoff becomes a warning on the
   design that says why and what the primary gate found; the progress event is a pass with no cutoff to contradict
   it.
2. **The driver gives the `api` scorer the real flanks when a transcript is configured.**
   `fetch_flanking_sequence` asks Ensembl for the transcript's own strand (`fetch_genomic_context` returns the
   forward strand only, which is wrong for a minus-strand gene), locates the window by exact match and requires
   exactly one, takes 5 kb either side, and N-fills where the genome ends (then `complete` is false and the verdict
   is advisory). The wait is bounded (2 attempts of 45 s unless the operator set their own) and cached. **Any**
   failure (no transcript, Ensembl unreachable, the window not found, an applied variant) leaves the scorer on
   N-padding, and the run manifest says so (`splice_context`).
3. **Unchanged:** MaxEntScan stays the primary hard splice gate; `real` mode (the HPC scorer) is untouched; the 0.2
   cutoff is untouched.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Raise the cutoff for the padded case | Twenty unrelated inserts score 0.29-0.77, so a threshold that passes them is above 0.8 and means nothing |
| Normalise by insert length | An invented statistic with no validation behind it |
| Drop the padded scorer (back to ADR 0011's "dormant") | Throws away a number that real context makes useful, and the `api` bundle was written to show it |
| Always bind, and fetch context so it rarely matters | Turns an unreachable Ensembl into a failed design for a check that is advisory by decision |
| Fetch context in every mode | HPC compute nodes usually have no egress and `real` mode's scorer is the raw one; the change stays where the problem was |
| Reuse `fetch_genomic_context` for the sequence | Forward strand only |

## Consequences

* The CTCF cassette run is no longer stopped by a number that cannot discriminate. With real context it is 0.159
  and passes as a binding verdict; if Ensembl is down it passes with a visible warning.
* A SpliceAI rejection now means something: it was computed in the real sequence around the edit.
* Both outcomes are reported, never silent: the gate result and the manifest name the context, and a high
  advisory number is a warning in the "check before you order" list.
* The first design for a transcript that is not in the Ensembl cache pays one fetch (tens of seconds for a large
  gene); later runs read the cache. A window carrying a variant relative to GRCh38, or in the wrong orientation, is
  not located and falls back to advisory.
* Tests: `tests/test_splice_context.py` (31) and `tests/test_structure_diagnostics.py` (15, the addendum): the walk's binding and advisory paths and the progress event, the
  scorer's padding and flank handling, the fetch (strand, margin, N-fill, ambiguity, refusal), the driver's
  fail-soft wiring end to end. Mutating the walk back to "always binds" fails six of them.

## Addendum: the next gate on the same cassette (folding)

With the splice gate fixed, the cassette run went on to fail every candidate at the folding gate. It looked like the
same story (a gate calibrated on 139-228 nt pegRNAs, now meeting 671 nt), so it was measured before anything was
changed (`scripts/analysis/probe_structure_gate.py`, real models, then re-checked locally with the same ViennaRNA).
**It is not the same story, and the gate is unchanged.**

* **One helix, the same in every candidate.** The first four bases of the rcPBS2 (`GGGG`) pair at p >= 0.996 with
  the first four of the RTT (`CCCC`): a 4-bp GC stem closing an 8-nt loop (3 bp with PBS1 = 19). The RTT starts
  with the reverse complement of the cargo's last four bases, and the cassette ends in a T7 promoter remnant
  (`...TCACTATAGGGG`). The rcPBS2's 5' end is fixed by the nick, so `GGGG` is not a choice.
* **Not length-driven.** Ten random inserts of each length in the same backbone: 5/10 rejected at 24 nt, 2/10 at 60,
  2/10 at 120, 1/10 at 240, 2/10 at 360, **0/10 at 501**. No trend with length (n = 10 each, so only that).
  The signal is local to the junction.
* **The INSERT lever did not run, and it could have helped.** ADR 0022's lever needs a clean in-frame CDS. This
  cargo is one: a 159-codon ORF (GGGGS, FLAG, GSG, P2A, the resistance gene) with no stop inside, ending in ONE real
  stop (`TGA`), followed by a 24-nt tail: that `TGA`, a T7 promoter (`TAATACGACTCACTATAGGG`) and one more G. Read in
  the same frame the promoter adds two more stop codons (`TAA`, `TAG`), which are never translated, but the lever's
  guard counts every stop in the whole cargo's translation, so it returned the original only. (An earlier version of
  this note said the cargo "translates with three stop codons" and that the stem was "in that non-coding tail anyway",
  as if the lever could not have helped. There is one stop, and the lever could have helped.) Measured on real folds
  (`scripts/analysis/probe_orf_only_recoding.py`): 8 of 30 random synonymous recodings of the **ORF alone**, tail kept
  verbatim, have no cross-region pair (the count is 0 or 4), so whether the 4-bp stem is called depends on the ORF's
  codons as well as the tail. Changing only the tail's last codon (`GGG` to `GGC` or `GGA`) or dropping two bases
  gives 0 pairs. The guard is stricter than it needs to be for a very common cargo shape (ORF, stop, non-coding
  tail). **Correction (ADR 0029): that experiment used random synonymous codons, which is not what the lever does.
  The lever's own generator returns no variant at all for this ORF (it fails from 16 residues), so relaxing the
  guard would not have rescued the design either.**
* **The other levers.** PBS1 length changes nothing (it is not at the junction). A different rcPBS2 length would
  (11, 13, 14, 15, 17, 18 give 0 pairs), but at this nick only 12 is inside the design's Tm band (37.4 C; 13 is 42.1
  against a limit of 42.0), so that lever would have to widen a deliberate band. Not built.
* **What clears it.** Dropping the cassette's last two bases (`...TATAGG`) gives 0 pairs; one base leaves a 3-bp
  stem. With those two bases dropped the whole pipeline passes on the real models in 114 s, no stand-in gates: splice
  0.162 on real context (a binding verdict), structure 0, docking **advisory only** because 669 nt is outside the
  calibrated 96-232 nt range (interface ipTM 0.21, not judged; issue #28).

What changed: the failure now says what paired and which folder said so (`cand1: cross-region pairing (ViennaRNA):
rcPBS2[0..3] GGGG x RTT[0..3] CCCC (4 bp, p >= 0.996)`), where it said "RiNALMo cross-region pairing" five times.
The count and the verdict are exactly as before.

Left open, on purpose:

* Relaxing the INSERT lever's guard to recode the ORF up to the first stop and keep the rest verbatim: done in
  ADR 0029, with locks for anything that is not protein and a report of every change. It would not have rescued this
  design, because the variant generator returns nothing for it (ADR 0029).
* Whether a stem this short should reject or only warn. Nothing in this repository measures what a 4-bp stem
  between the rcPBS2 and the RTT does to editing (ADR 0026 makes the same point about scaffold pairing), and the gate
  was accepted in ADR 0024 as a heuristic that removes the worst tail. That is a question for a wet-lab control, not
  for this change.

## References

- `scripts/analysis/probe_splice_gate.py`, `scripts/analysis/probe_structure_gate.py`,
  `scripts/analysis/probe_orf_only_recoding.py`; the failing runs `20261104-145156-8bb93b` and `20261104-145239-a8609d`
  on the DSI box.
- `src/pegrna_pipeline/hpc/scorer_factory.py` (`_PaddedSpliceAIScorer`), `recoded_tjpe_workflow.py` (the walk),
  `recoded_driver.py` (`_attach_splice_context`), `genome_ingest.py` (`fetch_flanking_sequence`).
- ADR 0011 (MaxEntScan primary, SpliceAI advisory), ADR 0020 (indel track alignment). Jaganathan et al. 2019 *Cell*
  176:535.
