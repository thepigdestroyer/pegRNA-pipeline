# 0027. Prime Assembly: why the prototype cannot work on real data, and what a designer needs

Date: 2026-11-03
Status: Accepted as a finding. The decision it asked for was made on 2026-11-04: option **A+**, below.
Design: `docs/superpowers/specs/2026-11-04-prime-assembly-design.md`; structure: `docs/STRUCTURE.md`.

## Context

Issue #29 reported Prime Assembly (PA) as "largely implemented on the legacy HPC path, 42 passing
tests", with the gate cascade, CLI/web exposure and an end-to-end run missing. That was a reading of
the code and the tests. This ADR is what happened when it was run against the real thing.

## Findings

**1. The prototype cannot pair real PRIDICT output.** The experiment is committed
(`scripts/analysis/probe_legacy_pa.py`, run on the DSI box). It builds the two flap windows the way
the legacy driver does (each flap inserted at the same site, the CTCF stop codon), runs real
PRIDICT2, and applies the legacy `filter_compatible_pairs`:

* PRIDICT returned **765 pegRNAs per flap**, and the same `(strand, Editing_Position)` set for both:
  Fw 0, 6, 11 and Rv 2, 17.
* `Editing_Position` is the edit's **distance from the nick**, not a coordinate. The filter reads it
  as a coordinate on a shared axis and needs `|up - down|` of 30-80. The largest possible difference
  is 17, so **no pair can ever pass**: 0 compatible pairs.
* The test fixtures use 50/58 (upstream) against 100/110 (downstream). No real run produces those,
  which is why 42 passing tests said nothing about whether it works.

**2. It is structural, not a threshold.** PRIDICT only proposes pegRNAs whose nick lies within about
17 bp of the edit, so for one insertion site two of them have nicks at most about 28 bp apart. A
Prime Assembly pair needs nicks that flank the region to be replaced. That spacing cannot come from
PRIDICT's insert-at-one-site designs, whatever the filter does.

**3. Its pegRNAs are the wrong kind.** PRIDICT's RT template for "insert this at S" encodes the
genomic gap from the nick to S, then the flap, then genomic homology after S. A Prime Assembly flap
is a stand-alone sequence that must match one end of the donor, so the donor of the prototype
(`flap + bulk + flap`) would not anneal to what such a pegRNA writes.

**4. The rest of the plan's checklist is missing too.** The "no spacer overlap" and "flap
non-overlap" pair rules do not exist; there is no gate (folding, docking, off-target, splice) on
either pegRNA or the donor; and nothing exposes it (the web app serves only `design-recoded`).

## What the published method specifies

Read: the bioRxiv preprint of Liu et al., "Prime assembly with linear DNA donors enables large
genomic insertions" (PMC12262537; published as *Nature* 654:806, 2026, not read in full here).

* twinPE-style: two pegRNAs (epegRNAs) on opposite strands, with an editor such as PE6c.
* Each makes a **3' flap of about 35 nt**. Flaps of 30 and 50 nt worked well and 13 nt "very weakly".
* The two flaps are **not complementary to each other** (as in twinPE); each matches one end of a
  linear donor, which supplies the insert (0.1 to 11 kb).
* The genomic sequence between the two nicks is replaced (a 1 kb region was replaced with a HaloTag
  at HEK3).
* Donors: PCR-made dsDNA, 3'-overhang dsDNA (lambda exonuclease), or ssDNA; several donors can
  overlap by 30 bp and assemble in the cell. A DNA-PK inhibitor (AZD-7648) roughly doubled efficiency.

## What a designer would need

It is a new designer, not a repaired filter, built the way the tj-PE flow already builds its
second nick: **choose the nicks itself and derive PBS and RTT from the genome**, using PRIDICT only
to rank spacers.

1. Enumerate NGG protospacers on both strands whose nicks flank the interval to replace, at a chosen
   spacing; screen them with `guide_sanity` and `offtarget_rejection` (default tolerance 0).
2. For each, a genome-derived PBS and an RTT equal to the reverse complement of a ~35 nt donor-end
   flap, assembled into a standard pegRNA.
3. Gates already built, per pegRNA: folding (ViennaRNA at 0.99, with the scaffold-pairing advisory)
   and the length-matched docking gate. A PA pegRNA is roughly 96 + 35 + 15 = 150 nt, well inside the
   calibrated docking range (96-232 nt).
4. Donor output whose ends are the two flaps, and the integrated product for the splice checks.

## The decision this is waiting on

What is PA meant to do for this project? It changes the design.

* **A. Replacement or insertion where losing the sequence between the nicks is acceptable**, as in
  the paper. Build the designer above.
* **B. Insertion with no net deletion** (for example a C-terminal tag, where the between-nick
  sequence is coding). PA replaces that region, so the donor would have to restore it, and the flaps
  would then copy genomic sequence. That is a different design problem and needs the biology worked
  out first.
* **C. Keep tj-PE for large inserts and drop PA from v1.0.** The tj-PE flow already does the
  zero-deletion case with a recoded arm, and runs real models for cargos the UI offers up to a few
  hundred bp.

Recommendation: C unless there is a concrete PA use case, and A if there is. B should not be
guessed at.

## Decision

Do not expose PA. Mark it a prototype where users meet it (README, and a warning from
`pegrna run --mode pa` before any cluster job is queued). Keep the evidence reproducible
(`scripts/analysis/probe_legacy_pa.py`) and pin the verified facts in
`tests/test_pa_real_pridict_shape.py`, which should be deleted, not "fixed", when a real designer
replaces the prototype.

## Consequences

* Nothing new can be wrongly believed: an empty PA result now comes with its reason.
* The 42 tests of the prototype are unchanged and still describe the prototype's own logic; they
  are not evidence that it works.
* No Prime Assembly feature ships until A, B or C is chosen.

## References

* Liu B. et al., *Prime assembly with linear DNA donors enables large genomic insertions*
  (bioRxiv 2025.06.16.659978; *Nature* 654:806, 2026).
* `scripts/analysis/probe_legacy_pa.py`, `tests/test_pa_real_pridict_shape.py`,
  `src/pegrna_pipeline/hpc/pair_scorer.py`, `pa_workflow.py`, `driver.py`.

## Update, 2026-11-04: the decision

The owner chose **A, with the donor restoring the DNA between the nicks** ("A+"). The designer is specified in
`docs/superpowers/specs/2026-11-04-prime-assembly-design.md`, and the structure it builds on, for tj-PE as built
and for Prime Assembly as designed, is in `docs/STRUCTURE.md`.

Two corrections to the text above:

- Option B was described as "a different design problem" that "needs the biology worked out first". That
  overstated it. Restoring the span between the nicks is the mechanism tj-PE already uses: the Replacement Arm
  recodes coding runs, keeps non-coding runs byte-exact, and leaves a net deletion of 0. A+ reuses it. Two real
  questions remain, and they are narrower: a restored stretch that lies under a flap is a copy of the strand the
  flap displaces, and a restoring donor brings back the PAMs that the paper's replacement removes. The spec
  handles both.
- Under "What a designer would need", item 3 said a PA pegRNA is about 150 nt and inside the calibrated docking
  range. That stands, and is stronger than stated: with a 50-nt flap and a 17-nt PBS the longest possible PA
  pegRNA is 166 nt, so docking gives a real verdict for any cargo size, which tj-PE cannot do past about 25 to
  60 bp.

The legacy prototype (`pa_workflow.py`, `pair_scorer.py`, the `pa` branch of `driver.py`) is removed when the new
designer ships, and `tests/test_pa_real_pridict_shape.py` goes with it, as this ADR said it should.

## Update, 2026-11-07: the prototype is gone

Prime Assembly shipped as `design-pa` ([ADR 0030](0030-prime-assembly-design-pa.md)), and the prototype this ADR describes
was removed: `pa_workflow.py`, `pair_scorer.py`, the `pa` branch of `driver.py`, `tests/test_pa_real_pridict_shape.py`, and
the probe `scripts/analysis/probe_legacy_pa.py` that produced the numbers above. The numbers stay here, which is their
record. The probe is not in the tree any more, so the promise above to "keep the evidence reproducible" is met through git
instead: `git log --diff-filter=D -1 --format=%H -- scripts/analysis/probe_legacy_pa.py` prints the commit that removed
it, and `git show <that commit>^:scripts/analysis/probe_legacy_pa.py` brings it back (it also needs `pair_scorer.py` from
the same parent, and a PRIDICT2 install). `pegrna run --mode pa` still refuses, with a pointer to `design-pa`.
