# 0029. Only protein-coding bases are ever recoded; every change from the input is reported

Date: 2026-11-05
Status: Accepted. Amends ADR 0022 (the INSERT lever) and corrects ADR 0028's addendum.

## Context

The CTCF cassette ends in a T7 promoter, read in frame 0 after the resistance gene's stop codon. Three things
were wrong or missing, one of them found by checking a claim of ours that turned out to be false.

1. **The lever's guard counted the promoter's own stop codons.** The cargo is one clean ORF of 159 codons with
   ONE real stop (`TGA`), followed by `TAATACGACTCACTATAGGGG`. Read in the same frame the promoter adds two more stop
   codons (`TAA`, `TAG`) that are never translated, but ADR 0022's guard refused any cargo with a stop in its
   translation, so the lever never ran. Relaxing it naively would have started recoding the promoter, whose
   consensus cannot change without breaking it. The requirement: nothing that is not protein is ever recoded, the user
   can say what else must not change, and the pipeline shows exactly what it did change.
2. **Nothing said when a design used a recoded cargo.** The walk returned the winner of a synonymous variant with no
   warning and no field: the pegRNA writes the recoded DNA, and it looked exactly like one built from the cargo the
   user supplied.
3. **The recoder's output was not checkable.** The report said "8 synonymous codon edits", not which bases.

**A correction to ADR 0028's addendum, which this work found.** The addendum said the guard was "the only thing that
kept the lever from rescuing this design", from a measurement that recoded the ORF with random synonymous codons.
That was the wrong experiment. The lever does not use random codons: it uses a best-first CAI search under hard
filters (no CpG, no restriction site, no homopolymer of 5 or more), and for this ORF that search returns **no
variant at all**. Measured on the cassette's own residues, asking for 6 variants:

| first N residues of the ORF | 5 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 24 | all 159 (ZeoR alone: 124) |
|---|---|---|---|---|---|---|---|---|---|---|
| variants returned | 6 | 6 | 6 | 6 | 6 | **0** | 0 | 0 | 0 | 0 |

It fails from 16 residues, at the Gly-rich linker and FLAG; a Leu/Lys/Glu repeat of 60 residues returns 6 at once. The
search explores every higher-scoring sequence that has a CpG or a homopolymer before it reaches a compliant one, and
the cap is hit first. So the lever has been inert for realistic cassettes whatever its guard said, relaxing the guard
alone would not have rescued this design, and the T7 element change below is what did.

## Decision

1. **A lock model, in one place** (`hpc/cargo_spec.py`). A cargo is a sequence plus *locked* bases. A *free codon*
   is a codon of a protein-coding part that is not a stop and touches no locked base (a codon that touches one is
   kept whole). The lever changes free codons to other synonymous codons and nothing else, and checks it
   (`assert_locks_hold`: a recoder bug raises; it does not ship). Locks come from, strongest first:
   1. `cargo_parts`: parts (name, DNA, `locked`), as the web builder and `pegrna validate` produce. When given, nothing
      is guessed. An unmarked part is locked. A part marked coding that is not a clean in-frame CDS is kept exactly,
      with a note.
   2. `cargo_locked`: `[first, last]` nucleotide ranges, 1-based inclusive, that must never change. They apply on top.
   3. Known non-protein elements, found by sequence in either orientation wherever they sit: the T7 promoter with its
      G run (or the C run of its reverse complement). Always locked, even in an unannotated cargo and even inside what
      looks like an ORF: the reverse T7 has no stop in frame, so only recognising it keeps the ORF from running
      through it.
   4. Otherwise the cargo is read from its first base in frame 0 and only the ORF up to the first stop codon is
      coding: when it is the whole cargo (a clean CDS, or a CDS and one terminal stop), or at least 30 codons (a wrong
      reading frame reaches that by chance about one time in ten). Everything after the stop is locked.
   An unannotated cargo therefore gets less recoding than before, never more.
2. **Every change is reported, exactly** (`ChangeReport`, on `result.changes`). Each changed nucleotide: scope (cargo
   or genomic arm), part, 1-based position, codon number, old and new codon, amino acid, old and new base. Also the
   stretches kept exactly and why, and any note. It is the report's "Changes from your input" section, one CSV row per
   nucleotide (`recoded_tjpe_changes.csv`), a `changes` object in the gates JSON, and a table in the web result. A
   recoded cargo is also a warning ("order the DNA in the report, not your original"). The genomic arm's changes are
   listed too (they were a count). The report ends "Everything not listed above is exactly as you specified."
3. **The T7 elements have two G's** (`insert_elements`): `T7-G` is `TAATACGACTCACTATAGG` and `rcT7-G` is
   `CCTATAGTGAGTCGTATTA`, 19 nt. The earlier `GGGG` (21 nt, a multiple of 3) was the `GGG` of the consensus plus a
   padding G; at the 3' end of the CTCF cassette it is the reverse complement of `CCCC`, which paired with the
   rcPBS2's `GGGG` and failed the folding gate on every candidate (ADR 0028 addendum). **The consequence is
   deliberate and visible: the element is not a multiple of 3, so mid-ORF it would read what follows one base off.**
   `frame_shift_after` says so; the builder warns on any such part with coding sequence after it, and the engine
   notes two coding parts read in different frames. At the 3' end of a cargo or after a stop codon, where the CTCF
   cassette has it, nothing shifts.
4. **The builder** (web UI): every part is "recodable" or "kept exactly", defaulted from the library (a protein is
   recodable; a promoter, the Kozak 5'UTR and stop codons are locked) and switchable per part; your own sections are
   locked unless ticked as protein-coding; pasted DNA takes excluded ranges. `/api/motifs` carries the `locked` flag
   (one rule, `library_part_is_coding`, shared with the agent), the design endpoints refuse bad lock settings with the
   engine's own message, and `pegrna validate` prints each part as LOCKED or coding and a paste-ready `cargo_parts`.

## Evidence

Real models on the DSI box (api mode, no stand-in gates), the CTCF cassette with the two-G T7 element
(499 nt):

| run | cargo given as | result |
|---|---|---|
| a | a plain string | every gate passes (splice 0.162 on real context, structure 0 pairs); report: "Cargo: used exactly as you gave it (499 nt, 0 changed)", stop codon (nt 478-480) and T7 promoter (nt 481-499) kept exactly, 8 arm nucleotides listed, `recoded_tjpe_changes.csv` has 8 rows |
| b | `cargo_parts` (ORF recodable; `stop` and `T7-G` locked) | identical design; the same lines under the user's own part names |
| c | the old 4-G tail, as a plain string | fails the folding gate as before, and now says why the lever did not help (below) |

For c the failure ends: "insert lever: no synonymous cargo variant could be generated (the quality filters cannot
be met within the search cap), so the cargo was tried only as given".

Tests: `tests/test_cargo_spec.py` (53), `test_construct_locks.py`, `test_web_cargo_locks.py`, and the T7 tests in
`test_insert_elements.py`. Mutation checks: switching motif recognition off fails 9; making the free-codon choice
ignore locked bases fails 2 (one of them uses a recoder that rewrites every codon it is allowed to, because the real
generator seldom touches the protected codons in these tests).

## Alternatives considered

| Alternative | Why not |
|---|---|
| Lock everything unless annotated | Breaks every existing plain-string config; the auto rule already gives an unannotated cargo less recoding, not more |
| Recognise more promoters now (SP6, T3, ...) | One tuple in `cargo_spec`; add them when a design needs them, with the same both-orientations rule |
| Partial-codon locks (change an unlocked base of a codon that has a locked one) | A codon with a locked base is kept whole: simpler to state, check and report |
| Fix the CAI generator in this change | A separate defect with its own design question (below); the lock model and the report are correct whatever the generator returns |

## Consequences and what is left open

* **The lever still cannot recode a realistic cassette.** A constraint-aware generator (beam or dynamic programming
  over local constraints: CpG at a junction, homopolymer, restriction-site windows) would fix it. Not done here; the
  report and the failure message say so when it matters.
* The 30-codon threshold for an unannotated cargo is a stop-frequency argument, not a measurement on real cargos.
* A design that used a recoded cargo now says so, but "order the DNA in the report" still relies on the reader.
* The legacy `run` modes (`tjpe`, `pe7`, the old `pa`) do not use the lock model; they are not the designer's path.
