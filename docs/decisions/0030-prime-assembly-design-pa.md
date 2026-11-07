# 0030. Prime Assembly (`design-pa`): how it works, what the real runs showed, and what is still untested

Date: 2026-11-07
Status: Accepted. Implements option A+ of ADR 0027 (the design is in
`docs/superpowers/specs/2026-11-04-prime-assembly-design.md`, which this ADR corrects in four places) and completes
ADR 0027's removal of the legacy prototype.

## Context

ADR 0027 found that the Prime Assembly prototype could not pair real PRIDICT output, and on 2026-11-04 the owner chose
**A+**: Prime Assembly as in Liu et al. (*Nature* 2026), with a donor that restores the DNA between the nicks, so the
net deletion is zero. The spec fixed the geometry and the order of the gates and left four choices to the owner. Building
it and running it on real models settled those choices and found three things the spec did not foresee:

1. **The folding gate rejected every pegRNA, on three base pairs the design makes inevitable.** The first real run
   found one guide pair left after the off-target screen and rejected both of its pegRNAs, at every PBS length.
2. **A change to a PAM's last G is not always a break.** SpCas9 still cuts NAG and NGA.
3. **The strict off-target default leaves very little.** On the project's own acceptance case (the CTCF window) it
   leaves one pair of 24.

This ADR records what was built, the decisions, the evidence, and what is still open. The status of the design on
paper is in the spec; its structure, drawn with numbers from a real run, is in `docs/STRUCTURE.md`.

## Decision

**What `design-pa` does** (`hpc/pa_geometry.py`, `hpc/recoded_pa_workflow.py`, `hpc/recoded_pa_driver.py`):

1. **Pairs.** Guide A is on the top strand and cuts at or before the insertion point S; guide B is on the bottom strand
   and cuts at or after it ("PAM-in": both PAMs lie between the nicks, so the two flaps grow toward each other). The nicks
   are 30 to 100 bp apart (hard floor 23) and the two PAM-bearing footprints do not overlap. The donor is
   `arm5 + cargo + arm3`; the product is `genome[:n1] + donor + genome[n2:]`, which is the reference plus the cargo at S
   plus protein-silent changes, so INV-0 holds. The cargo must be at least two flaps long (70 nt at the default flap of
   35): a shorter one is refused, with a pointer to `design-recoded`. The cargo itself is never recoded.
2. **The walk.** Geometry first: every PAM-in pair of the window. Then only the guides that belong to some pair go through
   the sequence rules (`guide_sanity`) and the off-target screen (one batched cas-offinder call). Each surviving pair gets
   its arms (the existing partition: coding runs recoded synonymously, non-coding runs byte-exact), its PAMs broken where
   possible, INV-0 asserted, and, for each side, a flap length chosen (30 to 50 from 35 outward) and a pegRNA per PBS
   length checked with `simulate_side`. The pairs are ranked. Then, best first, each pair is judged by ViennaRNA folding
   of both pegRNAs, MaxEnt and SpliceAI on the product, and Boltz-2 docking (A, then B), and the winner is confirmed by
   `simulate_pa`. Every per-pegRNA verdict is memoized by the pegRNA's sequence, so a pegRNA shared by many pairs is
   folded and docked once. The first pair that passes everything wins; if none does, the failure lists every rejection.
3. **Ranking** is `(advisories, genomic run under the two flaps, span, summed on-target score, cut A, cut B)`: fewest
   known risks first, then the shorter copy of the genome under the flaps, then the shorter span (fewer recoding edits),
   and the cuts last so the order is total and the same in every process.
4. **A PAM the donor brings back.** The paper's replacement removes both PAMs; a restoring donor regenerates them, so the
   guides could cut the product again. Where a PAM lies in coding sequence it is broken by one synonymous change to one of
   its two G's, only if that leaves no PAM. A change to NAG or NGA does not (below), so a C or a T is always preferred to
   an A, and a PAM that can only be taken to NAG or NGA is reported as *weakened*, not *broken*. A PAM that cannot be
   changed (non-coding DNA is kept exactly; or no synonymous change exists) is *regenerated*. Weakened and regenerated are
   advisories: they rank the pair lower and are shown, and neither rejects it. Every edit that lies in a PAM, whoever made
   it, is marked in the changes table; undoing the recode's own edit is not a break.
5. **The folding gate and the nick.** The donor copies the genome at each nick, so the RTT's last three bases are the
   reverse complement of the three bases just right of the nick, which are also the spacer's own last three. The
   spacer-PBS helix therefore runs on into the RTT: spacer base `17 + m` pairs with the RTT's last-but-`m` base, for
   `m` = 0, 1, 2. Those three position pairs, and no others, are not counted by the gate
   (`functional_region_pairing(expected_pairs=)`, `pa_geometry.nick_continuity_pairs`); how many of them folded is
   reported as `nick_pairing_bp`. tj-PE starts its insert at the nick and does not meet them; its gate is unchanged.
6. **TTTT is read on the transcript.** The RTT and the PBS are one RNA, so a TTTT is looked for in `RTT + PBS`: an RTT
   ending in TT followed by a PBS starting TT terminates Pol III although neither part holds TTTT alone. The RTT's 3' end is
   fixed by the donor's end, so a run formed there is the same at every flap length.
7. **The owner's choices** (the spec's four open points): the 3' tail is `UUU`, as in tj-PE; a PAM that cannot be recoded
   warns and ranks lower, and does not reject the pair; a TTTT in an RTT or PBS rejects; a genomic run under a flap above
   8 bp is an advisory and never a rejection. Nothing is recoded except protein-coding bases of the genomic arms, and every
   changed nucleotide is listed (ADR 0029). The cargo is never recoded here, so the lock inputs of ADR 0029
   (`cargo_parts`, `cargo_locked`) are accepted and validated and have nothing to protect.
8. **Where it is exposed.** `pegrna design-pa` writes `recoded_pa_report.txt`, `_components.csv`, `_changes.csv`,
   `_donor.fasta` and `_gates.json`. `/api/design` and `/api/design/stream` take `strategy: "pa"` (default `tjpe`), with
   caps a public server needs: at most 60 pairs to try (`PEGRNA_PA_MAX_PAIRS`) and 5 PBS lengths. The page offers the
   choice ("Inside the pegRNA" or "As a separate donor DNA"), advises but never switches by itself, shows an order sheet
   (both pegRNAs and the donor, drawn to scale, with the PAM of each guide in the genome and in the product), and explains a
   failed run by step. The notes on preparing the donor come from the server, from the paper (PMC12262537): a dsDNA PCR
   donor gave about 50% insertion at AAVS1 in HEK293T; a donor with 3' overhangs (lambda exonuclease, five phosphorothioates)
   gave 20 to 40%, slightly less but more precise; a DNA-PK inhibitor, AZD-7648, about doubled it; PE6c with 35-nt flaps.
   They say that none of these experiments used a donor that restores the DNA between the nicks.
9. **The legacy prototype is removed** (`pa_workflow.py`, `pair_scorer.py`, the `pa` branch of `driver.py`, their report
   and annotator functions and tests, two fixture directories, and `scripts/analysis/probe_legacy_pa.py`, the probe that
   produced ADR 0027's numbers). Git has them all: `git log --diff-filter=D -1 --format=%H -- scripts/analysis/probe_legacy_pa.py`
   prints the commit that removed them, and `git show <that commit>^:<path>` brings a file back. `pegrna run --mode pa` still
   refuses, with a pointer to `design-pa`.

## Evidence

**Real models** (api mode, no stand-in gate; DSI box; CTCF with the 499 nt GGGGS-FLAG-P2A-ZeoR-T7 cassette at S = 147; window
297 nt, 24 PAM-in pairs from 16 guides):

| run | what | result |
|---|---|---|
| first real run | strict default (no off-target site tolerated), before the exemption in 5 | 1 guide failed the sequence rules, 11 the off-target screen; **one pair left (101/164)**, and both its pegRNAs were rejected by folding: A on 1 base pair, B on 2, at every PBS length |
| `probe_pa_levers.py` | the 12 pairs left with off-target tolerance unlimited; flap 30 to 50; PBS 13, 15, 17 | **no pair ever folded clean**: for none of the 12 did any flap length give 0 cross-region pairs |
| `probe_pa_nick_helices.py` | 96 pegRNAs (12 pairs, both sides, flaps 30, 35, 42, 50, PBS 15), counted again independently of the exemption | **230 counted pairs, all 230 the predicted three**; none anywhere else. With the prediction emptied out the same script reports all 230 as "other" and exits 1 |
| after the exemption, strict default | pair 101/164 | design found in **175 s**; docking 0.72 (range 0.68 to 0.76) and 0.61 (0.55 to 0.66) against neutral 0.58 and non-binder 0.33 at 147 nt; MaxEnt 5 (< 8); SpliceAI Δ 0.162 on the real genomic context; folding 0 pairs (nick pairing 1 and 2 bp); 8 recoding edits; both PAMs regenerated; flap B shares a 19 bp run with the strand it displaces |
| one off-target site allowed per guide | 4 pairs; **61/149 won** | **155 s**; docking 0.59 (0.53 to 0.70) and 0.72 (0.65 to 0.77); A's PAM weakened (AGG to AAG), B's broken (CGG to CTA); 14 edits; SpliceAI 0.161 |
| through the deployed web service, and again from a browser on the live site | strict default | the same pair, the same docking values (0.72 and 0.61), a 562 nt donor, 140 s through the service alone and 279 s in the browser run, which shared the server's two slots with a tj-PE design (itself 5.6 min) |

Both pegRNAs of these designs are 147 nt (flap 35, PBS 13); the longest one possible is 166 nt (flap 50, PBS 17). Either way
it is inside the docking gate's calibrated range (96 to 232 nt), so docking is a verdict for any cargo size, where a tj-PE
pegRNA for the same cassette would be 650 to 700 nt.

**The weak PAMs.** NAG: SpCas9 cleaves it at about a fifth of NGG's efficiency (Hsu et al. 2013, Supplementary). In a
human-cell PAM-library screen (PAM-DOSE, Tang et al. 2019, PMC6833347) NAG made up 0.21 to 0.23 and NGA 0.08 to 0.14 of
all cleaved PAM reads at three targets. Wild-type SpCas9 gave more than 20% indels at 17 of 17 NGG sites and 3 of 19 NGA
sites, and at none of 18 NGT or 15 NGC sites (Nishimasu et al. 2018, HEK293FT). So an A in the PAM's last two places leaves
a PAM that still cuts, and a C or a T does not.

**What a reader can check.** About 450 tests cover it: `test_pa_geometry.py` (57), `test_pa_simulator.py` (29, including
the property test that every wrong strand operation on either pegRNA makes it raise), `test_recoded_pa_workflow.py` (214),
`test_pa_report_driver_cli.py` (80), `test_web_pa.py` (41), `test_recoder_determinism.py` (7, identical output across
`PYTHONHASHSEED` values in subprocesses), `test_nick_cut_convention.py` (9), and `test_pa_page_contract.py` (9). The last
reads the page's own regular expressions and numbers out of its TypeScript, builds 24 kinds of rejection with the real walk,
and requires each line to be claimed by exactly one step of the page, the right one, and the form's defaults and limits to
be the engine's: a reworded message or a changed limit on either side fails a test instead of mislabelling a failure. The
page's logic has 151 checks that run straight on its sources (`scripts/qa-frontend-logic.mjs`, a CI job of its own) and
84 browser checks (`scripts/qa-ui-quick.mjs`). A bug was put on purpose into the source and the suite required to fail,
for 53 (the contract), 71 (the form and pipeline logic), 14 (the stream reader), 17 (the components), 7 (the PAM labels),
4 (the diagnostics cap) and 13 (the legacy removal); all were caught, and the survivors found on the way (a check that
iterated the code's own list, a boundary with no case, a PAM span one base too long) each got a test.
`scripts/qa-browser.mjs pa` is the live end-to-end check, and `scripts/audit-production.sh` (31 checks) now covers Prime
Assembly.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Reject a pair whose PAM cannot be recoded | At the strict default the only pair of 24 that survives regenerates both PAMs, so the acceptance case would have no design. It is an advisory, ranks lower, and says so |
| Count NAG and NGA as broken | SpCas9 cuts both (above); calling the guide harmless would be false |
| Exempt the whole spacer-RTT region from the folding gate | It would hide real misfolds. Only the three position pairs the design makes inevitable are exempt, and the probe shows they are all the gate was counting |
| Let PRIDICT propose the pegRNAs | ADR 0027: its designs cannot have the nick spacing a pair needs, and its RTT is the wrong kind |
| Raise the flap length to escape the helix | It cannot: the RTT's 3' end is fixed by the donor's end, so the same three bases meet the spacer at every length (`probe_pa_levers.py`) |
| tevopreQ1 as the 3' tail, as in the paper's epegRNAs | Not built; `UUU` as in tj-PE. The tail is an argument of the engine and a one-line change |
| Multi-fragment donors | A non-goal; the paper chains fragments up to 11 kb, here the donor is `arm5 + cargo + arm3` |

## Consequences and what is left open

* **No restoring donor has been tested at the bench.** The paper replaced the region with non-genomic sequence. Every
  result carries a standing caveat saying so, and the page shows it as a callout of its own. All the verdicts above are
  in-silico verdicts.
* **The strict off-target default is thin on this locus.** One pair of 24 survives, with both PAMs regenerated and a 19 bp
  genomic run under flap B. The page's "Off-target hits allowed" setting is the lever, each tolerated hit is listed on the
  result, and the failure view says when the screen is what emptied the list.
* **A regenerated PAM means the guide can cut the product again.** That is reported, not prevented.
* **A stream that goes silent could not be recovered.** One real-model run through the public hostname stopped delivering
  data about a minute in, after the SpliceAI event. The server went on and wrote a complete result at +140 s. The client
  received nothing at all (not even the keepalive every 15 s) until its own 1900 s read timeout, 31 minutes later. Two later
  runs of the same request, one through the tunnel (170 s) and one on the box (191 s), delivered every keepalive and the
  result, so the path works; one connection did not, and the cause was not found. The page now gives up on a stream that has
  been silent for 90 s and says so (`readDesignStream`). **Not done:** the run id is only sent with the result, so a page
  that loses the stream cannot offer the finished files. Sending the id first, and a way to fetch a run's result by id,
  would make that case recoverable and not only explained.
* **A design takes 140 to 280 s** and the server has two slots; a third concurrent design is refused with 429.
* **The static page's `<title>` and meta description still say "C-terminal tagging"** (the other session's search
  text, left alone here). Prime Assembly is not a tagging method.
* The 500 bp range that the TJ-PE and PE7 report caveats state for TJ-PE was not touched or verified here.

## Corrections

* **NGA.** The commit that added the weak-PAM rule and the comment above `_WEAK_PAM_TAILS` gave "8-16%" for NGA. The source
  (Tang et al. 2019) says 0.08 to 0.14, and as a share of all cleaved PAM reads at three targets, not as a rate relative to
  NGG. The comment now says that.
* **The spec's ranking** listed "fewer restored base pairs (`arm5 + arm3`)" before the on-target score and "the smaller
  span" after it. The donor restores exactly what lies between the nicks, so the two are one number, and the second could
  never break a tie the first had left. The implementation has the one criterion; its order is what the spec's wording
  gave, and the spec now says so.
* **The spec and the module docstring** described the walk as "enumerate guides, drop the ones that fail the rules and the
  screen, then pair". It pairs first and screens only the guides that belong to a pair, which keeps the genome-wide scan
  to the guides that can matter. Both now say that.
* **The spec's step 4** said the RTT and the PBS "have no TTTT", which can be read part by part; the rule is on `RTT + PBS`
  as one transcript (item 6). Its step 6 had no exemption in the folding gate (item 5).
* **The commit message of `501e969`** says the page contract test builds 28 kinds of rejection. It builds 24. (A pushed
  commit is not rewritten.)
* **ADR 0027** promised that the legacy probe would stay reproducible. It does not, in the tree: it was removed with the
  prototype it probed (item 9). ADR 0027 gets a note saying where to find it.
