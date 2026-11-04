# Prime Assembly on the canonical workflow, with the donor restoring the DNA between the nicks — Design Spec

Date: 2026-11-04
Status: Proposed (awaiting user review)
Builds on: `docs/STRUCTURE.md` (the structure this spec relies on, read it first), ADR 0003 (arm partition),
ADR 0004 and 0014 (pegRNA architecture, strand convention), ADR 0021 (one canonical workflow), ADR 0024 to 0026
(calibrated gates, scaffold advisory), ADR 0027 (why the legacy PA prototype cannot work). Issues #29 (this
feature), #30 (recoder is not deterministic), #31 (reverse-strand nick coordinate).

## 1. Problem

The canonical workflow (`design-recoded`, tj-PE) carries the whole insert in the pegRNA's 3' extension. That is
right for tags and wrong for cassettes: past about 25 to 60 bp of cargo the pegRNA exceeds the 232 nt the docking
gate can judge, so the gate degrades to an advisory (41 of the 59 real designs measured for #28), and a 480 bp
cassette would need a 650 to 700 nt pegRNA. The project's own acceptance case for v1.0 is exactly such a cargo:
CTCF + GGGGS-FLAG-P2A-ZeoR (issue #2, M8).

Prime Assembly (Liu et al., *Nature* 2026) moves the insert out of the pegRNA: two extended guides write two
~35-nt flaps, the ends of a separate linear donor match them, and the DNA between the nicks is replaced by the
donor. The pegRNAs stay about 147 nt whatever the insert. The legacy PA code cannot do this (ADR 0027).

On 2026-11-04 the owner chose **A+**: build it as in the paper, and **the DNA lost between the nicks is added
back by the donor**, the way the Replacement Arm does it in tj-PE.

## 2. Goal and non-goals

**Goal.** `pegrna design-pa`, beside `design-recoded`, taking the same locus inputs plus a few PA parameters, and
returning two pegRNAs, one donor sequence, the predicted product allele and a report. Every gate is a real
model. If nothing survives it fails loudly with a per-candidate diagnostic, as tj-PE does.

**Acceptance.**
1. CTCF + GGGGS-FLAG-P2A-ZeoR (480 bp) on the DSI box with real models: a design is found; no stand-in gates;
   both pegRNAs are at most 166 nt; docking is a verdict, not an advisory; INV-0 holds; the simulator reproduces
   the intended product.
2. A property test over random genomes, strands, nick positions and cargos: `simulate_pa` returns the intended
   allele, and every single wrong strand operation (reverse, complement, identity where a reverse complement
   belongs) on either pegRNA's RTT or PBS makes it raise.
3. Stage 3 only: the live site offers it and the browser acceptance test passes with it.

**Non-goals.** Multi-fragment donors (the paper chains fragments with 30 bp overlaps up to 11 kb; here the donor
is `arm5 + cargo + arm3`). PCR primer design. PRIDICT2 efficiency prediction for flaps (an unvalidated
extrapolation of a model trained on short edits). Allele-specific (het) designs. Any change to tj-PE other than
the two stage 0 fixes.

## 3. What the design rests on

**From the paper** (preprint PMC12262537; the supplement with guide-design parameters was not available):
twin epegRNAs with PE6c; 35-nt flaps (30 and 50 robust, 13 nt very weak); the donor's ends match the flaps; the
donor's 5' ends must be displaced or resected for the flaps to anneal (dsDNA PCR product, or a 3'-overhang
version made with lambda exonuclease and five phosphorothioates); the DNA-PK inhibitor AZD-7648 about doubles
insertion; about 50% insertion of 0.8 kb at AAVS1 in HEK293T, junctions precise in 83 to 88% of reads. The region
between the nicks is replaced.

**From the code** (`docs/STRUCTURE.md` §1): the Replacement Arm machinery (`_partition_arm`, `recode_arm`), INV-0
(`assert_inv0`), the in-silico self-check pattern (`simulate_tprt`), the gates, and the walk shape.

**From the CTCF probe** (`scripts/analysis/probe_pa_geometry.py`; S = 147, the stop codon starts at S):

- 24 valid PAM-in pairs with S between the nicks and a span of 30 to 100 bp.
- The 5' arm is coding and recodable at every length tried (6 to 92 nt).
- The 3' arm is the stop codon and UTR: non-coding, so byte-exact. The genomic stretch under flap B is about the
  arm's length, capped at the flap.
- The best pair, A at 101 + B at 149 (span 48 = arm5 46 + arm3 2), has **both PAMs in coding sequence** and
  genomic runs under the flaps of only 8 and 5.
- A 24 bp FLAG cargo gives donors of 59 to 87 nt: the two flaps overlap, which is twinPE, not assembly.

**Two defects found on the way** (stage 0): the recoder is not deterministic across processes (#30), and
`enumerate_nick_sites` reports reverse-strand cuts one base off (#31).

## 4. Design

### 4.1 Geometry: `hpc/pa_geometry.py` (pure, no models)

Conventions are those of `docs/STRUCTURE.md` §0. In particular a cut is `b`, the first base right of the cut,
for both strands, defined here and tested independently of `enumerate_nick_sites`.

- **Guides.** Enumerate NGG-bearing guides on both strands. Guide **A** is on the top strand with cut `n1 <= S`;
  guide **B** is on the bottom strand with cut `n2 >= S`. Left-A, right-B is the only arrangement whose flaps
  fill the gap (PAM-in), so it is the only one built. A pair needs `min_span <= n2 - n1 <= max_span`, S inside,
  and the two PAM-bearing footprints not overlapping (`n1 + 6 <= n2 - 6`).
- **Arms.** `arm5 = genome[n1:S)` and `arm3 = genome[S:n2)`, each built by the existing partition: coding runs
  recoded synonymously with the longest stretch identical to the displaced strand at most `max_run`;
  non-coding runs byte-exact (Decision 15). Where a guide's PAM lies in a coding run, the recode is asked to
  break it too (a synonymous change at the PAM) when one exists.
- **Donor.** `arm5 + cargo + arm3`. **Product.** `genome[:n1] + donor + genome[n2:]`, which is the reference plus
  the cargo at S plus silent changes: INV-0 holds, net deletion 0.
- **Flaps and pegRNAs.** For a flap length `F` and PBS length `L`: guide A has `PBS = RC(genome[n1-L:n1])` and
  `RTT = RC(donor[:F])`; guide B has `PBS = genome[n2:n2+L]` and `RTT = donor[-F:]`. The pegRNA is the standard
  layout `Spacer · Scaffold · RTT · PBS · tail`, which is not tj-PE's (there is no rcPBS2; the second strand is
  captured by the donor).
- **Cargo rule.** `len(cargo) >= 2 * F`. Below that the flaps overlap each other, and the designer refuses and
  says to use `design-recoded`. It also guarantees each flap stays within its own arm plus the cargo.

### 4.2 The walk: `hpc/recoded_pa_workflow.py`

A ranked walk over guide pairs. The sides are nearly independent, since guide A's pegRNA depends on `(n1, S,
cargo)` and guide B's on `(n2, S, cargo)`, but not quite: a guide's PAM can lie inside the *other* guide's arm
(in the CTCF best pair, B's PAM is inside arm5), so breaking it is a pair-level edit to the donor. What keeps the
expensive models cheap is therefore not independence but **memoization: every per-pegRNA gate is cached by the
pegRNA's sequence**, so a pegRNA that appears in many pairs is folded and docked at most once, and a pair whose
edit does not touch a flap re-uses its neighbour's verdicts.

1. **Enumerate** guides on both strands; drop those failing `guide_sanity` (TTTT, GC) and the off-target screen
   (cas-offinder, batched for all guides in one call, tolerance 0 and exact duplicates never tolerated: ADR 0025).
2. **Pair and rank** (below). For each pair, in order:
3. **Cheap, deterministic, per pair.** Build the two arms (coding runs recoded, non-coding byte-exact), break the
   PAMs where a synonymous change exists, assemble the donor and the product, assert INV-0, and check that the two
   donor ends share no run longer than `max_run` (`longest_common_run(donor[:F], donor[-F:])`; if they did, the
   two flaps could pair with each other instead of with the donor).
4. **Per guide, choose the extension.** Scan the flap length `F` from 30 to 50 outward from 35 (35, 36, 34, 37,
   ...): the first `F` whose RTT and PBS have no TTTT (a Pol III terminator) and whose RTT does not start with C
   (it pairs with the scaffold); failing that, the first `F` with no TTTT; failing that, the pair is rejected. Then
   for each PBS length make the pegRNA.
5. **Models on the product, once per pair:** MaxEnt cryptic-splice, SpliceAI, and the AlphaFold importance check
   at S (computed once and reused), all unchanged.
6. **Models per pegRNA, memoized:** ViennaRNA folding (cross-region pairs at p >= 0.99, ADR 0024) with the
   scaffold-pairing advisory (ADR 0026), then Boltz-2 docking against the length-matched reference. A PA pegRNA is
   at most 166 nt (20 + 76 + 50 + 17 + 3), always inside the calibrated 96 to 232 nt, so docking is a verdict.
7. **First pair that passes everything wins**, after `simulate_pa` has confirmed the ordered sequences.

**Ranking** (a heuristic, as there is no wet-lab data for PA here): fewest applicable advisories; then the
smaller combined genomic run under the two flaps; then fewer restored base pairs (`arm5 + arm3`); then the higher
summed on-target score; then the smaller span; then the cut positions, so the order is total and deterministic.
It prefers pairs that avoid the two known risks (a flap that is a copy of the genome, a regenerated PAM).

### 4.3 The simulator: `simulate_pa`

The design rests on strand arithmetic that is easy to get silently wrong (ADR 0004 records how). So, like
`simulate_tprt`, an executable check rebuilds the outcome from what would be ordered and compares it with the
intent:

- from each pegRNA's extension, recover the flap it would write and check that it equals the donor end it must
  match (`flap A == donor[:F]`, `flap B == RC(donor[-F:])`);
- check each PBS pairs with the nicked strand's 3' end at its cut;
- splice the donor between the cuts and require the result to equal the intended product.

It is a self-consistency check, not an experimental proof of the strand logic (ADR 0004 says the same of
`simulate_tprt`).

### 4.4 Outputs, command line, API, UI

- **Result** (`RecodedPAResult`): both spacers, PBS and RTT, both pegRNAs, `F` per side, the donor with a part
  map (`arm5`, `cargo`, `arm3`, the two flap-matching ends), the product, the cuts, gate results, warnings,
  candidates tried, and the number of pairs considered.
- **Order sheet:** pegRNA A, pegRNA B and the donor, 5'→3', and the paper's notes on preparing the donor
  (dsDNA PCR, or a 3'-overhang version; AZD-7648; PE6c), worded as one locus in HEK293T.
- **Standing caveat on every result:** a restoring donor is this project's extension of the published method.
  The paper replaced regions with non-genomic sequence; nothing validates a donor whose arms copy the genome
  under the flaps.
- **CLI:** `pegrna design-pa --config ... --mode api|fake --output-dir ...`, writing
  `recoded_pa_report.txt`, `recoded_pa_components.csv`, `recoded_pa_donor.fasta` and `recoded_pa_gates.json`.
  Example config: `examples/ctcf_cterm/config_pa.json`.
- **API and UI (stage 3 only):** a `strategy` field on `/api/design` (`tjpe` by default, `pa`), and in the form a
  choice of how the insert is delivered: inside the pegRNA, or as a separate donor DNA. The UI is built last,
  after the acceptance run, as the plan for #29 already said.

### 4.5 Retiring the legacy prototype

`pa_workflow.py`, `pair_scorer.py`, the `pa` branch of `driver.py` and their tests (including
`tests/test_pa_real_pridict_shape.py`) are removed when `design-pa` ships. `pegrna run --mode pa` then points to
`design-pa`.

## 5. Defaults and where each comes from

| Setting | Default | Provenance |
|---|---|---|
| flap length `F` | 35, scanned 30 to 50 per side | the paper: 35 used, 30 and 50 robust, 13 very weak |
| PBS length `L` | 13, 15, 17 | tj-PE's PBS1 lengths |
| span `min` / `max` | 30 / 100 bp | the legacy plan's values; hard floor 23 (`pair_scorer`); geometric minimum 12 |
| minimum cargo | `2 * F` (70 bp) | below it the flaps overlap each other, and a flap could reach into the other side's arm |
| `max_run` | 8 | tj-PE's target (RecA MEPS about 8 bp; MMEJ onset 6 to 8) |
| advisory on a genomic run under a flap | above `max_run`, where it cannot be recoded | project value, **unvalidated for PA** |
| folding, scaffold advisory, docking | the existing `THRESHOLDS` | ADR 0024 (folding calibrated at 139, 191, 228 nt), 0026, 0024 (docking 96 to 232 nt) |
| off-target tolerance | 0, exact duplicates never | ADR 0025 |
| pegRNA 3' tail | `UUU` (PolyU) | tj-PE's; see decision 1 |

## 6. Testing (TDD, tests first at each stage)

- **Geometry** (`tests/test_pa_geometry.py`): the cut convention for both strands against an independent
  construction (build the genome around a known protospacer, compute the cut by hand); left-A right-B only; S
  inside; arms, donor, flap, RTT and PBS derivations; the product equals the intended allele.
- **Simulator** (`tests/test_pa_simulator.py`): the property test and the mutation tests of §2, acceptance 2.
- **Workflow** (`tests/test_recoded_pa_workflow.py`, fake scorers): each gate rejects with a diagnostic; each
  pegRNA is docked at most once however many pairs contain it (memoization by sequence); the
  minimum-cargo refusal; **identical output across `PYTHONHASHSEED` values in subprocesses** (this also pins #30).
- **CLI and report:** files, content, the standing caveat.
- **Real models:** `scripts/analysis` acceptance run on DSI (acceptance 1), then the browser test at stage 3.

## 7. Rollout (each stage its own commits, CI, and a DSI redeploy where code changes)

0. **Fixes first.** #30 (deterministic recoder, with the subprocess test) and #31 (one nick convention, with the
   independent test and a forward-strand design checked base by base). Real-model run of the CTCF demo before
   and after, because #30 changes which valid arm some inputs produce.
1. **`pa_geometry.py` and `simulate_pa`** with their tests. No models.
2. **`recoded_pa_workflow.py`, `design_pa`, the CLI and reports**, the example config, the DSI acceptance run,
   ADR 0029, the legacy removal.
3. **API, UI, browser acceptance, docs, deploy.**

## 8. Risks, and decisions that are the owner's

**Risks.**
- **No wet-lab validation of a restoring donor.** This is the largest. The mitigations (recode under the flaps,
  prefer short non-coding stretches, warn) are reasoning, not evidence.
- The paper's guide-design parameters (PBS lengths, spacing) were not in the text available; the defaults come
  from twinPE conventions and from tj-PE.
- Longer 5' arms need more edits (up to about 17 for 91 nt), so the ranking pushes toward short arms.
- Real-model verification needs the DSI box and so Zscaler to be running.

**Decisions for the owner** (defaults in brackets):
1. **3' tail:** PolyU like tj-PE [default], or tevopreQ1 as in the paper's epegRNAs (an option I would add).
2. **A PAM that cannot be recoded** (non-coding stretch): warn [default], or reject the pair.
3. **TTTT in an RTT or PBS:** reject [default, consistent with the spacer rule in ADR 0025], or warn for users
   who order synthetic pegRNAs.
4. **The advisory threshold for a genomic run under a flap** (8): keep, or change once there is bench data.
