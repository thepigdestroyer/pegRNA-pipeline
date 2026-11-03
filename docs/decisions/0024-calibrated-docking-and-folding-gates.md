# 0024. Calibrate the docking and folding gates against controls; fix the arm-side pool

Date: 2026-11-03
Status: Accepted

## Context

Once the off-target and docking gates ran real models (rather than stand-ins), the
shipped real-models demo returned **"No design"**. Two things were asked: why, and are
the gates necessary or just too harsh?

The honest way to answer the second question is to measure each gate against inputs
whose correct answer is known, instead of arguing from the code. This ADR records what
was measured and what changed. Every number below comes from a script in
`scripts/analysis/`.

## Findings

### 1. The demo failed on a pool-construction bug, not on a gate

The run reported `no candidate survived all gates after **0** tried`. Zero candidates
reached any gate. In the auto-spacer path PRIDICT's nicking guides are screened and the
survivors are injected as the *entire* nick pool; only an empty list falls back to the
genome-enumerated nicks. But the arm-side rule (`arm_side=5prime` needs a nick 5' of the
insertion site) was applied later, inside the walk. Real off-target screening removed
PRIDICT's 5' guides, leaving one guide 3' of S. The injected list was therefore
non-empty, no fallback fired, the walk discarded that one guide, and a clean nick (101,
one genome match, no near-matches) was never considered.

Fix: arm side is part of "admissible" when the pool is built
(`nick_search.nick_usable_for_arm`, used by `build_pridict_ngrnas`), so an unusable
PRIDICT list is an empty list. The regression test fails with `NoValidPegRNA` on the old
code and passes now.

### 2. Off-target screen (ngRNA): necessary, and not the cause

Whole-genome cas-offinder (hg38, up to 3 mismatches) on the six candidate nick guides:

| nick | guide | hits at 0 / 1 / 2 mm | verdict |
|---|---|---|---|
| 101 | `GCTGCCACAGATGCCCCCAA` | 1 / 0 / 0 | keep |
| 194 | `TTCTCTGGGCTGTGTTTAAA` | 1 / 0 / 0 | keep |
| 70 | `CGCAGAGGGAGAGGAAGAGG` | 1 / 2 / 23 | drop: 25 near-matches (low-complexity) |
| 232 | `TTCTCCCTTCTTTCTTTTTT` | 1 / 1 / 33 | drop: 34 near-matches, and a poly-T run |
| 67 | `GCCCGCAGAGGGAGAGGAAG` | 1 / 0 / 5 | drop |
| 61 | `TGCTGAGCCCGCAGAGGGAG` | 1 / 0 / 1 | drop: a single 2-mismatch site |

The rule discards the guides it should (dozens of near-matches). It is strict at the
margin: a random guide has ~0.6 expected hits at <=2 mismatches by chance alone, so
"zero" removes roughly half of ordinary guides, and nick 61 fails on one hit. For a
nick-only guide a small tolerance is defensible, but ADR 0010 chose zero deliberately and
nothing here shows it causing harm, so the default is unchanged. It is the one gate where a
case for relaxing exists; it is not what failed the demo. (Follow-up, ADR 0025: the tolerance
is now an explicit, configurable, always-disclosed setting, default still zero. It also found
that the existing `max_offtargets` parameter could never take effect above 0.)

### 3. Spacer sanity is applied inconsistently (resolved in ADR 0025)

The PRIDICT path rejects any spacer with a homopolymer run >= 5, which includes the
`CCCCC` in nick 101's guide; the enumerated path never applies the check. Only `TTTT` (a
Pol III terminator) is a hard biological rule. The same guide is accepted by one path and
rejected by the other. ADR 0025 puts both paths on one definition: `TTTT` and GC are hard
rules, a homopolymer is an advisory.

### 4. Folding gate: 0.5 was too harsh for ViennaRNA

The gate counts cross-region base pairs with probability >= threshold. The default 0.5 is
`THRESHOLDS.rinalmo_pair_prob`, a pair-calling boundary for a different folder. For the CTCF
example's spacer the pass rate was **0% of 324 candidates at 0.5 and 100% of 54 at 0.99**.
The null (300 random pegRNAs per size, same spacer and scaffold, random extension):

| size | p>=0.5 | p>=0.8 | p>=0.9 | p>=0.99 |
|---|---|---|---|---|
| 139 nt (typical PE) | 35% | 62% | 73% | 88% |
| 191 nt (tjPE) | 9% | 37% | 54% | 80% |
| 228 nt (tjPE) | 12% | 37% | 54% | 84% |

At 0.5 the gate rejects ~90% of random tjPE-sized designs, so it screens out almost
everything rather than the worst of them. At 0.99 it removes the worst 15-20% tail.

### 5. Docking gate: the absolute floor was invalid

Boltz-2 NIM, SpCas9 + RNA, interface ipTM, mean of 5 diffusion samples. This is the FIRST
measurement (1-5 sequences per length); the Addendum below extends it to 122 measurements and
supersedes these numbers, which it found to be close.

| total RNA | neutral extension | designed pegRNA | scaffold-sequestering (positive ctrl) | non-binder |
|---|---|---|---|---|
| 96 nt | 0.693 (canonical sgRNA) | | | 0.448 |
| 156 nt | 0.615 (n=3) | | 0.325 | |
| 194 nt | 0.480 (n=5) | **0.544** | 0.276 | 0.272 (n=2) |
| 232 nt | 0.432 (n=5) | **0.432** | | 0.262 (n=2) |

- ipTM falls with RNA length regardless of design quality, so the AlphaFold-Multimer floor
  of 0.6 is unreachable for any tjPE-sized pegRNA. Neutral extensions of 194/232 nt average
  0.48/0.43.
- The designed pegRNAs sit at or above neutral extensions of the same length. They were
  rejected for being long, not for being bad.
- The gate has real sensitivity: an extension built to sequester the scaffold scores
  0.28-0.33, 0.2-0.3 below neutral, and the folding gate deliberately ignores every pair touching
  the scaffold. It should be kept. (The Addendum qualifies this: it sees extensive sequestration,
  not localized pairing, so it is not "the only check that can see" the failure in general.)
- A single sample is too noisy to threshold: sd 0.03-0.11. A scrambled non-binder drew
  0.634 once; one draw of a 194-nt design reached 0.628.
- The service is **deterministic**: identical requests return identical scores (sd 0.000
  over replicates). Earlier text in this repo called docking "stochastic" and said the same
  inputs could pass on another run. That was wrong. `seed` is not accepted by the API (HTTP
  422); the sampling noise is real but fixed per request, so a verdict is reproducible while
  still being one draw from a noisy distribution unless several are averaged.
- The pipeline's own apo baseline (scaffold only, 76 nt) averages 0.632 with samples from
  0.596 to 0.680: it straddles the floor it was supposed to be compared against.

### 6. An outage was reported as a verdict

The NIM client retried a 429 twice (5 s, 10 s), then raised an ordinary exception, which the
pool walk caught per candidate and skipped. A rate limit therefore surfaced as "no candidate
survived all gates". (A calibration burst got 21 of 26 calls answered with HTTP 429.) This is
the mirror image of a stand-in being mistaken for a pass: here an outage is mistaken for a
fail.

### 7. The other gates

INV-0 and the recode are invariants of the design, not heuristics. The cryptic-splice (MaxEnt)
and SpliceAI checks are cheap and pass comfortably on the example (MaxEnt -2.1 against a
threshold of 8; SpliceAI 0.06-0.16 against 0.2). The AlphaFold importance gate is close to
constant for a C-terminal tag (the terminus is exposed). None of these is harsh.

### 8. A design could take the whole machine

Running the demo on the shared 64-core host, one PRIDICT2 worker ran at ~36 cores (3,607% CPU,
355 threads) and had used 280 CPU-minutes after five minutes; the load average reached 63.
`--cores 3` limits PRIDICT's worker *processes*, but PyTorch and TensorFlow start a thread per
core inside each. Oversubscription also made it slower, not faster (it was still running after
nine minutes; capped, it finishes in under 45 s). Separately, one whole-genome cas-offinder scan
took 100 s of wall time but 69 CPU-minutes. With a public site, every visitor's design would do
this to a machine other people use.

## Decision

1. Arm side is evaluated when the candidate pool is built (finding 1).
2. The ViennaRNA scorer carries its own calibrated default, `vienna_pair_prob = 0.99`. The
   workflow uses an explicit value if given, else the scorer's, else the generic 0.5.
3. Scorers that return per-sample ipTMs are judged by a **length-calibrated** docking gate: a
   design passes when its 5-sample mean interface ipTM lies at least halfway from the
   length-matched non-binder reference to the length-matched neutral reference
   (`gate_thresholds.DOCKING_REFERENCE`, linear interpolation). The apo baseline is no
   longer used, saving one NIM request per design. Scorers without per-sample output keep the
   legacy absolute-floor gate.
4. Outside the calibrated length range (96-232 nt; originally typed as 240, now derived from the
   data, see the second addendum) the docking gate is **advisory**: it passes,
   says so in the message, sets `metrics['advisory']`, and the design carries a warning. A
   safety gate that cannot be validated must say so rather than silently accept or reject.
5. Everything a design starts is resource-capped: the pipeline process gets `PEGRNA_THREADS`
   (default 4) for its TensorFlow/BLAS pools and starts under `nice -n 10`; PRIDICT2 workers get
   `PRIDICT2_THREADS` (default 2) each; cas-offinder gets `PEGRNA_CASOFFINDER_THREADS` (default
   16) via `POCL_MAX_PTHREAD_COUNT`. Off-target results are cached per guide (keyed on the genome
   file, PAM and mismatch settings), so repeat runs of the same guides do not rescan. Measured
   after the change: pipeline at nice 10, cas-offinder at ~5 cores instead of ~42, load average
   3.5 instead of 63.
6. A docking-service failure is a `DockingServiceError`: retried with `Retry-After`-aware
   backoff, and if it persists the run fails as what it is. It is never recorded as a
   candidate rejection. The calibrated gate also refuses to run if the service stops
   returning `pair_chains_iptm_scores`, since the references were measured on that metric.
   Each call has a 300 s socket timeout and a 600 s overall budget (a request once hung for
   over 30 minutes). `NVCF-POLL-SECONDS` is deliberately not sent: on this endpoint it turns a
   slow job into a 504 whose request id is not pollable (the status route answers 404).
7. A docking verdict within 15 points of the cutoff is flagged BORDERLINE (Addendum): still
   decided by the rule, reported as not robust, and a borderline pass adds a warning.

## Out-of-sample check

The references above were fitted to 26 five-sample measurements and the gate then
reproduces them (`tests/test_docking_gate_calibration.py`), which is circular. To get an
honest estimate, 16 **new** sequences were scored whose random seeds, shuffles and designs
were not used to derive anything (`scripts/analysis/`, 5 diffusion samples each), and judged
by the implemented gate:

| class (new sequences) | n | outcome |
|---|---|---|
| neutral random extension (194 and 232 nt) | 9 | **9 accepted**, 0 false rejects |
| scrambled backbone + extension (non-binder) | 4 | **4 rejected**, 0 false accepts |
| extension built to pair with the scaffold | 3 | **2 rejected, 1 accepted** |

The accepted one is `RC(scaffold[30:76])` + pad: it scored 0.398 against a cutoff of 0.377
(a margin of +0.021). It pairs with the scaffold's 3' end. Whether that actually impairs
binding is not known: the 3' stem loops tolerate disruption far better than the
repeat:anti-repeat core, so the label "should be rejected" on that design is itself an
assumption, and the gate may be right. The two designs pairing with the core region, and
all four scrambled backbones, were rejected.

Read this with its size in mind. Zero false rejects in 9 neutral designs is compatible with
a true rate as high as ~28% (95% upper bound); pooling the 11 in-sample neutral sequences
gives ~14%. Margins are thin at the longest length: the lowest neutral design at 232 nt
cleared its cutoff by 0.021, and at 232 nt neutral and non-binder references are only 0.17
apart. The honest summary is "no evidence of a high false-reject rate, from a small sample,
with the gate working least well where the signal is smallest".

## Addendum: extended calibration (issue #23)

The first reference rested on 1-5 sequences per length and had one open miss. It was extended to
**122 five-sample measurements**. The raw records are committed
(`scripts/analysis/data/docking_calibration.jsonl`), the table is now *derived* from them
(`scripts/analysis/docking_reference.py`), and `tests/test_docking_reference.py` fails if
`DOCKING_REFERENCE` stops matching them or if any figure below stops being true.

**The table** (neutral / non-binder; independent sequences behind each):

| total RNA | neutral | non-binder | n |
|---|---|---|---|
| 96 nt | 0.660 | 0.392 | 5 / 5 |
| 126 nt | 0.602 | 0.343 | 8 / 6 |
| 156 nt | 0.571 | 0.324 | 9 / 6 |
| 194 nt | 0.499 | 0.294 | 16 / 12 |
| 232 nt | 0.419 | 0.252 | 15 / 8 |

The first table was close. The largest changes: the 126 nt non-binder, which had been
interpolated at 0.386, measures 0.343; and the 96 nt non-binder moves from 0.448 to 0.392. Every
earlier fixture is judged the same way.

**Validation.**

* *Leave-one-out over 90 sequences* (each judged by a table rebuilt without it): **53/53 neutral
  accepted, 37/37 non-binders rejected**. Zero errors in that many judgements bounds each error
  rate below about 6% and 8% respectively (95%, rule of three), not at zero.
* On a scale where 0 is the non-binder reference and 1 the neutral one, neutral designs sit at
  1.00 +- 0.21 (lowest 0.51) and non-binders at 0.00 +- 0.18 (highest 0.42). Any cutoff from 0.45 to
  0.50 makes no error; 0.40 and below accepts 2 of 37 non-binders, 0.55 rejects 1 of 53 neutral
  designs (and 0.70 rejects 4). **0.5 is kept**: it is the midpoint, and the data give no reason to
  move it.
* *A second spacer* (`TTCTCTGGGCTGTGTTTAAA`, the CTCF example's nick 194, AT-richer than the first)
  was never used to build the table: **6/6 neutral accepted, 6/6 non-binders rejected**.
* The pipeline's own designs sit at 1.22 and 1.08 on that scale.

**Borderline verdicts.** The closest calls on either side (lowest neutral 0.51, highest non-binder
0.42) are within 0.15 of the cutoff, so a verdict inside that band is not robust. The gate now says
so: BORDERLINE in the message and the metrics, and a warning on a borderline pass. It does not
change what the rule decides. The band flags 3.8% of neutral designs (still accepted) and 5.4% of
non-binders (still rejected); 96% of verdicts are clear.

**What the gate sees of scaffold sequestration** (this settles the 3'-end question). Extensions
reverse-complementary to a segment of the scaffold, 24 nt at six positions, two sequences each, plus
four longer complements:

| complement of scaffold[a:b] | ipTM position (0 = non-binder, 1 = neutral) | verdict |
|---|---|---|
| 0:24 (repeat, loop) | 0.54, 0.46 | accepted, rejected |
| 10:34 | 0.61, 0.60 | both accepted |
| 20:44 (anti-repeat, nexus) | 0.46, 0.30 | both rejected |
| 30:54 (stem loop 1-2) | 0.64, 0.67 | both accepted |
| 40:64 | 1.07, 0.43 | accepted, rejected |
| 52:76 (3' hairpin) | 0.93, 0.98 | both accepted |
| 0:40, 0:60 (x2), 10:70 (40-60 nt) | 0.05, 0.00, -0.09, -0.08 | all rejected |
| 30:76 (46 nt, 3' end) | 0.51 | accepted |

So the gate rejects every large complement over the scaffold's 5' two-thirds and 4 of 12 24-nt
complements, and it is blind to the 3' hairpin. The earlier "miss" (`RC(scaffold[30:76])`) is the
model's insensitivity there, not a threshold problem: lowering the cutoff far enough to catch it
would reject normal designs. The earlier claim that this gate is "the only check that can see"
sequestration was too strong. It is a check that Cas9 can still bind. It is not a check of the
scaffold's folding.

That leaves a real gap, and the data point at the tool for it: in ViennaRNA's MFE structure,
30 neutral designs pair at most 13 scaffold nucleotides with the extension, and **all 12** sweep
constructs pair at least 22. Whether pairing at that level lowers editing is a biological question
this work cannot answer (there is no wet-lab control), and a gate that rejects on an unvalidated
signal repeats the mistake this ADR is about. So no rejecting rule is added; the gap and the
evidence are filed as an issue, with an advisory warning as the cautious first step.

## Second addendum: where the gate stops working, and what real designs do (issue #28)

The first addendum left the long end and the real-design question open. Both are now measured:
111 more measurements brought the file from 122 to 233 records (`scripts/analysis/data/`).

**The calibrated range is a finding, and it is 232 nt, not the 240 that had been typed.** Rows were
measured to 336 nt (6 neutral and 6 non-binder sequences each at 256, 296 and 336 nt), and the 232 nt
row was raised from 8 to 16 non-binders (23 neutral).

| total RNA | neutral | non-binder | gap | gap in between-sequence sd | n |
|---|---|---|---|---|---|
| 194 nt | 0.499 | 0.294 | 0.205 | 4.3 | 16 / 12 |
| 232 nt | 0.422 | 0.255 | 0.167 | 5.2 | 23 / 16 |
| 256 nt | 0.397 | 0.250 | 0.147 | **3.1** | 6 / 6 |
| 296 nt | 0.347 | 0.235 | 0.112 | 4.6, but a leave-one-out error | 6 / 6 |
| 336 nt | 0.287 | 0.232 | 0.055 | **2.4**, a leave-one-out error | 6 / 6 |

Boltz-2's ipTM keeps falling with length for normal and non-binding RNAs alike, and the two curves
converge until the gap is smaller than the model's own sequence-to-sequence noise. With the cutoff
midway a row is only worth trusting when the gap is at least 4 standard deviations (each class then
sits 2 sd from the cutoff). The range is therefore *derived*: it ends at the last row that meets
that and makes no leave-one-out error. `docking_reference.py --emit` produces it and
`tests/test_docking_reference.py` checks it, as for the table.

Inside the range the evidence got stronger: **61/61** neutral designs accepted and **45/45**
non-binders rejected, leave-one-out (up from 53 and 37), the same closest calls as before (lowest
neutral 0.51, highest non-binder 0.42), and any cutoff from 0.45 to 0.50 still makes no error. The
second spacer is still 12/12. The borderline band flags 3.3% of neutral designs and 4.4% of
non-binders. Past the range 2 of 36 sequences are misjudged, shown for the record and used nowhere.

**Real pipeline designs, measured.** The reference was built from random and scrambled extensions,
so 59 pegRNAs the pipeline itself produces were measured (CTCF: 8 cargos from the UI's motif
library, 2 spacers, all four 5' nicks; 194-316 nt).

* **Inside the range (18 designs): all 18 accepted**, none borderline. They sit at 0.89 to 1.68 on
  the 0 (non-binder) to 1 (neutral) scale, median 1.36. Real designs score **above** the
  random-extension baseline, so that baseline is conservative for them, not flattering.
* The first designs production would try (walk order): 9 of 9 accepted.
* **The old 0.6 floor would have passed 1 of the 59.** This is the measured answer to whether it
  was too harsh.
* **Beyond the range (41 designs): reported, not judged.** For information only, against the
  references of their own length they sit at 0.78 to 1.83, none below 0.5. That is not a verdict;
  it is a reason not to expect trouble.

**The consequence to say plainly.** Most realistic designs fall outside the range: 41 of 59. A
cargo up to about 60 bp stays inside with the nearest second nick, and the farther nicks the 40-90 bp
window allows cut that to about 25 bp (a 24 bp FLAG tag fits almost everywhere). A longer cargo gets
an advisory in place of a docking verdict. That is the honest state of this model, not a setting to
tune; a fixed-length truncation probe (docking only the part of the extension next to the scaffold)
is a possible way to recover a verdict and is recorded in #28 as an experiment, not built.

**Found on the way: the walk evaluated each design up to six times, and `rtt_aa` does nothing.** The
CTCF pool of 18 candidates (6 synonymous-codon "RTT variants" x 3 PBS lengths) is 3 distinct
pegRNAs, because the walk uses only a variant's PBS length. Five different anchor residues give
byte-identical pegRNAs. A design that failed therefore re-ran the docking gate for every repeat,
each a ~40 s paid request, and the UI claimed the pipeline "re-encodes" the anchor residue with a
synonymous codon, which it never did. The walk now evaluates each distinct design once, `rtt_aa` is
optional, and the UI says it does not change the design (`tests/test_walk_evaluates_each_design_once.py`).

## Consequences and limitations

- **The docking reference is an empirical null from 90 sequences**, derived from committed data
  and tested leave-one-out and on a held-out spacer. It is no longer a table typed from a few
  measurements, but it is still an estimate: the 232 nt non-binder column rests on 8 sequences and
  is the thinnest. Extend it with `scripts/analysis/calibrate_docking_gate.py` (it resumes, and
  skips anything already measured).
- It was measured with one scaffold and two spacers, on synthetic extensions (random or
  shuffled). Real RTTs have structure of their own. It says nothing about wet-lab editing
  efficiency: ipTM is a model's confidence in an interface, not an affinity. Neutral designs spread
  by about a fifth of the neutral-to-non-binder span, which is why borderline verdicts exist.
- Past 232 nt the gate reports and does not decide. This was first a guess (the dynamic range was
  0.17 at 232 nt and unmeasured further out); it was then measured to 336 nt, and the second
  addendum shows the gap closing into the model's own noise, so the range is a finding, not a
  placeholder.
- Findings 2 and 3 (ngRNA off-target tolerance, inconsistent spacer sanity) were left
  unchanged here and are settled in ADR 0025.
