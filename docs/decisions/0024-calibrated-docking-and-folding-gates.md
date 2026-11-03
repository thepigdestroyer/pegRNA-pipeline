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
nothing here shows it causing harm, so it is unchanged. It is the one gate where a case
for relaxing exists; it is not what failed the demo.

### 3. Spacer sanity is applied inconsistently (noted, not changed)

The PRIDICT path rejects any spacer with a homopolymer run >= 5, which includes the
`CCCCC` in nick 101's guide; the enumerated path never applies the check. Only `TTTT` (a
Pol III terminator) is a hard biological rule. The same guide is accepted by one path and
rejected by the other.

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

Boltz-2 NIM, SpCas9 + RNA, interface ipTM, mean of 5 diffusion samples:

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
  0.28-0.33, 0.2-0.3 below neutral. The folding gate deliberately ignores every pair touching
  the scaffold, so this is the only check that can see that failure. It should be kept.
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
4. Outside the calibrated length range (96-240 nt) the docking gate is **advisory**: it passes,
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

## Consequences and limitations

- **The docking reference is an empirical null from small samples** (n = 1-5 sequences per
  length). It was fitted and then tested on the same data, so its fit to those sequences
  proves little; the holdout above is the honest estimate. Extend it with
  `scripts/analysis/calibrate_docking_gate.py --per-length 8` before relying on a result near
  the decision boundary.
- It was measured on four extension lengths with one spacer and one scaffold, and on
  synthetic controls. It says nothing about wet-lab editing efficiency. ipTM is a model's
  confidence in an interface, not an affinity.
- The non-binder reference at 126 nt is interpolated, not measured.
- Past 240 nt the dynamic range between neutral and non-binder is narrow (0.17 at 232 nt) and
  unmeasured further out, which is why the gate steps down to advisory rather than extrapolate.
- Off-target zero-tolerance (finding 2) and the inconsistent spacer sanity check (finding 3)
  are unchanged and remain open questions.
