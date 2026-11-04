# The structure of a tj-PE design and of a Prime Assembly design

Reference, 2026-11-04. Checked against the code (`recoded_tjpe_workflow.py`, `strand_orientation.py`,
`tprt_sim.py`), ADR 0003, 0004 and 0014, and a real run of the CTCF demo. **tj-PE below is what is built.
Prime Assembly below is the design in
[the spec](superpowers/specs/2026-11-04-prime-assembly-design.md): it is not built yet**, and says so.

## 0. Conventions that apply to every design

- The target is a window of the **sense** strand, indexed from 0. Everything is written 5'→3'. Antiparallel
  relationships use `reverse_complement`; bare `reverse` and bare `complement` are never correct.
- **S** is the insertion site: the cargo goes between `genome[S-1]` and `genome[S]`. S never moves.
- A **cut** is written `b`, the index of the first base to the right of the cut: `genome[:b] | genome[b:]`,
  for a guide on either strand.
- SpCas9 nicks the PAM strand 3 bp 5' of the PAM. Prime editing extends the nicked 3' end, so the **flap grows
  from the nick toward the PAM side**.

## 1. tj-PE as built: one pegRNA and one plain nicking guide

Template-jumping prime editing (Zheng et al. 2023, *Nat Commun* 14:3369). The CTCF demo, with S = 147, the
second nick at 101, a reverse-strand pegRNA, a forward-strand ngRNA and a 24 bp FLAG cargo:

```
genome    ──────────┬──── lost span, 46 bp ────┬──────────
                ngRNA nick (101)         pegRNA nick = S (147)

product   ──────────┬── Replacement Arm, 46 bp ──┬── cargo, 24 bp ──┬──────────
                     recoded copy of the lost span                 FLAG, inserted at S
                     same protein; no run of more than 8 bp identical to the genome

pegRNA, 5'→3', 196 nt
  Spacer 20 · Scaffold 76 · rcPBS2 12 · RTT 70 [ RC(cargo) 24 · RC(arm) 46 ] · PBS1 15 · UUU 3
```

**What is lost, and what puts it back.** The flap replaces the DNA between the two cuts: `genome[n:S)` when the
second nick `n` is 5' of S, `genome[S:n)` when it is 3'. That span is not deleted for good. The **Replacement
Arm** re-templates it as a recoded copy of itself, so the product is the reference plus the cargo at S plus
silent changes, and the net deletion is 0 (`deletion_size = 0`, enforced by the INV-0 gate). The arm is built
per run of the span (`_partition_arm`): coding runs are recoded synonymously so that the longest stretch
identical to the displaced strand is at most `max_run` (8 bp); non-coding runs are kept byte-exact (Decision
15). The sense flap is `arm + cargo` for a nick 5' of S and `cargo + arm` for a nick 3' of S
(`_assemble_edit`).

**The second nick.** The ngRNA is a plain guide with no extension; it only nicks, on the strand opposite the
pegRNA, within `[ngrna_distance_min, ngrna_distance_max]` of S (the code default is 40 to 70 bp; the CTCF example
config widens it to 90).

| Segment (5'→3') | What it is | Derived from |
|---|---|---|
| Spacer, 20 | targets the pegRNA's nick at S | PRIDICT2 in auto-spacer mode, or supplied |
| Scaffold, 76 | `STANDARD_SGRNA_SCAFFOLD`, binds nCas9 | constant |
| **rcPBS2**, 10 to 18 | reverse complement of PBS2; sits **5' of the RTT**, right after the scaffold | `RC(sense bases just 5' of the ngRNA cut)`, length scanned to a Tm of 37 to 42 °C |
| **RTT** | the reverse-transcription template | `RC(sense flap)` = `RC(cargo) · RC(arm)` for a 5' arm |
| **PBS1**, 13 / 15 / 17 | primes off the 3'-OH the pegRNA's nick makes | `RC(genome[S-L:S])` for a reverse-strand pegRNA |
| PolyU, 3 | La / PE7 tail | `UUU` |

**PBS2 is not on the pegRNA.** The RT copies rcPBS2 last, so PBS2 = `RC(rcPBS2)` ends up on the 3' end of the new
flap, where it captures the strand the ngRNA's nick exposed. That capture is the "template jump" (ADR 0014).
rcPBS2 must be 5' of the RTT because the RT extends away from PBS1 toward the pegRNA's 5' end, so nothing 3' of
PBS1 is ever copied.

**Self-check.** `simulate_tprt` rebuilds the integrated allele from the designed RTT and fails loudly if any
strand operation is wrong. It is a self-consistency check, not an experimental proof of the strand logic.

## 2. Prime Assembly: the design (not built yet)

Prime Assembly (Liu et al., *Nature* 2026; preprint PMC12262537): **twin epegRNAs write two 35-nt flaps, and the
ends of a separate linear donor match those flaps**, so the donor is stitched in without a double-strand break.
The DNA between the two nicks is replaced by the donor (the paper replaces a 1 kb HEK3 region with a 1.3 kb
Halo-tag). The pegRNAs stay about 147 nt whatever the insert's size; the insert lives in the donor.

**The owner's requirement (A+).** As in tj-PE, the span between the nicks must not be lost: **the donor restores
it**, recoded by the same rules as the Replacement Arm.

```
genome   ── U ──┬───────────── replaced span, n1 … n2 ─────────────┬── D ──
          nick n1: pegRNA A, top strand            nick n2: pegRNA B, bottom strand
          PAM-in: both PAMs lie between the nicks

donor    [ arm5 = recoded copy of genome[n1:S) ] [ cargo ] [ arm3 = genome[S:n2), recoded if coding ]
          first F nt  = flap A                              last F nt = flap B (on the other strand)

product  U + donor + D  =  reference + cargo at S + silent changes only, net deletion 0

each pegRNA, 5'→3'
  Spacer 20 · Scaffold 76 · RTT F · PBS L · tail 3  ≈ 147 nt
  A (top):    PBS = RC(genome[n1-L:n1])   RTT = RC(donor[:F])
  B (bottom): PBS = genome[n2:n2+L]       RTT = donor[-F:]
```

**Why PAM-in is forced, not chosen.** The flap grows toward the PAM side of its nick. For the two flaps to fill
the gap, guide A (top strand) must be on the left and guide B (bottom strand) on the right, so both PAMs sit
inside the span. Derived from the SpCas9 rule above, not quoted from the paper.

**Two consequences of restoring the span.** Restored DNA that lies under a flap is a copy of the strand the flap
displaces, so the flap can re-anneal where it came from instead of on the donor: coding stretches are recoded
(at most 8 bp identical), non-coding stretches cannot be (byte-exact, Decision 15) and are warned about, and the
designer prefers nicks that keep that stretch short. And the paper's replacement removes both PAMs, whereas a
restoring donor brings them back, so a PAM in a coding stretch is recoded where a synonymous change exists and
one in a non-coding stretch is warned about.

**What the CTCF probe found** (`scripts/analysis/probe_pa_geometry.py`, S = 147, the stop codon TGA starts at S):

| | |
|---|---|
| NGG-bearing guides in the 297 nt window | 18 top strand, 24 bottom strand |
| valid PAM-in pairs, S between the nicks, span 30 to 100 bp | **24** |
| 5' arm (coding) | recodable at every length tried (6 to 92 nt); 0 to about 17 silent edits, never a run over 8 |
| 3' arm (stop codon and UTR) | non-coding, so byte-exact; the genomic run under flap B is about the arm's length, capped at F |
| best pair | A at 101 + B at 149: span 48 = arm5 46 + arm3 2; **both PAMs in coding sequence**; runs under the flaps 8 and 5 |
| a 24 bp FLAG cargo | donor 59 to 87 nt: the flaps overlap, which is twinPE, not assembly. PA needs a cargo of at least 2 × F |
| the GGGGS-FLAG-P2A-ZeoR cargo (480 bp) | donors of 515 to 543 nt |

## 3. tj-PE and Prime Assembly side by side

| | tj-PE | Prime Assembly |
|---|---|---|
| guides | 1 pegRNA + 1 plain ngRNA | 2 pegRNAs, both extended |
| where the insert is | in the pegRNA's RTT | in a separate donor DNA |
| pegRNA length | 196 nt for FLAG; 650 to 700 nt for a 480 bp cargo | about 147 nt for any cargo |
| docking gate (calibrated to 232 nt) | advisory above about 25 to 60 bp of cargo | a real verdict at any cargo size |
| the lost span | re-templated by the Replacement Arm, inside the pegRNA | restored by the donor's arms |
| what the flap is | arm + cargo | the donor's end |

## 4. Statements elsewhere in the repo that this supersedes

- "PBS2 sits on the pegRNA" / "PBS2 is the identity copy for a reverse-strand pegRNA" (ADR 0004,
  `DESIGN_HISTORY.md` §3, the README's structure line). ADR 0014 corrected it: the pegRNA carries **rcPBS2 =
  RC(PBS2)**, 5' of the RTT.
- "tj-PE loses no DNA". It does lose the span between the nicks; the Replacement Arm restores it.
