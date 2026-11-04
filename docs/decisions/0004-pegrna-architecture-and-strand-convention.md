# 0004. pegRNA architecture and strand convention

Date: 2026-06-18
Status: Accepted. The PBS2 clause ("PBS2 occupies the homology-arm slot", identity copy for a
reverse-strand pegRNA) was superseded by ADR 0014: the pegRNA carries rcPBS2 = RC(PBS2), 5' of the RTT.
See `docs/STRUCTURE.md`.

## Context

The recoded TJ-PE designer assembles a multi-segment pegRNA and reasons about two
antiparallel nicks. Early implementations carried latent, catastrophic bugs from
conflating `complement`, `reverse`, and `reverse_complement`, and from deriving PBS1
incorrectly. These bugs produce a pegRNA that looks correct in a report but cannot
prime or cannot template-jump. The segment order and the strand handling needed to be
locked and made auditable.

## Decision

The authoritative TJ-pegRNA structure, written 5'→3', is
`Spacer · Scaffold · PBS2 · RTT · PBS1 · PolyU`, assembled by
`assemble_tj_pegrna_v2`. PBS2 occupies the homology-arm slot (it captures the
ngRNA-nicked strand); PolyU is the 2-3 U La/PE7 tail.

Strand rules (all sequences stored 5'→3'; antiparallel relationships use
`reverse_complement`, never bare `reverse` or `complement`):

- RTT tail = `reverse_complement(sense_flap)` for a reverse-targeting pegRNA
  (`orient_rtt_tail`).
- PBS2 = identity copy of the captured-strand window for a reverse-targeting pegRNA
  (`orient_pbs2`); the complement-of-a-complement collapses to identity.
- PBS1 = `reverse_complement(genome[S-L:S])` for a reverse-targeting pegRNA (identity
  for forward), derived from the genome at the pegRNA nick — never a slice of the
  spacer (`_derive_pbs1`).
- Default `pegrna_strand="Rv"` to enable C-terminal recode-into-CDS.

Every assembled candidate is run through `simulate_tprt`, which reverse-transcribes
from the oriented RTT and asserts reconstruction of the intended sense allele; any
orientation error raises `TPRTMismatch` and rejects the candidate.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| PBS1 = slice of the spacer | Produces a non-priming, dead pegRNA that passes every other check (high-severity review bug). |
| PBS2 = reverse_complement for Rv | Two reviewers proposed this; the complement-of-complement reasoning plus the `simulate_tprt` self-consistency check show identity is correct for Rv. |
| Trust orientation by code review alone | Strand bugs are silent; an executable self-consistency net (`simulate_tprt`) is required. |

## Consequences

- Positive: orientation errors fail loudly at design time, not in the wet lab.
- Positive: segment order is fixed for both assembly and the RiNALMo region map.
- Negative: `simulate_tprt` is a self-consistency check, not an experimental strand
  proof; a wet-lab-validated golden pegRNA is still needed to fully close the loop.

## References

- Zheng et al. 2023, *Nat Commun* (TJ-PE).
- `src/pegrna_pipeline/hpc/strand_orientation.py`, `tprt_sim.py`, `tjpe_workflow.py`.
- `docs/DESIGN_HISTORY.md` §3.
