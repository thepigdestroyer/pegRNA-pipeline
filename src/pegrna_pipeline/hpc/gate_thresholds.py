"""Literature-derived default thresholds for the recoded TJ-PE gates.

These replace the earlier "wet-lab-pending" placeholders with values grounded in
published cutoffs / in-silico conventions. They remain *defaults*: a wet-lab
calibration library can still override them, but they are now rationally
defensible rather than arbitrary. Every value carries a citation in
:data:`GateThresholds.citations`.

Direction of conservatism is chosen per gate's purpose (most are SAFETY gates
that should err toward catching a problem):

* splice  -> high-recall cutoff (catch disruptions)
* docking -> flag falling out of the "confident" band
* importance -> standard buried/confident definitions
* microhomology -> break runs at the efficient-homology onset (recoder still
  minimises the actual run; smaller preferred, mirroring nick distance)

Nick distance is intentionally NOT here: it stays 40-70 bp with smaller
preferred, per the experimental design decision.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

__all__ = ["GateThresholds", "THRESHOLDS", "DOCKING_REFERENCE", "DOCKING_CALIBRATED_RANGE"]

#: Boltz-2 NIM interface ipTM for SpCas9 + RNA, by TOTAL RNA length in nt, as the mean
#: over 5 diffusion samples. Columns: (length, NEUTRAL, NON-BINDER).
#:
#: NEUTRAL    = spacer + scaffold + a neutral (random / shuffled) 3' extension: what a
#:              pegRNA that does nothing wrong scores at this length.
#: NON-BINDER = a SCRAMBLED backbone (plus the same kind of extension): an RNA that has
#:              lost Cas9.
#:
#: THIS TABLE IS DERIVED, NOT TYPED. It is the output of
#:     python scripts/analysis/docking_reference.py --emit
#: on scripts/analysis/data/docking_calibration.jsonl, and tests/test_docking_reference.py fails
#: when the two disagree. To change a number, add measurements; do not edit it. (Scaffold-
#: sequestering extensions and the pipeline's own designs are deliberately NOT in the table:
#: they are what it is tested on.)
#:
#: Why it exists. ipTM falls steadily with RNA length for reasons unrelated to binding quality,
#: so the AlphaFold-Multimer rule-of-thumb floor of 0.6 is out of reach for any tjPE-sized
#: pegRNA: neutral extensions of 194 / 232 nt average 0.50 / 0.42. The informative question is
#: whether a design looks like a neutral extension of ITS OWN length or like a non-binder of
#: that length, and the gap between them (0.17-0.27) is several times the sequence-to-sequence
#: spread. ADR 0024 has the validation: leave-one-out on 90 sequences and a held-out spacer.
#:
#: Independent sequences behind each row (neutral / non-binder): 96 nt 5/5, 126 nt 8/6,
#: 156 nt 9/6, 194 nt 16/12, 232 nt 15/8. It is an empirical null, not a theory; the 232 nt
#: non-binder column is the thinnest. Measure more with scripts/analysis/calibrate_docking_gate.py.
DOCKING_REFERENCE = (
    (96, 0.660, 0.392),
    (126, 0.602, 0.343),
    (156, 0.571, 0.324),
    (194, 0.499, 0.294),
    (232, 0.419, 0.252),
)

#: Inside this total-RNA-length range the reference above is interpolated; outside it
#: the gate cannot judge a design and reports advisory-only rather than inventing a
#: number. 96 nt is spacer + scaffold with no extension, the shortest possible pegRNA.
DOCKING_CALIBRATED_RANGE = (96, 240)


@dataclass(frozen=True)
class GateThresholds:
    # --- splice (SpliceAI/Pangolin Δ) ---
    splice_max_delta: float = 0.2

    # --- AlphaFold structural importance ---
    plddt_min: float = 70.0
    rsa_buried: float = 0.25
    contact_max: int = 18

    # --- scaffold-nCas9 docking (Boltz-2 / AF-multimer ipTM) ---
    docking_min_iptm: float = 0.6
    docking_max_drop: float = 0.2
    # Length-calibrated gate (used whenever the scorer returns per-sample ipTMs): a design
    # passes when it lies at least this far (0.5 = midway) from the length-matched
    # non-binder reference toward the length-matched neutral reference. The legacy
    # absolute-floor / apo-drop pair above remains for scorers that cannot sample.
    docking_samples: int = 5
    docking_pass_fraction: float = 0.5
    # A docking verdict this close to the cutoff (as a fraction of the neutral-to-non-binder
    # span, so 0.15 = 15% of it) is flagged BORDERLINE: still decided by the rule, but reported
    # as not robust. See the citation for where the number comes from.
    docking_borderline_margin: float = 0.15

    # --- microhomology recode ceiling (bp); recoder minimises below this ---
    microhomology_max_run: int = 8

    # --- full-pegRNA 2D base-pair probability decision boundary ---
    # The boundary is a property of the FOLDER, not of the gate: RiNALMo emits a
    # sigmoid score with one calibration, ViennaRNA's partition function another.
    rinalmo_pair_prob: float = 0.5
    rinalmo_min_accessibility: float = 0.5
    vienna_pair_prob: float = 0.99

    # --- PBS2 (second-strand primer-binding site) selection ---
    # PBS2 is genome-determined at the ngRNA nick, so its LENGTH is scanned to hit
    # a Tm target (Biopython nearest-neighbour, default salt/conc) rather than
    # fixed and then rejected. Band + length range below.
    pbs2_tm_lo: float = 37.0
    pbs2_tm_hi: float = 42.0
    pbs2_gc_lo: float = 0.40
    pbs2_gc_hi: float = 0.70
    pbs2_len_min: int = 10
    pbs2_len_max: int = 18

    citations: Dict[str, str] = field(
        default_factory=lambda: {
            "splice_max_delta": (
                "Jaganathan et al. 2019 Cell 176:535 — SpliceAI Δ cutoffs 0.2 (high "
                "recall), 0.5 (recommended), 0.8 (high precision). Safety gate uses "
                "0.2 to maximise detection of splice disruption."
            ),
            "plddt_min": (
                "Jumper et al. 2021 Nature 596:583; Tunyasuvunakool et al. 2021; "
                "AlphaFold DB bands: >90 very high, 70-90 confident, 50-70 low, "
                "<50 disordered. 70 = lower edge of 'confident'."
            ),
            "rsa_buried": (
                "Levy 2010 J Mol Biol 403:660 — interior/core relative ASA < 0.25 "
                "(surface > 0.25); Tien et al. 2013 PLoS ONE max-ASA normalisation."
            ),
            "contact_max": (
                "Yeh et al. 2014 Mol Biol Evol; Shih et al. 2012 — weighted contact "
                "number correlates with mutational sensitivity / evolutionary rate; "
                ">18 heavy-atom contacts <5 A marks a densely packed core residue."
            ),
            "docking_min_iptm": (
                "Evans et al. 2022 (AlphaFold-Multimer); Abramson et al. 2024 Nature "
                "(AF3) — ipTM >0.8 high, 0.6-0.8 medium, <0.6 likely incorrect; "
                "0.6 = floor for a credible interface."
            ),
            "docking_pass_fraction": (
                "Midway between the length-matched NON-BINDER and NEUTRAL references "
                "(DOCKING_REFERENCE): 'does this look more like a neutral extension of its "
                "own length, or like an RNA that has lost Cas9?'. Measured on the Boltz-2 NIM "
                "with 5 diffusion samples per call. Neutral and non-binder differ by 0.17 "
                "(232 nt) to 0.27 (96 nt); the between-sequence sd is 0.03-0.05 for neutral "
                "designs and 0.01-0.07 for non-binders. Leave-one-out over 90 sequences "
                "(each judged by a table built without it): 53/53 neutral accepted, 37/37 "
                "non-binders rejected. On a 0 (non-binder) to 1 (neutral) scale neutral designs "
                "sit at 1.00 +- 0.21 (lowest 0.51) and non-binders at 0.00 +- 0.18 (highest "
                "0.42), so 0.5 separates them, and any cutoff from 0.45 to 0.50 makes no "
                "error. A second spacer never used to build the table was judged correctly "
                "12/12. Single draws are far noisier (sd 0.03-0.11; one scrambled non-binder "
                "draw reached 0.634), which is why 5 samples are averaged. It is a check that "
                "Cas9 can still bind, not of the scaffold's folding: it rejects every "
                "complement of 40+ nt over the scaffold's 5' two-thirds but only 4 of 12 "
                "24-nt complements, and none that target the 3' hairpin (ADR 0024)."
            ),
            "docking_borderline_margin": (
                "0.15 of the neutral-to-non-binder span either side of the cutoff. The closest "
                "calls in the leave-one-out data (lowest neutral 0.51, highest non-binder 0.42 "
                "on the 0..1 scale) both lie inside +-0.15 of the 0.5 cutoff. The band flags "
                "3.8% of neutral designs (still accepted) and 5.4% of non-binders (still "
                "rejected), leaving 96% of verdicts clear. A flagged verdict is decided by the "
                "rule but is not robust, and is reported as such."
            ),
            "docking_max_drop": (
                "Derived from the AF-Multimer confidence bands (one full 0.2-wide "
                "band): a >0.2 ipTM drop vs the apo scaffold-Cas9 baseline indicates "
                "the extension degraded docking."
            ),
            "microhomology_max_run": (
                "RecA/RAD51 minimal efficient processing segment (MEPS) ~8 bp "
                "(Shen & Huang 1986; Hsieh et al. 1992); MMEJ efficiency onset ~6-8 bp "
                "(Sfeir & Symington 2015). Break perfect runs reaching this onset; "
                "the recoder minimises the actual longest run below the ceiling."
            ),
            "rinalmo_pair_prob": (
                "Sigmoid base-pair-probability decision boundary 0.5 for a confident "
                "predicted pair; RiNALMo (Penic et al. 2024) val-tunes 0.01-0.30 for "
                "max F1 on natural families — 0.5 used here to avoid over-blocking on "
                "low-probability pairs."
            ),
            "vienna_pair_prob": (
                "Calibrated against a null, not borrowed from RiNALMo. 300 random "
                "pegRNAs per size (same spacer + scaffold, random rcPBS2/RTT/PBS1), "
                "fraction with ZERO cross-region pairs at p>=0.5/0.8/0.9/0.99: "
                "139 nt 35/62/73/88%; 191 nt 9/37/54/80%; 228 nt 12/37/54/84%. At 0.5 a "
                "gate rejects ~90% of random tjPE-sized designs and so screens out almost "
                "everything; at 0.99 it removes the worst ~15-20% tail. Seeded; "
                "reproducible with scripts/analysis/calibrate_structure_gate.py."
            ),
            "rinalmo_min_accessibility": (
                "A functional region is flagged occluded when mean unpaired "
                "probability < 0.5 (more paired than not) — conservative structural-"
                "interference flag pending efficiency calibration."
            ),
        }
    )


#: Singleton used as the default source of truth across the gates.
THRESHOLDS = GateThresholds()
