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

__all__ = ["GateThresholds", "THRESHOLDS"]


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

    # --- microhomology recode ceiling (bp); recoder minimises below this ---
    microhomology_max_run: int = 8

    # --- RiNALMo full-pegRNA 2D base-pair probability decision boundary ---
    rinalmo_pair_prob: float = 0.5
    rinalmo_min_accessibility: float = 0.5

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
            "rinalmo_min_accessibility": (
                "A functional region is flagged occluded when mean unpaired "
                "probability < 0.5 (more paired than not) — conservative structural-"
                "interference flag pending efficiency calibration."
            ),
        }
    )


#: Singleton used as the default source of truth across the gates.
THRESHOLDS = GateThresholds()
