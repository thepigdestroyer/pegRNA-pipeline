"""Regression: the shipped real-models demo failed with "0 candidates tried".

With real off-target screening, the only PRIDICT nicking guide that survived for
the CTCF example sat 3' of the insertion site. The example uses arm_side=5prime,
which can only use a nick 5' of the site, so the walk threw that single guide away
and tried nothing -- while a perfectly good genome-enumerated nick (nick 101, one
exact genome match and no near-matches) was never considered, because a non-empty
injected list suppresses the fallback.

The fix is at pool construction (arm side is now part of "admissible"), so an
unusable PRIDICT list is an EMPTY list and the fallback engages.
"""

import json
from pathlib import Path

from pegrna_pipeline.hpc.data_types import PegRNADesign
from pegrna_pipeline.hpc.nick_search import locate_guide
from pegrna_pipeline.hpc.recoded_driver import design_recoded

CFG = json.loads(Path("examples/ctcf_cterm/config_autospacer_flag.json").read_text())
SITE_S = CFG["site_s"]


class _PridictOffersOnlyADownstreamNick:
    """A PRIDICT stand-in whose only nicking guide is 3' of S (distance +56)."""

    def run_batch(self, rows):
        out = {}
        for row in rows:
            peg = PegRNADesign(
                sequence_name=f"{row['sequence_name']}_pegRNA1",
                spacer="TCAGCATGATGGACCGGTGA",
                pbs="TCGACGCTTAGGAACAC",
                rt="ACGTACGTACGTACGTACGTACGTACGTAC",
                pbs_length=13,
                rt_length=30,
                target_strand="Fw",
                editing_position=50,
                pe2_hek293t=0.7,
                pe2_k562=0.7,
                pegrna_full="",
            )
            nicking = [
                {
                    "Nicking-Protospacer": "TTCTCTGGGCTGTGTTTAAA",
                    "Nicking-Position-to-edit": "56",
                    "Nicking-PAMdisrupt": "Nicking_PAM_disrupt",
                    "Target_Strand": "Fw",
                    "DeepCas9score": "0.7",
                }
            ]
            out[row["sequence_name"]] = {"pegrnas": [peg], "nicking": nicking}
        return out


PRIDICT_GUIDE = "TTCTCTGGGCTGTGTTTAAA"  # the real PRIDICT guide that triggered the failure


def test_the_premise_holds_the_only_pridict_guide_is_really_3_prime_of_s():
    """The designer places a guide by its sequence, so this has to be true of the real window for the
    regression to be about what it says: the guide is in the window, once, and 3' of S."""
    located = locate_guide(CFG["genome"], PRIDICT_GUIDE)
    assert (
        located is not None
    ), "the guide is not in the window: the fallback below would be vacuous"
    strand, cut = located
    assert strand == "Fw" and cut > SITE_S


def _design():
    return design_recoded(
        rtt_aa=CFG["rtt_aa"],
        cargo=CFG["cargo"],
        genome=CFG["genome"],
        site_s=SITE_S,
        cds_mask=CFG["cds_mask"],
        pbs_lengths=CFG["pbs_lengths"],
        ngrna_distance_min=CFG["ngrna_distance_min"],
        ngrna_distance_max=CFG["ngrna_distance_max"],
        pegrna_strand=CFG["pegrna_strand"],
        arm_side=CFG["arm_side"],
        max_run=CFG["max_run"],
        pridict_runner=_PridictOffersOnlyADownstreamNick(),
        mode="fake",
    )


def test_unusable_pridict_nick_falls_back_to_a_usable_enumerated_one():
    assert CFG["arm_side"] == "5prime"  # the premise of the failure
    result = _design().result
    assert result.nick_position < SITE_S, "arm_side=5prime needs a nick 5' of the insertion site"
    assert result.candidates_tried >= 1


def test_the_downstream_pridict_guide_is_not_the_one_used():
    result = _design().result
    assert result.nick_position != SITE_S + 56
