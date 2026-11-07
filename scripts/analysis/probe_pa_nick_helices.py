"""Which cross-region pairs does the folding gate count on Prime Assembly pegRNAs, and are they the ones a donor that
copies the genome at the nick makes inevitable?  (real ViennaRNA; run on the host with the models)

    ~/pegrna/bin/dsi_run_py.sh scripts/analysis/probe_pa_nick_helices.py examples/ctcf_cterm/config_pa.json

Evidence for ADR 0030. The first real-model Prime Assembly run of the CTCF cassette found one guide pair left after
the off-target screen, and the folding gate rejected both of its pegRNAs on a single base pair (guide A) and two
(guide B), at every flap length and every PBS length. ``probe_pa_levers.py`` showed it was not that pair: no pair of
the window ever folded clean.

The prediction tested here. The donor copies the genome at each nick, so the RTT's last three bases are the reverse
complement of the three bases just right of the nick. SpCas9 cuts 3 bp from its PAM and the protospacer is continuous
across the nick, so those three bases are the spacer's own last three. The spacer-PBS helix therefore runs on into the
RTT: spacer base ``17 + m`` pairs with the RTT's last-but-``m`` base, for ``m`` = 0, 1, 2. A donor that restores the
DNA between the nicks has this at every pegRNA; tj-PE starts its insert at the nick and does not.

For every PAM-in pair of the window (no off-target filter, so all of them), both pegRNAs, at flap lengths 30, 35, 42
and 50 (PBS 15), this folds the pegRNA and recounts the pairs the gate counts (the same rule as
``functional_region_pairing``, written out again here so that the count does not depend on the exemption under
test), then splits them into the three predicted position pairs and everything else. If the prediction is right,
"other" is zero.

Measured on the DSI box, the CTCF cassette window, ViennaRNA at the gate's threshold (p >= 0.99): 96 pegRNAs, 230
counted pairs, all 230 of them the predicted ones, none in any other place; 0 pegRNAs with anything else.
"""

import json
import sys

from pegrna_pipeline.hpc.full_structure import _region_of
from pegrna_pipeline.hpc.pa_geometry import (
    SPACER_BASES_RIGHT_OF_NICK,
    assemble_pa_pegrna,
    derive_extension,
    pa_regions,
)
from pegrna_pipeline.hpc.recoded_driver import _real_scorer_kwargs
from pegrna_pipeline.hpc.recoded_pa_workflow import plan_pairs
from pegrna_pipeline.hpc.recoded_tjpe_workflow import resolve_structure_threshold
from pegrna_pipeline.hpc.scorer_factory import make_scorers

FLAPS = (30, 35, 42, 50)
PBS = 15
# The flap length whose pegRNAs are listed one by one, for the first SHOW_PAIRS pairs.
SHOW_FLAP = 35
SHOW_PAIRS = 2


def main(config_path: str) -> int:
    cfg = json.load(open(config_path))
    kw = _real_scorer_kwargs(
        uniprot_acc=cfg["uniprot_acc"],
        cas9="nCas9",
        genome_fasta=None,
        boltz_cache=None,
        rinalmo_weights=None,
        rnaformer_model_dir=None,
        pam="NGG",
        max_mismatches=3,
    )
    scorers = make_scorers(mode="api", **kw)
    planning = plan_pairs(
        cargo=cfg["cargo"],
        genome=cfg["genome"],
        site_s=cfg["site_s"],
        cds_mask=cfg["cds_mask"],
        offtarget_scorer=scorers.offtarget_scorer,
        pbs_lengths=(13, 15, 17),
        flap_length=35,
        flap_min=30,
        flap_max=50,
        min_span=30,
        max_span=100,
        max_run=8,
        tail="UUU",
        max_offtargets=10**6,
        min_ontarget=0.5,
    )
    fold = scorers.full_structure_scorer
    threshold = resolve_structure_threshold(None, fold)
    total = {"pegrnas": 0, "counted": 0, "continuity": 0, "other": 0, "with_other": 0}

    for n, plan in enumerate(planning.plans):
        for side, guide in (("A", plan.pair.a), ("B", plan.pair.b)):
            for flap in FLAPS:
                ext = derive_extension(side, cfg["genome"], guide.cut, plan.donor, flap, PBS)
                if ext is None:
                    continue
                peg = assemble_pa_pegrna(spacer=guide.spacer, rtt=ext.rtt, pbs=ext.pbs)
                regions = pa_regions(guide.spacer, ext.rtt, ext.pbs, "UUU")
                rtt_end = regions["rtt"][1]
                predicted = {
                    (regions["spacer"][1] - SPACER_BASES_RIGHT_OF_NICK + m, rtt_end - 1 - m)
                    for m in range(SPACER_BASES_RIGHT_OF_NICK)
                }
                counted, continuity, other = [], [], []
                for (i, j), p in fold.pair_probs(peg).items():
                    if p < threshold:
                        continue
                    ri, rj = _region_of(i, regions), _region_of(j, regions)
                    if (
                        ri is None
                        or rj is None
                        or "scaffold" in (ri, rj)
                        or ri == rj
                        or {ri, rj} == {"spacer", "pbs1"}
                    ):
                        continue  # not a pair the gate counts
                    row = (i, j, ri, rj, round(p, 3))
                    counted.append(row)
                    (continuity if (i, j) in predicted else other).append(row)
                total["pegrnas"] += 1
                total["counted"] += len(counted)
                total["continuity"] += len(continuity)
                total["other"] += len(other)
                total["with_other"] += 1 if other else 0
                if flap == SHOW_FLAP and n < SHOW_PAIRS:
                    label = f"pair {plan.pair.a.cut}/{plan.pair.b.cut} {side} F={flap}"
                    print(
                        f"{label}: counted {len(counted)} = continuity {len(continuity)} + other {len(other)}"
                    )
                    print("    continuity:", continuity)
                    print("    other     :", other)

    pairs = len(planning.plans)
    print(
        f"\nTOTAL over {total['pegrnas']} pegRNAs ({pairs} pairs x 2 sides x flaps {FLAPS}, PBS {PBS}):",
        total,
    )
    predicted_only = total["counted"] > 0 and total["other"] == 0
    print(
        "every counted pair is one of the predicted nick-continuity pairs"
        if predicted_only
        else "the prediction does NOT hold: see 'other' above"
    )
    return 0 if predicted_only else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__.split("\n\n")[0] + "\n\nusage: probe_pa_nick_helices.py CONFIG.json")
    sys.exit(main(sys.argv[1]))
