"""Where do real Prime Assembly designs fail, and which lever (off-target tolerance, flap length) rescues them?

The first real-model run of the CTCF cassette found one PAM-in pair left after the off-target screen and both of its
pegRNAs rejected by the folding gate. Before touching a lever this measures, with the real scorers (cas-offinder and
ViennaRNA), for EVERY PAM-in pair the window offers:

* each guide's off-target count, so the tolerance needed to keep it is visible;
* for each side, the number of confident cross-region pairs (what the folding gate rejects on) at every flap length
  30-50, taking the best PBS length (13, 15, 17), so a flap length that rescues a pegRNA shows up.

Run on a host with the real models (the DSI box), as the web service runs the pipeline:

    ~/pegrna/bin/dsi_run_py.sh scripts/analysis/probe_pa_levers.py examples/ctcf_cterm/config_pa.json

Nothing here changes a design: it only reads the walk's own planning and the gate's own functions.
"""

from __future__ import annotations

import json
import sys
import time

from pegrna_pipeline.hpc.full_structure import functional_region_pairing
from pegrna_pipeline.hpc.pa_geometry import assemble_pa_pegrna, derive_extension, pa_regions
from pegrna_pipeline.hpc.recoded_driver import _real_scorer_kwargs
from pegrna_pipeline.hpc.recoded_pa_workflow import plan_pairs
from pegrna_pipeline.hpc.recoded_tjpe_workflow import resolve_structure_threshold
from pegrna_pipeline.hpc.scorer_factory import make_scorers

FLAPS = range(30, 51)
PBS = (13, 15, 17)


def main(config_path: str) -> None:
    cfg = json.load(open(config_path))
    kwargs = _real_scorer_kwargs(
        uniprot_acc=cfg.get("uniprot_acc", ""),
        cas9="nCas9",
        genome_fasta=None,
        boltz_cache=None,
        rinalmo_weights=None,
        rnaformer_model_dir=None,
        pam="NGG",
        max_mismatches=3,
    )
    scorers = make_scorers(mode="api", **kwargs)
    genome, site_s, cargo = cfg["genome"], cfg["site_s"], cfg["cargo"]
    started = time.time()
    planning = plan_pairs(
        cargo=cargo,
        genome=genome,
        site_s=site_s,
        cds_mask=cfg["cds_mask"],
        offtarget_scorer=scorers.offtarget_scorer,
        pbs_lengths=PBS,
        flap_length=35,
        flap_min=30,
        flap_max=50,
        min_span=30,
        max_span=100,
        max_run=8,
        tail="UUU",
        max_offtargets=10**6,  # let every guide through: the counts are what is being measured
        min_ontarget=0.5,
    )
    print(
        f"{time.time() - started:.0f} s: {planning.pairs_geometric} PAM-in pairs, "
        f"{planning.pairs_found} after the sequence rules, {len(planning.plans)} plannable at any tolerance"
    )
    fold = scorers.full_structure_scorer
    threshold = resolve_structure_threshold(None, fold)
    print(
        f"folding gate: a pegRNA is rejected if any confident (p >= {threshold}) cross-region pair exists\n"
    )
    header = "flap " + " ".join(f"{f:2d}" for f in FLAPS)
    summary = []
    for plan in planning.plans:
        a, b = plan.pair.a, plan.pair.b
        print(
            f"pair {a.cut:3d}/{b.cut:3d}  span {plan.pair.span:3d}  off-targets A={plan.offtargets[0]} "
            f"B={plan.offtargets[1]}  advisories={len(plan.advisories)}  PAM A={plan.pam_outcomes[0].status} "
            f"B={plan.pam_outcomes[1].status}"
        )
        print(f"      {header}")
        ok_by_side = {}
        for side, guide in (("A", a), ("B", b)):
            cells = []
            passing = []
            for flap in FLAPS:
                best = None
                for pbs_len in PBS:
                    ext = derive_extension(side, genome, guide.cut, plan.donor, flap, pbs_len)
                    if ext is None:
                        continue
                    peg = assemble_pa_pegrna(spacer=guide.spacer, rtt=ext.rtt, pbs=ext.pbs)
                    pairing = functional_region_pairing(
                        peg,
                        regions=pa_regions(guide.spacer, ext.rtt, ext.pbs, "UUU"),
                        scorer=fold,
                        threshold=threshold,
                        pairs=fold.pair_probs(peg),
                    )
                    n = pairing.cross_region_pairs
                    best = n if best is None else min(best, n)
                if best is None:
                    cells.append(" .")
                else:
                    cells.append(f"{best:2d}")
                    if best == 0:
                        passing.append(flap)
            ok_by_side[side] = passing
            print(f"      {side:>4} " + " ".join(cells) + f"   passes at flap {passing or 'none'}")
        summary.append((a.cut, b.cut, plan.offtargets, ok_by_side))
        print()
    print("SUMMARY (a side 'rescuable' = some flap length 30-50 gives no cross-region pair):")
    for a_cut, b_cut, offs, ok in summary:
        both = bool(ok["A"]) and bool(ok["B"])
        print(
            f"  pair {a_cut:3d}/{b_cut:3d}  off-targets {offs}  A passes at {len(ok['A']):2d} flaps, "
            f"B at {len(ok['B']):2d}  -> {'folds at some flap length' if both else 'never folds clean'}"
        )


if __name__ == "__main__":
    main(sys.argv[1])
