"""What does the pegRNA-folding gate actually reject, and why?

    python scripts/analysis/scan_structure_gate.py CONFIG.json SPACER [THRESHOLD]


Runs the real pool walk with the REAL ViennaRNA scorer and a docking stub that
always fails (so the walk visits every candidate instead of stopping at the
first winner), and records for each candidate which region pairs the gate
counted and at what base-pair probability.
"""

import collections
import json
import sys

import pegrna_pipeline.hpc.recoded_tjpe_workflow as wf
from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.full_structure import _region_of
from pegrna_pipeline.hpc.scorer_factory import make_scorers
from pegrna_pipeline.hpc.vienna_structure import ViennaRNAStructureScorer

cfg = json.load(open(sys.argv[1]))
spacer = sys.argv[2]
threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5

fake = make_scorers(mode="fake")
vienna = ViennaRNAStructureScorer()
detail = []
_orig = wf.functional_region_pairing


def spy(seq, *, regions, scorer, threshold):
    res = _orig(seq, regions=regions, scorer=scorer, threshold=threshold)
    pairs = scorer.pair_probs(seq)
    rows = []
    for (i, j), p in pairs.items():
        if p < threshold:
            continue
        ri, rj = _region_of(i, regions), _region_of(j, regions)
        if ri is None or rj is None or "scaffold" in (ri, rj) or ri == rj:
            continue
        if {ri, rj} == {"spacer", "pbs1"}:
            continue
        rows.append((ri, rj, round(p, 2)))
    detail.append({"len": len(seq), "n": res.cross_region_pairs, "rows": rows})
    return res


wf.functional_region_pairing = spy

try:
    wf.run_recoded_tjpe(
        spacer=spacer,
        rtt_aa=cfg["rtt_aa"],
        cargo=cfg["cargo"],
        genome=cfg["genome"],
        site_s=cfg["site_s"],
        cds_mask=cfg["cds_mask"],
        pbs_lengths=cfg["pbs_lengths"],
        ngrna_distance_min=cfg["ngrna_distance_min"],
        ngrna_distance_max=cfg["ngrna_distance_max"],
        pegrna_strand=cfg["pegrna_strand"],
        arm_side=cfg["arm_side"],
        max_run=cfg["max_run"],
        splice_scorer=fake.splice_scorer,
        full_structure_scorer=vienna,
        importance_scorer=fake.importance_scorer,
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.10),  # always fails
        offtarget_scorer=fake.offtarget_scorer,
        maxent_scorer=fake.maxent_scorer,
        protein_length=cfg.get("protein_length"),
        structure_threshold=threshold,
    )
except Exception as exc:  # noqa: BLE001 - expected: every candidate is rejected
    print(f"(walk ended: {type(exc).__name__})")

n = len(detail)
bad = [d for d in detail if d["n"] > 0]
print(f"spacer={spacer}  threshold p>={threshold}  candidates reaching the structure gate: {n}")
print(f"  rejected by it : {len(bad)}   ({100 * len(bad) / max(n, 1):.0f}%)")
if n:
    print(f"  pegRNA length  : {min(d['len'] for d in detail)}-{max(d['len'] for d in detail)} nt")
pair_ct = collections.Counter()
pmax = collections.defaultdict(float)
for d in bad:
    for ri, rj, p in d["rows"]:
        k = "-".join(sorted((ri, rj)))
        pair_ct[k] += 1
        pmax[k] = max(pmax[k], p)
for k, c in pair_ct.most_common():
    print(f"  {k:<14} in {c} candidate-pair hits, max p={pmax[k]}")
if bad:
    print("  example offending rows:", bad[0]["rows"][:6])
