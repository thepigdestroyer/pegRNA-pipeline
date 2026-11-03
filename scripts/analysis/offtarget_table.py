"""Per-candidate off-target table for the nick (ngRNA) screen.

    python scripts/analysis/offtarget_table.py CONFIG.json /path/to/hg38.2bit

Needs the cas-offinder binary on PATH. One whole-genome scan (~45-100 s).


Runs the pipeline's own CasOffinderScorer plumbing once for every candidate and
prints what cas-offinder actually found, split by mismatch count and by contig
type, so we can see WHY the screen empties the pool instead of guessing.
"""

import collections
import json
import os
import sys
import tempfile
import time

from pegrna_pipeline.hpc.nick_search import enumerate_nick_sites
from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer

cfg = json.load(open(sys.argv[1]))
genome2bit = sys.argv[2]

cands = enumerate_nick_sites(
    cfg["genome"],
    site_s=cfg["site_s"],
    pegrna_strand=cfg.get("pegrna_strand", "Rv"),
    distance_min=cfg.get("ngrna_distance_min", 40),
    distance_max=cfg.get("ngrna_distance_max", 90),
)
print(f"window={len(cfg['genome'])} bp  site_s={cfg['site_s']}  strand={cfg.get('pegrna_strand')}")
print(
    f"{len(cands)} candidate nick guides in the {cfg.get('ngrna_distance_min')}-"
    f"{cfg.get('ngrna_distance_max')} bp window\n",
    flush=True,
)

sc = CasOffinderScorer(genome=genome2bit, pam="NGG", max_mismatches=3)
protos = sorted({c.protospacer for c in cands})
rows = collections.defaultdict(list)
t0 = time.time()
with tempfile.TemporaryDirectory() as tmp:
    inp, outp = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
    sc._write_input(inp, protos, spacer_len=20)
    sc._run(inp, outp)
    for line in open(outp):
        p = line.split()
        if len(p) < 6:
            continue
        guide = p[0][: len(p[0]) - 3]
        rows[guide].append((p[1], int(p[2]), int(p[5])))
print(f"cas-offinder: {time.time() - t0:.0f}s\n", flush=True)


def primary(chrom):
    return not any(t in chrom for t in ("_alt", "_random", "chrUn", "_fix"))


print(
    f"{'nick':>5} {'dist':>4} {'strand':>6}  {'protospacer':<20}  "
    f"{'mm0':>3} {'mm1':>3} {'mm2':>3} {'mm3':>3}  {'mm0(primary)':>12}  verdict(prod rule)"
)
survivors = 0
for c in cands:
    hits = rows.get(c.protospacer, [])
    n = collections.Counter(mm for _, _, mm in hits)
    prim = collections.Counter(mm for ch, _, mm in hits if primary(ch))
    exact, near = n[0], n[1] + n[2]
    off = max(0, exact - 1) + near
    ok = exact >= 1 and off == 0
    survivors += ok
    why = (
        "KEEP"
        if ok
        else (
            "drop: no exact match in genome"
            if exact == 0
            else f"drop: {max(0, exact - 1)} duplicate exact + {near} hits at <=2 mismatches"
        )
    )
    print(
        f"{c.nick_position:>5} {c.distance_to_edit:>4} {c.target_strand:>6}  {c.protospacer:<20}  "
        f"{n[0]:>3} {n[1]:>3} {n[2]:>3} {n[3]:>3}  {prim[0]:>12}  {why}"
    )

print(f"\nproduction rule keeps {survivors}/{len(cands)}")
print("\nlocations of every hit at <=2 mismatches (the ones that disqualify):")
for c in cands:
    for ch, pos, mm in sorted(rows.get(c.protospacer, []), key=lambda r: r[2]):
        if mm <= 2:
            print(f"  nick {c.nick_position:>4} {c.protospacer}  {ch}:{pos}  mm={mm}")
