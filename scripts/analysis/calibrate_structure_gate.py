"""Null distribution for the pegRNA-folding gate.

Reproduces the evidence behind ``THRESHOLDS.vienna_pair_prob``. Run from the repo root:

    python scripts/analysis/calibrate_structure_gate.py

Result (300 random designs per size, seed 7): at p>=0.5 only 9-12% of tjPE-sized designs
pass, at p>=0.99 80-84% do.


Same architecture and region lengths as a real tjPE design, but with random
3'-extension sequence. Reports, per threshold, the fraction of random designs
that PASS (zero cross-region pairs). A gate that fails ~all random sequences is
not discriminating between good and bad designs.
"""

import random

from pegrna_pipeline.hpc.full_structure import functional_region_pairing
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD
from pegrna_pipeline.hpc.vienna_structure import ViennaRNAStructureScorer

rng = random.Random(7)
vienna = ViennaRNAStructureScorer()
SPACER = "TCAGCATGATGGACCGGTGA"
N = 300


def rand(n, gc=0.5):
    return "".join(rng.choice("GC" if rng.random() < gc else "AT") for _ in range(n))


def build(rcpbs2, rtt, pbs1):
    scaf = STANDARD_SGRNA_SCAFFOLD
    i1 = len(SPACER)
    i2 = i1 + len(scaf)
    i3 = i2 + len(rcpbs2)
    i4 = i3 + len(rtt)
    i5 = i4 + len(pbs1)
    seq = SPACER + scaf + rcpbs2 + rtt + pbs1
    regions = {
        "spacer": (0, i1),
        "scaffold": (i1, i2),
        "rcpbs2": (i2, i3),
        "rtt": (i3, i4),
        "pbs1": (i4, i5),
    }
    return seq, regions


for label, (r2, rt, p1) in {
    "short extension (typical PE: PBS13 + RTT~30)": (0, 30, 13),
    "tjPE-sized (rcPBS2 12 + RTT 70 + PBS 13)": (12, 70, 13),
    "tjPE-sized, long (rcPBS2 17 + RTT 98 + PBS 17)": (17, 98, 17),
}.items():
    # Per-sequence probabilities depend on the whole sequence, so cache one fold per
    # design and read it at every threshold.
    designs = []
    for _ in range(N):
        seq, regions = build(rand(r2), rand(rt), rand(p1))
        designs.append((seq, regions))
    print(f"\n{label}   ({len(designs[0][0])} nt, n={N} random designs)")
    for th in (0.5, 0.8, 0.9, 0.99):
        passed = sum(
            functional_region_pairing(s, regions=rg, scorer=vienna, threshold=th).cross_region_pairs
            == 0
            for s, rg in designs
        )
        print(f"   p>={th:<4}  random designs passing: {100 * passed / N:5.1f}%")
