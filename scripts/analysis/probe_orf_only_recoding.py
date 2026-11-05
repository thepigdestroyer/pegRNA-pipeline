"""Would recoding only the ORF of a cargo (its tail kept verbatim) change what the folding gate sees?

The CTCF cassette (GGGGS-FLAG-GSG-P2A-resistance gene, then a T7 promoter) is rejected by the folding gate on one
4-bp stem between the first four bases of the rcPBS2 and the first four of the RTT, which are the reverse
complement of the cassette's last four bases. ADR 0022's INSERT lever swaps the cargo for synonymous variants, but
refuses a cargo that has any stop codon in its translation, and this one has three: the real terminator after the
resistance gene, and two more from reading the T7 promoter that follows it in the same frame. Those two are never
translated. This asks what the lever would have done had its guard let it recode the ORF alone.

Run where ViennaRNA is installed (the DSI box, or any env with the pipeline):

    python scripts/analysis/probe_orf_only_recoding.py /path/to/config.json [--spacer SEQ] [--nick 101] [--n 30]

``--spacer`` is the pegRNA#1 spacer PRIDICT chose for the window (the config leaves it to PRIDICT), and the default
is the one it chooses for the CTCF example. Output: how many random synonymous recodings of the ORF clear the gate,
and what changing only the tail does.
"""

from __future__ import annotations

import argparse
import json
import random

from pegrna_pipeline.hpc.full_structure import functional_region_pairing
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    _assemble_edit,
    _derive_pbs1,
    _partition_arm,
    _regions_for,
)
from pegrna_pipeline.hpc.strand_orientation import orient_rtt_tail
from pegrna_pipeline.hpc.tjpe_workflow import (
    STANDARD_SGRNA_SCAFFOLD,
    assemble_tj_pegrna_v2,
    derive_rc_pbs2,
)
from pegrna_pipeline.hpc.vienna_structure import ViennaRNAStructureScorer

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # standard code, TCAG order
CODE = {
    a + b + d: _AAS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, d in enumerate(_BASES)
}
SYNONYMS: dict = {}
for _codon, _aa in CODE.items():
    SYNONYMS.setdefault(_aa, []).append(_codon)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--spacer", default="TCAGCATGATGGACCGGTGA")
    ap.add_argument("--nick", type=int, default=101)
    ap.add_argument("--rcpbs2-len", type=int, default=12)
    ap.add_argument("--pbs1-len", type=int, default=13)
    ap.add_argument("--n", type=int, default=30)
    args = ap.parse_args()

    cfg = json.load(open(args.config))
    genome, site_s, cds_mask = cfg["genome"].upper(), cfg["site_s"], cfg["cds_mask"]
    cargo = cfg["cargo"].upper()

    scorer = ViennaRNAStructureScorer()
    arm = _partition_arm(genome, lo=args.nick, hi=site_s, cds_mask=cds_mask, max_run=8)
    rcpbs2 = derive_rc_pbs2(genome, ngrna_nick_position=args.nick, pbs2_length=args.rcpbs2_len)
    pbs1 = _derive_pbs1(genome, site_s=site_s, pbs_length=args.pbs1_len, pegrna_strand="Rv")

    def pairs_for(c: str) -> int:
        flap, _edited, _lo, _hi = _assemble_edit(
            genome=genome, site_s=site_s, nick_position=args.nick, arm_sense=arm.sense_arm, cargo=c
        )
        rtt = orient_rtt_tail(flap, pegrna_strand="Rv")
        peg = assemble_tj_pegrna_v2(
            spacer=args.spacer,
            pbs2=rcpbs2,
            rtt=rtt,
            pbs1=pbs1,
            polyu="UUU",
            scaffold=STANDARD_SGRNA_SCAFFOLD,
        )
        regions = _regions_for(args.spacer, rcpbs2, rtt, pbs1)
        return functional_region_pairing(
            peg, regions=regions, scorer=scorer, threshold=0.99
        ).cross_region_pairs

    first_stop = next(i for i in range(0, len(cargo) - 2, 3) if CODE[cargo[i : i + 3]] == "*")
    orf, tail = cargo[:first_stop], cargo[first_stop:]
    print(f"cargo {len(cargo)} nt = ORF {len(orf)} nt ({len(orf) // 3} codons) + tail {len(tail)} nt: {tail}")
    print(f"as given: {pairs_for(cargo)} cross-region pairs\n")

    rng = random.Random(20261105)
    recoded = [
        "".join(rng.choice(SYNONYMS[CODE[orf[i : i + 3]]]) for i in range(0, len(orf), 3))
        for _ in range(args.n)
    ]
    counts = [pairs_for(r + tail) for r in recoded]
    print(f"{args.n} random synonymous recodings of the ORF only, tail kept verbatim:")
    print(
        f"   with no cross-region pair: {sum(1 for c in counts if c == 0)}/{args.n}"
        f"    pair counts: {sorted(set(counts))}\n"
    )

    print("the tail's own bases, ORF untouched:")
    for label, t in (
        ("as given              ", tail),
        ("last codon GGG -> GGC ", tail[:-3] + "GGC"),
        ("last codon GGG -> GGA ", tail[:-3] + "GGA"),
        ("drop the last 2 bases ", tail[:-2]),
    ):
        print(f"   {label} ...{t[-9:]}  ->  {pairs_for(orf + t)} pairs")


if __name__ == "__main__":
    main()
