"""What would a Prime Assembly designer find at the CTCF C-terminal locus?  (read-only probe)

    python scripts/analysis/probe_pa_geometry.py

Evidence for docs/superpowers/specs/2026-11-04-prime-assembly-design.md. It builds no pegRNA and calls no
model. It uses the repo's locus config and its own recoder to answer, on real sequence:

  1. Do PAM-in guide pairs exist with the insertion site S between the two nicks?
  2. How much genomic DNA would the donor have to restore (arm5 + arm3), and can it be recoded?
  3. How much of each 35-nt flap is genomic? That stretch is a copy of the strand the flap displaces,
     so it competes with the donor for annealing; the recoder exists to break exactly this.
  4. Where do the two PAMs fall? Coding PAMs can be changed silently; non-coding ones are byte-exact.

Coordinates. Every cut is "b": the index of the FIRST base to the right of the cut on the sense window,
genome[:b] | genome[b:], for a guide on either strand. (``nick_search.enumerate_nick_sites`` reports a
different convention for reverse-strand guides: the last base LEFT of the cut. That is a separate bug,
see the spec; it is not used here.)

PAM-in is not a choice but a consequence of SpCas9: the nick is 3 bp 5' of the PAM and the flap extends 3'
from the nick, toward the PAM. For two flaps to fill the gap between the nicks, the top-strand guide must
lie on the left and the bottom-strand guide on the right, so both PAMs sit inside the replaced span.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from Bio.Seq import Seq  # noqa: E402

from pegrna_pipeline.agent.motif_library import MotifLibrary  # noqa: E402
from pegrna_pipeline.hpc.guide_sanity import assess_spacer  # noqa: E402
from pegrna_pipeline.hpc.microhomology import longest_common_run  # noqa: E402
from pegrna_pipeline.hpc.recoded_tjpe_workflow import _partition_arm  # noqa: E402

FLAP = 35  # nt; the paper: 30 and 50 robust, 13 very weak
PBS_MAX = 13  # room needed beside a cut for a 13-nt PBS
SPAN_MIN, SPAN_MAX = (
    30,
    100,
)  # nick-to-nick, bp (project floor 23: pair_scorer.MIN_NICK_SEPARATION_BP)
PAM_CLASH = 12  # the two PAM-bearing footprints overlap below this span (geometric minimum)


def rc(s):
    return str(Seq(s).reverse_complement())


cfg = json.loads((REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json").read_text())
G = cfg["genome"].upper()
S = cfg["site_s"]
MASK = cfg["cds_mask"]
MAX_RUN = cfg.get("max_run", 8)


def coding(p):
    return bool(MASK[p]) if p < len(MASK) else False


def find_guides():
    """Every NGG-bearing guide in the window, on both strands, with its cut as defined above."""
    top, bottom = [], []
    for i in range(20, len(G) - 2):  # top strand: protospacer G[i-20:i], PAM G[i:i+3] = NGG
        if G[i + 1 : i + 3] == "GG":
            top.append(dict(cut=i - 3, spacer=G[i - 20 : i], pam=G[i : i + 3], pam_pos=(i, i + 3)))
    for j in range(0, len(G) - 22):  # bottom strand shows CCN at j; protospacer = RC(G[j+3:j+23])
        if G[j : j + 2] == "CC":
            bottom.append(
                dict(
                    cut=j + 6,
                    spacer=rc(G[j + 3 : j + 23]),
                    pam=rc(G[j : j + 3]),
                    pam_pos=(j, j + 3),
                )
            )
    return top, bottom


def pam_class(g):
    flags = [coding(p) for p in range(*g["pam_pos"])]
    return "coding" if all(flags) else "non-coding" if not any(flags) else "mixed"


def arm5(g):  # genome[n1:S): the 5' stretch the donor must restore
    return (
        _partition_arm(G, lo=g["cut"], hi=S, cds_mask=MASK, max_run=MAX_RUN)
        if g["cut"] < S
        else None
    )


def arm3(g):  # genome[S:n2): the 3' stretch
    return (
        _partition_arm(G, lo=S, hi=g["cut"], cds_mask=MASK, max_run=MAX_RUN)
        if g["cut"] > S
        else None
    )


def main():
    top, bottom = find_guides()
    five = [g for g in top if PBS_MAX + 1 <= g["cut"] <= S]
    three = [g for g in bottom if S <= g["cut"] <= len(G) - PBS_MAX]
    pairs = [
        (a, b, b["cut"] - a["cut"])
        for a in five
        for b in three
        if SPAN_MIN <= b["cut"] - a["cut"] <= SPAN_MAX and b["cut"] - a["cut"] >= PAM_CLASH
    ]
    print(f"window {len(G)} nt; S = {S}, stop codon {G[S:S + 3]} starts at S; flap {FLAP} nt")
    print(f"NGG-bearing guides: {len(top)} top strand, {len(bottom)} bottom strand")
    print(
        f"usable 5' guides (top strand, nick <= S): {len(five)}   3' guides (bottom, nick >= S): {len(three)}"
    )
    print(f"PAM-in pairs with S between the nicks, span {SPAN_MIN}-{SPAN_MAX} bp: {len(pairs)}")

    print("\n5' guides (top strand), nearest S first")
    print("   cut  arm5   PAM         recoding the restored coding arm      sane")
    for a in sorted(five, key=lambda g: -g["cut"]):
        built = arm5(a)
        what = (
            "no arm (nick at S)"
            if a["cut"] == S
            else (
                f"{built.n_edits:2d} silent edits, longest identical run {built.max_run_after}"
                if built
                else "NOT RECODABLE"
            )
        )
        print(
            f"   {a['cut']:3d}  {S - a['cut']:3d} nt  {a['pam']} {pam_class(a):10s}  {what:42s} "
            f"{assess_spacer(a['spacer']).ok}"
        )

    print("\n3' guides (bottom strand), nearest S first")
    print("   cut  arm3   PAM         restored 3' arm                       sane")
    for b in sorted(three, key=lambda g: g["cut"]):
        built = arm3(b)
        what = (
            "no arm (nick at S)"
            if b["cut"] == S
            else (
                "non-coding: byte-exact, cannot be recoded"
                if built and built.locked_span and not built.recode_span
                else "partly coding" if built else "NOT BUILDABLE"
            )
        )
        print(
            f"   {b['cut']:3d}  {b['cut'] - S:3d} nt  {b['pam']} {pam_class(b):10s}  {what:42s} "
            f"{assess_spacer(b['spacer']).ok}"
        )

    lib = MotifLibrary.load_full()
    cargos = {
        "FLAG, 24 bp": lib.lookup("FLAG").dna,
        "GGGGS-FLAG-P2A-ZeoR, 503 bp": "".join(
            lib.lookup(n).dna for n in ("GGGGS", "FLAG", "P2A", "ZeocinR")
        ),
    }
    print("\ndonors for the pairs with the least restored DNA (arm5 + arm3)")
    ranked = []
    for a, b, span in pairs:
        a5, a3 = arm5(a), arm3(b)
        if (a["cut"] < S and a5 is None) or (b["cut"] > S and a3 is None):
            continue
        ranked.append((S - a["cut"] + b["cut"] - S, span, a, b, a5, a3))
    for _, span, a, b, a5, a3 in sorted(ranked, key=lambda t: t[:2])[:6]:
        left = a5.sense_arm if a5 else ""
        right = a3.sense_arm if a3 else ""
        for name, cargo in cargos.items():
            donor = left + cargo + right
            run_a = longest_common_run(donor[:FLAP], G[a["cut"] : a["cut"] + FLAP])
            run_b = longest_common_run(donor[-FLAP:], G[b["cut"] - FLAP : b["cut"]])
            note = (
                "  <- the two flaps overlap: this is twinPE, not assembly"
                if len(donor) < 2 * FLAP
                else ""
            )
            print(
                f"  A@{a['cut']:3d} ({pam_class(a)[:3]} PAM) + B@{b['cut']:3d} ({pam_class(b)[:3]} PAM)  "
                f"span {span:3d} = {len(left):2d}+{len(right):2d} | {name:27s} donor {len(donor):4d} nt | "
                f"genomic run under flap A {run_a:2d}, B {run_b:2d}{note}"
            )


if __name__ == "__main__":
    main()
