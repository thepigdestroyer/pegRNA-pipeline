"""Does PRIDICT's nicking-guide position follow the same convention for both strands?  (issue #31)

Needs a working PRIDICT2 install, so run it where the pipeline's real models are (the DSI box):

    set -a; . ~/.config/pegrna/pegrna.env; set +a
    python scripts/analysis/probe_pridict_nick_convention.py

The PRIDICT path of the designer places the second nick at ``site_s + Nicking-Position-to-edit``. Whether
that equals the true cut for a guide on the reverse strand as well as the forward one is a fact about
PRIDICT, not something the code can assume. This runs real PRIDICT2 on the CTCF example window and, for
every nicking guide it returns, finds the protospacer in the window, computes the true cut independently
(SpCas9 nicks the PAM strand 3 bp 5' of the PAM), and prints the difference. A cut is the index of the
first base RIGHT of it (``genome[:b] | genome[b:]``).
"""

import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from Bio.Seq import Seq  # noqa: E402

from pegrna_pipeline.hpc.input_formatter import format_tjpe_insertion_editseq  # noqa: E402
from pegrna_pipeline.hpc.scorer_factory import make_scorers  # noqa: E402

MIN_FLANK = 100


def rc(s):
    return str(Seq(s).reverse_complement())


def main() -> None:
    cfg = json.loads((REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json").read_text())
    genome, site_s = cfg["genome"].upper(), cfg["site_s"]
    runner = make_scorers(
        mode="api",
        uniprot_acc=cfg.get("uniprot_acc", "P49711"),
        genome=os.environ.get("PEGRNA_GENOME"),
        nvidia_api_key=None,
    ).pridict_runner
    if runner is None:
        raise SystemExit("no PRIDICT runner: set PRIDICT2_DIR / PRIDICT2_PYTHON (see GO-LIVE.md)")

    editseq = format_tjpe_insertion_editseq(genome, site_s, cfg["cargo"][:40], min_flank=MIN_FLANK)
    entry = runner.run_batch([{"sequence_name": "peg1", "editseq": editseq}]).get("peg1") or {}
    rows = entry.get("nicking") or []
    print(f"window {len(genome)} nt, S = {site_s}; PRIDICT returned {len(rows)} nicking guides\n")
    print("  label  found-on   PRIDICT distance   true cut - S   difference (PRIDICT - true)")

    by_side = {}
    for row in rows:
        proto = row["Nicking-Protospacer"].upper()
        dist = int(row["Nicking-Position-to-edit"])
        i, j = genome.find(proto), genome.find(rc(proto))
        if i >= 0:
            side, cut = "top", i + 17
        elif j >= 0:
            side, cut = "bottom", j + 3
        else:
            print(f"  {row.get('Target_Strand', '?'):5s}  (protospacer not found in the window)")
            continue
        diff = dist - (cut - site_s)
        by_side.setdefault((row.get("Target_Strand", "?"), side), set()).add(diff)
        print(f"  {row.get('Target_Strand', '?'):5s}  {side:8s}  {dist:+9d}        {cut - site_s:+9d}      {diff:+d}")

    print("\nsummary (PRIDICT label, strand found on) -> set of differences:")
    for key, diffs in sorted(by_side.items()):
        print(f"  {key}: {sorted(diffs)}")
    agree = all(d == {0} for d in by_side.values())
    print("\nPRIDICT's position IS the true cut, on every strand seen." if agree else
          "\nPRIDICT's position differs from the true cut for at least one strand: normalise it.")


if __name__ == "__main__":
    main()
