"""Can the legacy Prime Assembly path pair REAL PRIDICT output?  (It cannot: issue #29, ADR 0027.)

Needs a working PRIDICT2 install, so run it where the pipeline's real models are (the DSI box):

    set -a; . ~/.config/pegrna/pegrna.env; set +a
    python scripts/analysis/probe_legacy_pa.py

It builds the two flap windows the way the legacy HPC driver does (each flap inserted at the same
site, here the CTCF stop codon), runs real PRIDICT2 once for both, and applies the legacy
``filter_compatible_pairs``.

What it shows, and why the hand-written fixtures hid it. The legacy filter treats
``editing_position`` as a coordinate on a shared axis and requires ``|up - down|`` to be 30-80.
Real PRIDICT reports the edit's distance FROM THE NICK (here 0-17), the same set for both runs,
so the distance is at most 17 and no pair can ever pass. The fixtures used 50/58 against
100/110, which no real run produces.

The problem is structural, not a threshold: PRIDICT only proposes pegRNAs whose nick lies within
~17 bp of the edit, so for one insertion site the two nicks are at most ~28 bp apart, while a
Prime Assembly pair needs nicks that flank the region to be replaced (the published design,
Liu et al., Nature 2026, uses ~35 nt flaps from two pegRNAs on opposite strands). A Prime
Assembly designer has to choose the nicks itself.
"""

import json
import os
from pathlib import Path

from pegrna_pipeline.hpc.input_formatter import format_tjpe_insertion_editseq
from pegrna_pipeline.hpc.pair_scorer import filter_compatible_pairs
from pegrna_pipeline.hpc.scorer_factory import make_scorers

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json"
MIN_FLANK = 100

# 30 nt each: the first and last flap of a longer construct (FLAG-GS ... a P2A-like tail).
FLAP_UP = "GACTACAAGGACGACGATGACAAGGGAAGC"
FLAP_DN = "GCTACTAACTTCAGCCTGCTGAAGCAGGCT"


def main() -> None:
    cfg = json.loads(CONFIG.read_text())
    genome, site_s = cfg["genome"], cfg["site_s"]
    runner = make_scorers(
        mode="api",
        uniprot_acc=cfg.get("uniprot_acc", "P49711"),
        genome=os.environ.get("PEGRNA_GENOME"),
        nvidia_api_key=None,
    ).pridict_runner
    if runner is None:
        raise SystemExit("no PRIDICT runner: set PRIDICT2_DIR / PRIDICT2_PYTHON (see GO-LIVE.md)")

    rows = [
        {
            "sequence_name": name,
            "editseq": format_tjpe_insertion_editseq(genome, site_s, flap, min_flank=MIN_FLANK),
        }
        for name, flap in (("pa_up", FLAP_UP), ("pa_dn", FLAP_DN))
    ]
    batch = runner.run_batch(rows)
    designs = {n: (batch.get(n) or {}).get("pegrnas") or [] for n in ("pa_up", "pa_dn")}

    for name, ds in designs.items():
        print(f"{name}: {len(ds)} pegRNAs")
    observed = {
        name: sorted({(d.target_strand, d.editing_position) for d in ds})
        for name, ds in designs.items()
    }
    print("distinct (strand, editing_position) per run:", observed)

    pairs = filter_compatible_pairs(
        upstream=designs["pa_up"], downstream=designs["pa_dn"], stop_codon_position=site_s
    )
    print(f"compatible pairs under the legacy filter: {len(pairs)}")
    largest = max(
        (p for ds in designs.values() for d in ds for p in (d.editing_position,)), default=0
    )
    print(f"largest editing_position seen: {largest}  (the filter needs a difference of >= 30)")


if __name__ == "__main__":
    main()
