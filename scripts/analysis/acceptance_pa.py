"""Acceptance check for a real-model Prime Assembly run (spec section 2, plan task 2.5).

Run ``pegrna design-pa --mode api`` (or ``real``) on the host with the models, then point this at the config and the
output directory:

    python scripts/analysis/acceptance_pa.py examples/ctcf_cterm/config_pa.json OUT_DIR

It reads the FILES a reader would order from (the FASTA, the CSVs, the JSON) and the config, and checks the design
from them, not from the pipeline's own objects:

* every output file is there and the run did not fail;
* no gate was served by a stand-in, and the run used real models;
* every gate of the walk ran, in order, and the docking gates gave a verdict (not an advisory);
* both pegRNAs are at most 166 nt, end in UUU and have no TTTT in the spacer or in the RTT followed by the PBS (the
  scaffold's own GTTTTAG is the standard one), and the flaps are the donor's ends;
* the executable strand check (``simulate_pa``) rebuilds the predicted product from the ORDERED pegRNAs and donor;
* INV-0: the product is the window with the cargo inserted at S and nothing else changed except the bases listed in
  the changes file, each of which is protein-coding and synonymous; the cargo is exactly the config's. (A change can
  only lie between the nicks: the product outside them is the window by construction.)

Exit status 0 only if every check passes.
"""

from __future__ import annotations

import csv
import io
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.Seq import Seq

from pegrna_pipeline.hpc.pa_geometry import find_guides
from pegrna_pipeline.hpc.pa_simulator import PAMismatch, simulate_pa

FILES = {
    "report": "recoded_pa_report.txt",
    "components": "recoded_pa_components.csv",
    "fasta": "recoded_pa_donor.fasta",
    "changes": "recoded_pa_changes.csv",
    "gates": "recoded_pa_gates.json",
}
GATE_ORDER = [
    "offtarget",
    "recode",
    "inv0",
    "strand",
    "cryptic",
    "splice",
    "importance",
    "structure_a",
    "structure_b",
    "docking_a",
    "docking_b",
]
MAX_PEGRNA_NT = 166

Result = Tuple[bool, str]


def _fasta(text: str) -> Dict[str, str]:
    """``{first word of the header: sequence}`` for each record."""
    out: Dict[str, str] = {}
    for record in text.split(">")[1:]:
        lines = record.strip().splitlines()
        out[lines[0].split()[0]] = "".join(lines[1:])
    return out


def run_checks(config: dict, out_dir: Path) -> List[Result]:
    results: List[Result] = []

    def check(ok: bool, text: str) -> bool:
        results.append((bool(ok), text))
        return bool(ok)

    missing = [name for name in FILES.values() if not (out_dir / name).exists()]
    if not check(
        not missing,
        "all five output files are present" + (f" (missing: {missing})" if missing else ""),
    ):
        return results
    gates = json.loads((out_dir / FILES["gates"]).read_text())
    if not check(
        not gates.get("failed"),
        "the run did not fail" + (f": {gates.get('error')}" if gates.get("failed") else ""),
    ):
        return results

    check(
        gates.get("stand_in_gates") == [],
        f"no gate was served by a stand-in (stand_in_gates = {gates.get('stand_in_gates')})",
    )
    check(
        gates.get("scorer_mode") in ("api", "real"),
        f"real models were used (scorer_mode = {gates.get('scorer_mode')})",
    )
    report = (out_dir / FILES["report"]).read_text()
    check("NOT A VALIDATED DESIGN" not in report, "the report carries no stand-in banner")
    check("extension of Prime Assembly" in report, "the report carries the standing caveat")

    info = dict(gates["gate_results"])
    check([k for k, _ in gates["gate_results"]] == GATE_ORDER, "every gate ran, in pipeline order")
    for gate in ("docking_a", "docking_b"):
        text = info.get(gate, "")
        check(
            "interface ipTM" in text and "advisory" not in text.lower(),
            f"{gate} is a verdict, not an advisory: {text[:90]}",
        )

    genome, site_s, cargo = config["genome"], config["site_s"], config["cargo"].upper()
    records = _fasta((out_dir / FILES["fasta"]).read_text())
    heads = list(records)
    check(
        len(records) == 3
        and heads[0] == "pegRNA_A"
        and heads[1] == "pegRNA_B"
        and heads[2] == "donor",
        f"the FASTA holds pegRNA A, pegRNA B and the donor (got {heads})",
    )
    peg_a, peg_b, donor = (
        records.get("pegRNA_A", ""),
        records.get("pegRNA_B", ""),
        records.get("donor", ""),
    )
    pa = gates["pa"]
    guide_a, guide_b = pa["guides"]
    for side, peg, g in (("A", peg_a, guide_a), ("B", peg_b, guide_b)):
        check(
            len(peg) <= MAX_PEGRNA_NT, f"pegRNA {side} is {len(peg)} nt (at most {MAX_PEGRNA_NT})"
        )
        check(
            peg.endswith("UUU") and peg[:20] == g["spacer"],
            f"pegRNA {side} starts with its spacer and ends in UUU",
        )
        # the parts the design chooses, as transcribed (the scaffold's own GTTTTAG is the standard one)
        chosen = (g["spacer"], g["rtt"] + g["pbs"])
        check(
            all("TTTT" not in part for part in chosen),
            f"pegRNA {side}: no TTTT (a Pol III terminator) in the spacer, or in the RTT followed by the PBS",
        )

    parts = {p["name"]: p["sequence"] for p in pa["donor"]["parts"]}
    arm5, arm3 = parts.get("arm5", ""), parts.get("arm3", "")
    check(parts.get("cargo") == cargo, "the cargo in the donor is exactly the config's cargo")
    check(donor == arm5 + cargo + arm3, "the donor is arm5 + cargo + arm3")

    n1, n2 = pa["span"]
    product = genome[:n1] + donor + genome[n2:]
    check(
        len(product) == len(genome) + len(cargo),
        f"the product is the window plus the cargo ({len(product)} nt)",
    )
    check(
        product[site_s : site_s + len(cargo)] == cargo,
        f"the cargo sits at the insertion site S ({site_s}) in the product",
    )
    by_key = {(g.strand, g.cut, g.spacer): g for g in find_guides(genome)}
    found_a = by_key.get((guide_a["strand"], guide_a["cut"], guide_a["spacer"]))
    found_b = by_key.get((guide_b["strand"], guide_b["cut"], guide_b["spacer"]))
    if check(
        found_a is not None and found_b is not None,
        "both guides are real NGG guides of the window, at the stated cuts",
    ):
        try:
            rebuilt = simulate_pa(
                genome=genome,
                a=found_a,
                b=found_b,
                donor=donor,
                peg_a=peg_a,
                peg_b=peg_b,
                flap_a=guide_a["flap_len"],
                pbs_a=guide_a["pbs_len"],
                flap_b=guide_b["flap_len"],
                pbs_b=guide_b["pbs_len"],
                tail_len=len(guide_a["tail"]),
            )
            check(
                rebuilt == product,
                "simulate_pa rebuilds the product from the ordered pegRNAs and donor",
            )
        except PAMismatch as exc:
            check(False, f"simulate_pa rejected the ordered pegRNAs and donor: {exc}")

    # INV-0: with the cargo taken out again, the product is the window except for the listed, synonymous changes.
    without_cargo = product[:site_s] + product[site_s + len(cargo) :]
    differs = {i for i, (x, y) in enumerate(zip(without_cargo, genome)) if x != y}
    rows = list(csv.DictReader(io.StringIO((out_dir / FILES["changes"]).read_text())))
    listed = {int(r["position"]) - 1 for r in rows if r["scope"] == "genomic arm"}
    check(
        len(without_cargo) == len(genome),
        "nothing is deleted: the product less the cargo is as long as the window",
    )
    check(
        differs == listed,
        f"every changed base is listed in the changes file ({len(differs)} changed, {len(listed)} listed)",
    )
    mask = list(config["cds_mask"])
    cds_start = mask.index(True) if True in mask else None
    bad = []
    for i in sorted(differs):
        if cds_start is None or not mask[i]:
            bad.append((i + 1, "not protein-coding"))
            continue
        c = i - ((i - cds_start) % 3)
        if str(Seq(genome[c : c + 3]).translate()) != str(
            Seq(without_cargo[c : c + 3]).translate()
        ):
            bad.append((i + 1, "changes the amino acid"))
    check(
        not bad,
        "every changed base is protein-coding and synonymous"
        + (f" (violations: {bad})" if bad else ""),
    )
    return results


def main(argv: List[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2
    config = json.loads(Path(argv[1]).read_text())
    results = run_checks(config, Path(argv[2]))
    for ok, text in results:
        print(("PASS  " if ok else "FAIL  ") + text)
    failed = [text for ok, text in results if not ok]
    print(f"\n{len(results) - len(failed)}/{len(results)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
