"""What does the folding gate see in a pegRNA that carries a 501-bp cassette, and is it evidence?

After ADR 0028 the CTCF cassette run gets past the splice gate and is then rejected on every candidate by the
folding gate ("cross-region pairing"). That gate counts base pairs of probability >= 0.99 (ViennaRNA's partition
function, ADR 0024) between two DIFFERENT regions of the pegRNA: spacer, rcPBS2, RTT, PBS1. ADR 0024 calibrated it
on pegRNAs of 139-228 nt, where it already rejects 12-20 % of random extensions. This pegRNA is ~690 nt: the RTT
alone is 547. A confident pair somewhere between two regions gets more likely the longer the longest region is.

Run where the pipeline's real models are (the DSI box):

    dsi_run_py.sh probe_structure_gate.py /path/to/run/config.json [--n 10]

Parts: (1) CAPTURE the folds the real pipeline path makes for this config, (2) DESCRIBE the confident cross-region
helices it found, (3) SCALE: random inserts of growing length in the same pegRNA backbone, how often does the gate
reject, (4) NULL: random synonymous recodings of this cassette, and the shipped 24-bp FLAG as a control.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import random
import sys
import time
from pathlib import Path

THRESHOLD_FALLBACK = 0.99


def rc(s: str) -> str:
    from Bio.Seq import Seq

    return str(Seq(s).reverse_complement())


def region_name(i: int, regions: dict) -> str:
    for name, (lo, hi) in regions.items():
        if lo <= i < hi:
            return name
    return "tail"


def helices(pairs: dict, regions: dict, threshold: float):
    """Confident cross-region pairs grouped into stacked helices: (region_a, a0, a1, region_b, b0, b1, bp, min_p)."""
    keep = sorted(
        (i, j, p)
        for (i, j), p in pairs.items()
        if p >= threshold
        and region_name(i, regions) not in ("scaffold", "tail")
        and region_name(j, regions) not in ("scaffold", "tail")
        and region_name(i, regions) != region_name(j, regions)
        and {region_name(i, regions), region_name(j, regions)} != {"spacer", "pbs1"}
    )
    out, used = [], set()
    index = {(i, j): p for i, j, p in keep}
    for i, j, p in keep:
        if (i, j) in used:
            continue
        run = [(i, j)]
        used.add((i, j))
        while (run[-1][0] + 1, run[-1][1] - 1) in index and (run[-1][0] + 1, run[-1][1] - 1) not in used:
            nxt = (run[-1][0] + 1, run[-1][1] - 1)
            run.append(nxt)
            used.add(nxt)
        ps = [index[x] for x in run]
        a0, b0 = run[0]
        a1, b1 = run[-1]
        out.append((region_name(a0, regions), a0, a1, region_name(b0, regions), b1, b0, len(run), min(ps)))
    return out


def capture(config_path: str, out_dir: Path):
    from pegrna_pipeline.agent import cli as cli_mod
    from pegrna_pipeline.hpc import recoded_tjpe_workflow as wf

    folds = []
    real = wf.functional_region_pairing

    def recorder(seq, *, regions, scorer, threshold, pairs=None):
        got = pairs if pairs is not None else scorer.pair_probs(seq)
        folds.append(
            {
                "seq": seq,
                "regions": dict(regions),
                "threshold": threshold,
                "pairs": {k: v for k, v in got.items() if v >= 0.5},
                "scorer": scorer,
            }
        )
        return real(seq, regions=regions, scorer=scorer, threshold=threshold, pairs=pairs)

    wf.functional_region_pairing = recorder
    buf, code = io.StringIO(), 0
    try:
        with contextlib.redirect_stdout(buf):
            cli_mod.design_recoded_cmd.main(
                ["--config", config_path, "--mode", "api", "--output-dir", str(out_dir), "--progress"],
                standalone_mode=False,
            )
    except SystemExit as exc:
        code = int(exc.code or 0)
    finally:
        wf.functional_region_pairing = real
    progress = []
    for line in buf.getvalue().splitlines():
        if line.startswith("@@P "):
            with contextlib.suppress(ValueError):
                progress.append(json.loads(line[4:]))
    return folds, progress, code


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--out", default=str(Path.home() / "pegrna" / "dev-runs" / "structure_probe"))
    ap.add_argument(
        "--flag-config",
        default=str(Path.home() / "pegrna" / "prime-editing" / "examples" / "ctcf_cterm" / "config_autospacer_flag.json"),
    )
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    cargo = cfg["cargo"].upper()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print(f"config: cargo {len(cargo)} nt, pbs_lengths {cfg.get('pbs_lengths')}, arm_side {cfg.get('arm_side')}")
    print("\n[1] CAPTURE: the folds the real pipeline path makes")
    folds, progress, code = capture(args.config, out_dir)
    print(f"    {len(folds)} fold(s) reached the gate; exit code {code}")
    for p in progress:
        if p.get("gate") == "structure":
            print(f"    cand{p['cand']}: passed={p['passed']}  {p['info']}")
    if not folds:
        print("    the folding gate was never reached")
        return
    from pegrna_pipeline.hpc.full_structure import functional_region_pairing
    from pegrna_pipeline.hpc.recoded_tjpe_workflow import _regions_for

    first = folds[0]
    seq, regions, thr, scorer = first["seq"], first["regions"], first["threshold"], first["scorer"]
    print(f"    pegRNA {len(seq)} nt; regions { {k: (v[1] - v[0]) for k, v in regions.items()} }; cutoff p >= {thr}")

    print("\n[2] DESCRIBE: the confident cross-region helices (p >= cutoff), per candidate")
    for n, f in enumerate(folds, 1):
        hs = helices(f["pairs"], f["regions"], f["threshold"])
        print(f"    cand{n} ({len(f['seq'])} nt): {len(hs)} helix/helices, {sum(h[6] for h in hs)} bp")
        for ra, a0, a1, rb, b1, b0, bp, mp in hs[:6]:
            sa, sb = f["seq"][a0 : a1 + 1], f["seq"][b1 : b0 + 1]
            print(f"        {ra}[{a0 - f['regions'][ra][0]}..{a1 - f['regions'][ra][0]}] x "
                  f"{rb}[{b1 - f['regions'][rb][0]}..{b0 - f['regions'][rb][0]}]  {bp} bp, min p {mp:.3f}   {sa} / {sb}")

    # backbone: spacer, scaffold, rcpbs2, [rtt = RC(cargo) + RC(arm)], pbs1, tail
    a, b = regions["rtt"]
    rtt = seq[a:b]
    if rtt[: len(cargo)].replace("U", "T") != rc(cargo):
        print("    NOTE: the RTT does not start with RC(cargo); the scaling below keeps the whole RTT as the backbone")
        keep_tail = ""
        head_len = 0
    else:
        keep_tail, head_len = rtt[len(cargo) :], len(cargo)
    spacer, rcpbs2, pbs1 = seq[: regions["spacer"][1]], seq[regions["rcpbs2"][0] : regions["rcpbs2"][1]], seq[regions["pbs1"][0] : regions["pbs1"][1]]
    scaffold = seq[regions["scaffold"][0] : regions["scaffold"][1]]
    tail = seq[regions["pbs1"][1] :]

    def judge(new_insert: str):
        """Cross-region pairs for this backbone with ``new_insert`` (sense strand) as the cargo."""
        new_rtt = rc(new_insert) + keep_tail
        peg = spacer + scaffold + rcpbs2 + new_rtt + pbs1 + tail
        regs = _regions_for(spacer, rcpbs2, new_rtt, pbs1)
        pairs = scorer.pair_probs(peg)
        score = functional_region_pairing(peg, regions=regs, scorer=scorer, threshold=thr, pairs=pairs)
        return score.cross_region_pairs, len(peg)

    rng = random.Random(20261104)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    print(f"\n[3] SCALE: {args.n} random inserts of each length in this backbone (same spacer, rcPBS2, PBS1, arm)")
    print("    insert (nt)  pegRNA (nt)   designs with >= 1 confident cross-region pair")
    if head_len:
        for k in (24, 60, 120, 240, 360, 501):
            res = [judge(rand(k)) for _ in range(args.n)]
            bad = sum(1 for c, _ in res if c > 0)
            print(f"    {k:11d}  {res[0][1]:11d}   {bad}/{args.n} rejected  (pairs: {sorted(c for c, _ in res)})")

    print("\n[4] NULL and controls")
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[1].forward_table
    by_aa: dict = {}
    for codon, aa in table.items():
        by_aa.setdefault(aa, []).append(codon)

    def shuffle(cds):
        return "".join(rng.choice(by_aa[table[cds[i : i + 3]]]) if cds[i : i + 3] in table else cds[i : i + 3] for i in range(0, len(cds) - len(cds) % 3, 3))

    if head_len:
        syn = [judge(shuffle(cargo)) for _ in range(args.n)]
        print(f"    random synonymous recodings of this cassette: {sum(1 for c, _ in syn if c > 0)}/{args.n} rejected (pairs {sorted(c for c, _ in syn)})")
        mine, _ = judge(cargo)
        print(f"    this cassette, as built                       : {mine} cross-region pair(s)")
        if Path(args.flag_config).exists():
            flag = json.loads(Path(args.flag_config).read_text())["cargo"].upper()
            c, n = judge(flag)
            print(f"    the shipped {len(flag)}-bp FLAG, same backbone ({n} nt)   : {c} cross-region pair(s)")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
