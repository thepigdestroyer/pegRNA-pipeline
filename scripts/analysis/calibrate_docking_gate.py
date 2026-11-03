"""Measure the Boltz-2 ipTM references behind the docking gate.

    python scripts/analysis/calibrate_docking_gate.py --out docking_calibration.jsonl
    python scripts/analysis/calibrate_docking_gate.py --per-length 8 --lengths 60,98,136

Spends NVIDIA NIM credits: roughly one request per sequence, each taking 30-60 s. Needs
NVIDIA_API_KEY, or the key file at ~/.config/pegrna/nvidia.key.

For each 3'-extension length it scores, with 5 diffusion samples per request,
  * NEUTRAL      the functional backbone (spacer + scaffold) + a RANDOM extension
  * NON-BINDER   a scrambled backbone + a random extension
  * SEQUESTERING the functional backbone + an extension built to pair with the scaffold
                 (the positive control: if ipTM cannot see this, it cannot see the failure
                 the gate exists to catch)
and prints per-length means and between-sequence sd, which are the rows of
``gate_thresholds.DOCKING_REFERENCE``.

Why this exists: the original gate used the AlphaFold-Multimer rule-of-thumb floor, ipTM >=
0.6. Measured here, ipTM falls with total RNA length regardless of design quality, so no
tjPE-sized pegRNA can reach it. See ADR 0024.
"""

import argparse
import concurrent.futures as cf
import json
import os
import random
import statistics as st
import time

from pegrna_pipeline.hpc.complex_dock_nim import Boltz2NIMDockScorer
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

SPACER = "TCAGCATGATGGACCGGTGA"
BACKBONE = SPACER + STANDARD_SGRNA_SCAFFOLD  # 96 nt


def rc(x: str) -> str:
    return x.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def rnd(n: int, seed: int) -> str:
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def scrambled(seed: int) -> str:
    s = list(BACKBONE)
    random.Random(seed).shuffle(s)
    return "".join(s)


class Client(Boltz2NIMDockScorer):
    """Boltz-2 NIM client returning every sample, with polite retry on 429/5xx."""

    def samples_for(self, rna: str, protein: str) -> list:
        # _predict_samples already retries transient failures with Retry-After backoff.
        return self._predict_samples(rna, protein)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="docking_calibration.jsonl")
    ap.add_argument("--lengths", default="30,60,98,136", help="extension lengths, nt")
    ap.add_argument("--per-length", type=int, default=3, help="independent sequences per cell")
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--workers", type=int, default=2, help="keep low: the NIM rate-limits bursts")
    ap.add_argument("--seed", type=int, default=0, help="offset for the random sequences")
    args = ap.parse_args()

    key = (
        os.environ.get("NVIDIA_API_KEY")
        or open(os.path.expanduser("~/.config/pegrna/nvidia.key")).read().strip()
    )
    cas9 = Client(api_key=key)._resolve_cas9("nCas9")
    lengths = [int(x) for x in args.lengths.split(",")]

    jobs = [
        ("control: canonical sgRNA (no extension)", 0, BACKBONE),
        ("control: scrambled backbone (no extension)", 0, scrambled(1)),
    ]
    for L in lengths:
        for k in range(args.per_length):
            jobs.append((f"neutral ext {L}", k, BACKBONE + rnd(L, args.seed + 1000 * L + k)))
            jobs.append(
                (
                    f"non-binder ext {L}",
                    k,
                    scrambled(args.seed + 50 * L + k) + rnd(L, args.seed + 2000 * L + k),
                )
            )
        # the positive control: pairs with the scaffold, so it should tie it up
        pad = max(0, L - 60)
        jobs.append(
            (
                f"sequestering ext {L}",
                0,
                BACKBONE + rc(STANDARD_SGRNA_SCAFFOLD[:60])[:L] + rnd(pad, args.seed + 3000 * L),
            )
        )

    def run(job):
        label, rep, rna = job
        c = Client(api_key=key, samples=args.samples, timeout=1200)
        t0 = time.time()
        vals = c.samples_for(rna, cas9)
        rec = {
            "label": label,
            "rep": rep,
            "len": len(rna),
            "samples": vals,
            "mean": st.fmean(vals),
            "seconds": round(time.time() - t0, 1),
        }
        with open(args.out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(
            f"{label:44} #{rep}  len={len(rna):3}  mean={rec['mean']:.3f}  "
            f"samples={[round(v, 2) for v in vals]}",
            flush=True,
        )
        return rec

    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        recs = list(ex.map(run, jobs))

    print("\nreference table  (length, mean, between-sequence sd, n)")
    cells = {}
    for r in recs:
        cells.setdefault((r["label"], r["len"]), []).append(r["mean"])
    for (label, L), v in sorted(cells.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        sd = f"{st.pstdev(v):.3f}" if len(v) > 1 else "  n/a"
        print(f"  {L:>4} nt  {label:44} {st.fmean(v):.3f}  sd={sd}  n={len(v)}")


if __name__ == "__main__":
    main()
