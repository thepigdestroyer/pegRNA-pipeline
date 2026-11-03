"""Measure the Boltz-2 interface-ipTM references behind the docking gate.

    python scripts/analysis/calibrate_docking_gate.py --dry-run             # plan only, no API calls
    python scripts/analysis/calibrate_docking_gate.py --per-length 6        # neutral + non-binder grid
    python scripts/analysis/calibrate_docking_gate.py --controls 4          # firm up the 96 nt row
    python scripts/analysis/calibrate_docking_gate.py --sweep               # which scaffold regions matter
    python scripts/analysis/calibrate_docking_gate.py --spacer GCTGCCACAGATGCCGCAAC \\
        --lengths 98,136 --per-length 3 --seed 500000 --source spacer2     # a held-out second spacer

Spends NVIDIA NIM credits: one request per RNA, 30-100 s each. Needs NVIDIA_API_KEY, or the key
file ~/.config/pegrna/nvidia.key.

Every RNA is appended to ``data/docking_calibration.jsonl`` as one record (class, spacer,
extension length, the sequence, the 5 per-sample interface ipTMs, and how many scaffold
nucleotides ViennaRNA pairs with the extension). A run SKIPS any RNA already in that file, so an
interrupted run resumes and a repeat costs nothing. ``docking_reference.py`` turns the file into
``gate_thresholds.DOCKING_REFERENCE`` and tests it; neither is meaningful without the other.

Classes, all built on the functional backbone (spacer + 76 nt scaffold = 96 nt):
  neutral       backbone + a random 3' extension: what a pegRNA that does nothing wrong scores
  nonbinder     a SCRAMBLED backbone + a random extension: an RNA that has lost Cas9
  sequestering  backbone + the reverse complement of one scaffold segment + random padding: an
                extension built to pair with, and tie up, that part of the scaffold (--sweep)
  control-*     no extension: canonical sgRNAs with different spacers, and scrambled backbones

Why: the original gate used the AlphaFold-Multimer rule-of-thumb floor, ipTM >= 0.6. Measured
here, ipTM falls with total RNA length regardless of design quality, so no tjPE-sized pegRNA can
reach it. See ADR 0024.
"""

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import random
import statistics as st
import sys
import threading
import time
from pathlib import Path

from pegrna_pipeline.hpc.complex_dock import DockingServiceError
from pegrna_pipeline.hpc.complex_dock_nim import Boltz2NIMDockScorer
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "data" / "docking_calibration.jsonl"
PRIMARY_SPACER = "TCAGCATGATGGACCGGTGA"
SCAFFOLD = STANDARD_SGRNA_SCAFFOLD  # 76 nt
SPACER_LEN = 20

#: 24-nt windows tiling the scaffold: repeat + GAAA loop, anti-repeat, the nexus / stem loop 1,
#: stem loop 2, and the terminator hairpin (stem loop 3). Each is reverse-complemented into the
#: extension to see which parts of the scaffold the model's interface score is sensitive to.
SWEEP_WINDOWS = ((0, 24), (10, 34), (20, 44), (30, 54), (40, 64), (52, 76))
SWEEP_EXT_LEN = 98


def rc(x: str) -> str:
    return x.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def rng(seed: int, *parts) -> random.Random:
    """An independent, reproducible stream per (seed, parts): no accidental overlap between
    the sequences of different classes or lengths."""
    digest = hashlib.sha256(repr((seed, *parts)).encode()).hexdigest()
    return random.Random(int(digest[:12], 16))


def rnd(n: int, seed: int, *parts) -> str:
    r = rng(seed, "ext", *parts)
    return "".join(r.choice("ACGT") for _ in range(n))


def shuffled(s: str, seed: int, *parts) -> str:
    x = list(s)
    rng(seed, "shuffle", *parts).shuffle(x)
    return "".join(x)


def random_spacer(seed: int, *parts) -> str:
    """A plausible 20-mer: GC 40-65%, no TTTT, no homopolymer run >= 5."""
    r = rng(seed, "spacer", *parts)
    while True:
        s = "".join(r.choice("ACGT") for _ in range(SPACER_LEN))
        gc = (s.count("G") + s.count("C")) / SPACER_LEN
        runs = any(c * 5 in s for c in "ACGT")
        if 0.40 <= gc <= 0.65 and "TTTT" not in s and not runs:
            return s


def scaffold_paired_to_extension(rna: str, spacer_len: int = SPACER_LEN):
    """Scaffold nucleotides that ViennaRNA's MFE structure pairs with the 3' extension, or None
    when ViennaRNA is not installed. Tells us whether a 'sequestering' design really does fold
    onto the part of the scaffold it was built for."""
    try:
        import RNA  # noqa: PLC0415
    except ImportError:
        return None
    structure, _ = RNA.fold(rna.replace("T", "U"))
    stack, partner = [], {}
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            partner[i], partner[j] = j, i
    backbone_end = spacer_len + len(SCAFFOLD)
    return sum(1 for i in range(spacer_len, backbone_end) if partner.get(i, -1) >= backbone_end)


def build_jobs(args) -> list:
    spacer = args.spacer.upper()
    backbone = spacer + SCAFFOLD
    seed = args.seed
    jobs = []

    def add(cls, ext_len, rep, rna, *, segment=None, spacer_=spacer):
        jobs.append(
            {
                "source": args.source,
                "class": cls,
                "spacer": spacer_,
                "ext_len": ext_len,
                "rep": rep,
                "segment": segment,
                "rna": rna,
            }
        )

    # No-extension controls: canonical sgRNAs, each with a different random spacer (so the 96 nt
    # reference does not rest on a single spacer), and scrambled backbones. ``--spacer`` itself is
    # not re-measured here: it is already the first control, from the original calibration.
    for k in range(1, args.controls + 1):
        sp = random_spacer(seed, "control", k)
        add("control-canonical", 0, k, sp + SCAFFOLD, spacer_=sp)
        add("control-scrambled", 0, k, shuffled(backbone, seed, "control", k))

    for length in [int(x) for x in args.lengths.split(",") if x]:
        for k in range(args.per_length):
            add("neutral", length, k, backbone + rnd(length, seed, "neutral", length, k))
            add(
                "nonbinder",
                length,
                k,
                shuffled(backbone, seed, "nonbinder", length, k)
                + rnd(length, seed, "nonbinder", length, k),
            )

    if args.sweep:
        for a, b in SWEEP_WINDOWS:
            for k in range(args.sweep_reps):
                pad = SWEEP_EXT_LEN - (b - a)
                add(
                    "sequestering",
                    SWEEP_EXT_LEN,
                    k,
                    backbone + rc(SCAFFOLD[a:b]) + rnd(pad, seed, "sweep", a, b, k),
                    segment=f"{a}:{b}",
                )
    return jobs


def load_done(path: Path, min_samples: int) -> set:
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                if r.get("rna") and len(r.get("samples", [])) >= min_samples:
                    done.add(r["rna"])
    return done


def api_key() -> str:
    key = os.environ.get("NVIDIA_API_KEY")
    if key:
        return key
    return Path(os.path.expanduser("~/.config/pegrna/nvidia.key")).read_text().strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--lengths", default="30,60,98,136", help="extension lengths, nt")
    ap.add_argument("--per-length", type=int, default=6, help="independent sequences per class")
    ap.add_argument(
        "--controls", type=int, default=0, help="ADDITIONAL no-extension controls per class"
    )
    ap.add_argument("--sweep", action="store_true", help="add the scaffold-segment sweep")
    ap.add_argument("--sweep-reps", type=int, default=2)
    ap.add_argument("--spacer", default=PRIMARY_SPACER)
    ap.add_argument("--samples", type=int, default=5, help="diffusion samples per request")
    ap.add_argument("--workers", type=int, default=2, help="keep low: the NIM rate-limits bursts")
    ap.add_argument("--seed", type=int, default=0, help="change it to draw different sequences")
    ap.add_argument("--source", default="calibration", help="tag stored in every record")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if len(args.spacer) != SPACER_LEN or set(args.spacer.upper()) - set("ACGT"):
        ap.error(f"--spacer must be {SPACER_LEN} nt of ACGT")

    jobs = build_jobs(args)
    done = load_done(args.out, args.samples)
    todo = [j for j in jobs if j["rna"] not in done]
    print(f"{len(jobs)} RNAs planned, {len(jobs) - len(todo)} already measured, {len(todo)} to run")
    by_class = {}
    for j in todo:
        by_class[(j["class"], j["ext_len"], j["segment"])] = (
            by_class.get((j["class"], j["ext_len"], j["segment"]), 0) + 1
        )
    for (cls, ext, seg), n in sorted(by_class.items(), key=lambda kv: (str(kv[0][1]), kv[0][0])):
        print(f"   {n:3d} x {cls:18s} ext {ext:>3}" + (f"  segment {seg}" if seg else ""))
    if args.dry_run or not todo:
        return 0

    key = api_key()
    cas9 = Boltz2NIMDockScorer(api_key=key)._resolve_cas9("nCas9")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    lock = threading.Lock()
    failures = []

    def run(job):
        client = Boltz2NIMDockScorer(api_key=key, samples=args.samples)
        t0 = time.time()
        try:
            vals = client._predict_samples(job["rna"], cas9)
            if client.last_metric_source != "pair_chains_iptm_scores":
                raise DockingServiceError(f"metric was {client.last_metric_source}")
        except DockingServiceError as exc:
            failures.append(job["rna"])
            print(f"  FAILED {job['class']} ext {job['ext_len']}: {str(exc)[:140]}", flush=True)
            return
        paired = (
            scaffold_paired_to_extension(job["rna"], SPACER_LEN)
            if job["class"] in ("neutral", "sequestering")
            else None
        )
        rec = {
            **job,
            "len": len(job["rna"]),
            "samples": [round(float(v), 4) for v in vals],
            "mean": round(st.fmean(vals), 4),
            "scaffold_paired": paired,
            "seconds": round(time.time() - t0, 1),
        }
        with lock:
            with open(args.out, "a") as fh:
                fh.write(json.dumps(rec, sort_keys=True) + "\n")
        seg = f" seg {job['segment']}" if job["segment"] else ""
        extra = f"  paired={paired}" if paired is not None else ""
        print(
            f"  {job['class']:18s} ext {job['ext_len']:>3}{seg:12s} len={rec['len']:3d} "
            f"mean={rec['mean']:.3f}{extra}  ({rec['seconds']:.0f}s)",
            flush=True,
        )

    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        list(ex.map(run, todo))

    print(f"\n{len(todo) - len(failures)} measured, {len(failures)} failed -> {args.out}")
    if failures:
        print("re-run the same command to retry only the failures")
    print("next: python scripts/analysis/docking_reference.py")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
