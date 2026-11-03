"""Measure the docking gate on REAL pipeline candidates, not synthetic ones.

    python scripts/analysis/sample_pool_candidates.py --dry-run     # what would be measured
    python scripts/analysis/sample_pool_candidates.py               # measure (spends NIM credits)

The docking reference was built from random and scrambled extensions. The question that matters
for "is the gate too harsh?" is what happens to the pegRNAs the pipeline actually produces. This
script runs the real candidate pool (recode, INV-0, rcPBS2, and the calibrated ViennaRNA folding
gate, with a docking gate that records every pegRNA it is shown and rejects it, so the walk visits
the whole pool instead of stopping at the first winner), then measures a sample with Boltz-2.

Records go to ``data/docking_calibration.jsonl`` with class ``candidate``. They are NOT part of the
reference table: they are what it is tested on. A run skips any RNA already measured.

Sampling. For each (spacer, cargo) pair: the first ``--head`` candidates in walk order, which is
what production judges first and so is the most direct evidence about the gate, plus a seeded
random spread of the rest for diversity. Cargos come from the motif library the UI offers, so the
pegRNAs span the lengths users will actually generate.

Limit, stated plainly: one locus (CTCF) and two spacers. Cargo, nick (so arm length and rcPBS2)
and PBS length vary; the genomic context does not. The sample includes nicks the real off-target
screen would drop for CTCF, standing in for loci where they are clean.
"""

import argparse
import concurrent.futures as cf
import json
import random
import statistics as st
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate_docking_gate as cal  # noqa: E402

from pegrna_pipeline.agent.motif_library import MotifLibrary  # noqa: E402
from pegrna_pipeline.hpc import recoded_tjpe_workflow as wf  # noqa: E402
from pegrna_pipeline.hpc.complex_dock import DockFake, DockingServiceError  # noqa: E402
from pegrna_pipeline.hpc.complex_dock_nim import Boltz2NIMDockScorer  # noqa: E402
from pegrna_pipeline.hpc.scorer_factory import make_scorers  # noqa: E402
from pegrna_pipeline.hpc.vienna_structure import ViennaRNAStructureScorer  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
LOCUS_CONFIG = REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json"
LOCUS = "CTCF-cterm"
SPACERS = ("TCAGCATGATGGACCGGTGA", "GCACAAGGCTCCGCCATCAC")

#: (name, motif names in order): realistic C-terminal cargos from the UI's library. All in frame.
CARGOS = (
    ("FLAG", ("FLAG",)),
    ("HA", ("HA",)),
    ("V5", ("V5",)),
    ("GGGGS-HA", ("GGGGS", "HA")),
    ("XTEN-FLAG", ("XTEN", "FLAG")),
    ("FLAG-GGGGS-V5", ("FLAG", "GGGGS", "V5")),
    ("FLAG-P2A", ("FLAG", "P2A")),
    ("V5-P2A", ("V5", "P2A")),
)
#: The second spacer is sampled on a subset (the first spacer carries the full cargo range).
SECOND_SPACER_CARGOS = ("FLAG", "V5", "FLAG-GGGGS-V5", "FLAG-P2A")
BACKBONE_LEN = 96


def build_cargo(parts) -> str:
    lib = MotifLibrary.load_full()
    return "".join(lib.lookup(p).dna for p in parts)


def collect_pool(cfg: dict, spacer: str, cargo: str) -> list:
    """Every distinct pegRNA the walk sends to the docking gate, in walk order."""
    seen: list = []

    def recording_gate(peg, cas9, *, scorer=None, max_drop=None, metrics=None):  # noqa: ARG001
        seen.append(peg)
        return False, "recorded"  # reject, so the walk continues through the whole pool

    fake = make_scorers(mode="fake")
    # Every second nick 5' of S is allowed through, not only the one the real off-target screen
    # keeps for CTCF (nick 101): they stand in for other loci where those positions are clean, and
    # give the pegRNAs real variety in arm length and rcPBS2. Off-target is not what is measured.
    offtarget = fake.offtarget_scorer

    original = wf.passes_docking_gate
    wf.passes_docking_gate = recording_gate
    try:
        try:
            wf.run_recoded_tjpe(
                spacer=spacer,
                rtt_aa=cfg["rtt_aa"],
                cargo=cargo,
                genome=cfg["genome"],
                site_s=cfg["site_s"],
                cds_mask=cfg["cds_mask"],
                pbs_lengths=tuple(cfg["pbs_lengths"]),
                ngrna_distance_min=cfg["ngrna_distance_min"],
                ngrna_distance_max=cfg["ngrna_distance_max"],
                pegrna_strand=cfg["pegrna_strand"],
                arm_side=cfg["arm_side"],
                max_run=cfg["max_run"],
                splice_scorer=fake.splice_scorer,
                full_structure_scorer=ViennaRNAStructureScorer(),
                importance_scorer=fake.importance_scorer,
                dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.1),
                offtarget_scorer=offtarget,
                maxent_scorer=fake.maxent_scorer,
                protein_length=cfg.get("protein_length"),
                max_insert_variants=1,  # no synonymous-cargo retries: we want THIS cargo's pool
            )
        except wf.NoValidPegRNA:
            pass  # expected: every candidate was rejected by the recording gate
    finally:
        wf.passes_docking_gate = original
    return list(dict.fromkeys(seen))


def choose(pool: list, head: int, total: int, seed_parts: tuple) -> list:
    """(index in walk order, sequence): the first ``head``, then a seeded spread of the rest,
    spread over the length range so the sample is not all one PBS length."""
    picked = list(range(min(head, len(pool))))
    rest = list(range(len(picked), len(pool)))
    rng = random.Random(cal.rng(0, *seed_parts).random())
    rng.shuffle(rest)
    picked += sorted(rest[: max(0, total - len(picked))])
    return [(i, pool[i]) for i in picked]


def build_jobs(args) -> list:
    cfg = json.loads(LOCUS_CONFIG.read_text())
    jobs = []
    cargos = dict((n, build_cargo(p)) for n, p in CARGOS)
    for spacer in SPACERS:
        for name, cargo in cargos.items():
            if spacer != SPACERS[0] and name not in SECOND_SPACER_CARGOS:
                continue
            pool = collect_pool(cfg, spacer, cargo)
            for rank, rna in choose(pool, args.head, args.per_combo, (spacer, name)):
                jobs.append(
                    {
                        "source": "pool",
                        "class": "candidate",
                        "spacer": spacer,
                        "ext_len": len(rna) - BACKBONE_LEN,
                        "rep": rank,
                        "segment": None,
                        "rna": rna,
                        "locus": LOCUS,
                        "cargo": name,
                        "cargo_len": len(cargo),
                        "pool_rank": rank,
                        "pool_size": len(pool),
                        "head": rank < args.head,
                    }
                )
    return jobs


def measure(jobs: list, out: Path, samples: int, workers: int) -> list:
    key = cal.api_key()
    cas9 = Boltz2NIMDockScorer(api_key=key)._resolve_cas9("nCas9")
    lock = threading.Lock()
    failures: list = []

    def run(job):
        client = Boltz2NIMDockScorer(api_key=key, samples=samples)
        t0 = time.time()
        try:
            vals = client._predict_samples(job["rna"], cas9)
            if client.last_metric_source != "pair_chains_iptm_scores":
                raise DockingServiceError(f"metric was {client.last_metric_source}")
        except DockingServiceError as exc:
            failures.append(job["rna"])
            print(f"  FAILED {job['cargo']} len {len(job['rna'])}: {str(exc)[:120]}", flush=True)
            return
        rec = {
            **job,
            "len": len(job["rna"]),
            "samples": [round(float(v), 4) for v in vals],
            "mean": round(st.fmean(vals), 4),
            "scaffold_paired": None,
            "seconds": round(time.time() - t0, 1),
        }
        with lock, open(out, "a") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        print(
            f"  {job['cargo']:14s} rank {job['pool_rank']:3d}/{job['pool_size']:<3d} "
            f"len={rec['len']:3d} mean={rec['mean']:.3f} ({rec['seconds']:.0f}s)",
            flush=True,
        )

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(run, jobs))
    return failures


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=cal.DEFAULT_OUT)
    ap.add_argument("--head", type=int, default=2, help="first N in walk order, per spacer+cargo")
    ap.add_argument("--per-combo", type=int, default=5, help="candidates per spacer+cargo")
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    jobs = build_jobs(args)
    done = cal.load_done(args.out, args.samples)
    todo = [j for j in jobs if j["rna"] not in done]
    lens = [len(j["rna"]) for j in jobs]
    print(
        f"{len(jobs)} candidates sampled ({len(jobs) - len(todo)} already measured, "
        f"{len(todo)} to run); lengths {min(lens)}-{max(lens)} nt"
    )
    by = {}
    for j in jobs:
        by.setdefault((j["spacer"][:8], j["cargo"]), []).append(len(j["rna"]))
    for (sp, cargo), ls in by.items():
        pool = next(j["pool_size"] for j in jobs if (j["spacer"][:8], j["cargo"]) == (sp, cargo))
        print(
            f"   {sp}..  {cargo:14s} pool {pool:3d} reaches docking; sampled {len(ls)}: "
            f"{min(ls)}-{max(ls)} nt"
        )
    if args.dry_run or not todo:
        return 0
    failures = measure(todo, args.out, args.samples, args.workers)
    print(f"\n{len(todo) - len(failures)} measured, {len(failures)} failed -> {args.out}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
