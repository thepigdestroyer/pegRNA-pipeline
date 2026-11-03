"""Derive ``gate_thresholds.DOCKING_REFERENCE`` from the measurements, and test the gate on them.

    python scripts/analysis/docking_reference.py            # tables and validation
    python scripts/analysis/docking_reference.py --emit     # the constant, ready to paste
    python scripts/analysis/docking_reference.py --check    # exit 1 when DOCKING_REFERENCE is stale

Input: ``data/docking_calibration.jsonl``, written by ``calibrate_docking_gate.py``. The table in
the source tree is meant to be a function of that file, and ``tests/test_docking_reference.py``
fails when it stops being one, so a number can never drift away from its evidence.

What is measured and how it is judged
-------------------------------------
Each record is one RNA scored by Boltz-2 (mean interface ipTM over 5 diffusion samples). The
reference table has one row per extension length with two columns:

  neutral    mean over ``neutral`` / ``neutral-shuffled`` / ``control-canonical`` sequences
  nonbinder  mean over ``nonbinder`` / ``control-scrambled`` sequences

``sequestering`` (an extension built to pair with part of the scaffold) and ``designed`` (the
pipeline's own pegRNAs) are NOT in the table: they are what the table is tested on.

Honest evaluation. A table built from some sequences and then used to judge those same sequences
proves little. So:

* leave-one-out: each neutral / non-binder sequence is judged by a table rebuilt WITHOUT it;
* sources listed in ``HELD_OUT_SOURCES`` (a second spacer) never enter the table at all.
"""

import argparse
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

from pegrna_pipeline.hpc.complex_dock import docking_reference
from pegrna_pipeline.hpc.gate_thresholds import (
    DOCKING_CALIBRATED_RANGE,
    DOCKING_REFERENCE,
    THRESHOLDS,
)

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "docking_calibration.jsonl"

#: Sources kept out of the table, so they test it instead of shaping it.
HELD_OUT_SOURCES = frozenset({"spacer2"})
NEUTRAL = frozenset({"neutral", "neutral-shuffled", "control-canonical"})
NONBINDER = frozenset({"nonbinder", "control-scrambled"})
BACKBONE_LEN = 96  # spacer (20) + scaffold (76); total RNA length = 96 + extension
#: A row is only trusted when the neutral-to-non-binder gap is at least this many between-sequence
#: standard deviations. With the cutoff midway, a gap of 4 sd puts each class 2 sd from it, so each
#: is misjudged about 2% of the time if its spread is roughly normal. A narrower row is not a gate.
MIN_SEPARATION_SD = 4.0
MIN_SEQUENCES_PER_CLASS = 5


def load(path: Path = DATA) -> list:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def kind_of(record: dict):
    if record["class"] in NEUTRAL:
        return "neutral"
    if record["class"] in NONBINDER:
        return "nonbinder"
    return None


def reference_records(records: list) -> list:
    """The records the table is built from."""
    return [
        r
        for r in records
        if r["source"] not in HELD_OUT_SOURCES and kind_of(r) and r.get("ext_len") is not None
    ]


def cells(records: list) -> dict:
    out = defaultdict(list)
    for r in records:
        out[(kind_of(r), r["ext_len"])].append(r["mean"])
    return out


def build_table(records: list) -> list:
    """``[(total_length, neutral, nonbinder), ...]``; a row needs both columns measured."""
    c = cells(records)
    exts = sorted({e for (k, e) in c if (("neutral", e) in c and ("nonbinder", e) in c)})
    return [
        (
            BACKBONE_LEN + e,
            round(st.fmean(c[("neutral", e)]), 3),
            round(st.fmean(c[("nonbinder", e)]), 3),
        )
        for e in exts
    ]


def row_quality(records: list) -> list:
    """Per table row: ``(length, span, sd, separation, n_neutral, n_nonbinder)``, where ``sd`` is
    the larger of the two classes' between-sequence spread and ``separation = span / sd``."""
    ref = reference_records(records)
    c = cells(ref)
    out = []
    for length, neutral, nonbinder in build_table(ref):
        e = length - BACKBONE_LEN
        n, b = c[("neutral", e)], c[("nonbinder", e)]
        sd = max(_sd(n), _sd(b))
        span = neutral - nonbinder
        out.append((length, span, sd, span / sd if sd > 0 else float("inf"), len(n), len(b)))
    return out


def calibrated_range(records: list, fraction: float = THRESHOLDS.docking_pass_fraction) -> tuple:
    """``(lowest, highest)`` total RNA length the gate may judge, derived from the data.

    Walks up the table and stops at the first row that is not trustworthy: too few sequences,
    a gap under ``MIN_SEPARATION_SD`` between-sequence sd, or any leave-one-out error. The range
    ends at the last trustworthy row; past it the gate reports and does not decide."""
    errors = defaultdict(int)
    for r, _f, ok in leave_one_out(records, fraction):
        errors[r["ext_len"] + BACKBONE_LEN] += ok != (kind_of(r) == "neutral")
    top = BACKBONE_LEN
    for length, _span, _sd_, sep, n_neutral, n_nonbinder in row_quality(records):
        trusted = (
            n_neutral >= MIN_SEQUENCES_PER_CLASS
            and n_nonbinder >= MIN_SEQUENCES_PER_CLASS
            and sep >= MIN_SEPARATION_SD
            and errors[length] == 0
        )
        if not trusted:
            break
        top = length
    return (BACKBONE_LEN, top)


def evaluate(record: dict, table, fraction: float):
    """(fraction of the way from non-binder to neutral, passes the gate) for one record."""
    ref = docking_reference(record["len"], table)
    f = (record["mean"] - ref.nonbinder) / ref.span
    return f, record["mean"] >= ref.cutoff(fraction) - 1e-9


def leave_one_out(records: list, fraction: float) -> list:
    """(record, f, passes) for every reference record, judged by a table built without it.
    Records alone in their cell cannot be evaluated this way and are skipped."""
    ref = reference_records(records)
    out = []
    for i, r in enumerate(ref):
        rest = ref[:i] + ref[i + 1 :]
        table = build_table(rest)
        if r["ext_len"] + BACKBONE_LEN not in {row[0] for row in table}:
            continue
        f, ok = evaluate(r, table, fraction)
        out.append((r, f, ok))
    return out


def _sd(values):
    return st.pstdev(values) if len(values) > 1 else float("nan")


def fraction_sweep(loo: list, fractions) -> list:
    """For each candidate pass fraction: neutral rejected, non-binder accepted (leave-one-out)."""
    neutral = [f for r, f, _ in loo if kind_of(r) == "neutral"]
    nonbinder = [f for r, f, _ in loo if kind_of(r) == "nonbinder"]
    return [
        (
            x,
            sum(f < x for f in neutral),
            len(neutral),
            sum(f >= x for f in nonbinder),
            len(nonbinder),
        )
        for x in fractions
    ]


def report(records: list, fraction: float) -> None:
    ref = reference_records(records)
    c = cells(ref)
    table = build_table(ref)

    print(
        f"{len(records)} records; {len(ref)} in the reference "
        f"(held out: {sorted(HELD_OUT_SOURCES)})\n"
    )
    print("REFERENCE TABLE   (total nt, neutral, non-binder  | n, between-sequence sd per column)")
    for length, neutral, nonbinder in table:
        e = length - BACKBONE_LEN
        nn, nb = c[("neutral", e)], c[("nonbinder", e)]
        print(
            f"  {length:4d}   {neutral:.3f}  {nonbinder:.3f}   span {neutral - nonbinder:.3f}"
            f"  | n={len(nn):2d} sd {_sd(nn):.3f}   n={len(nb):2d} sd {_sd(nb):.3f}"
        )

    print(
        "\nROW QUALITY   (gap between the classes, in between-sequence standard deviations; "
        f"a row is trusted from {MIN_SEPARATION_SD:g})"
    )
    for length, span, sd, sep, nn, nb in row_quality(records):
        print(
            f"  {length:4d}   span {span:.3f}   sd {sd:.3f}   separation {sep:4.1f}   "
            f"{'trusted' if sep >= MIN_SEPARATION_SD else 'NOT trusted'}   (n={nn}/{nb})"
        )
    derived = calibrated_range(records, fraction)
    print(
        f"  calibrated range derived from the data: {derived[0]}-{derived[1]} nt   "
        f"(in source: {DOCKING_CALIBRATED_RANGE[0]}-{DOCKING_CALIBRATED_RANGE[1]})"
    )

    loo_all = leave_one_out(records, fraction)
    in_range_top = derived[1]
    loo = [x for x in loo_all if x[0]["len"] <= in_range_top]
    beyond = [x for x in loo_all if x[0]["len"] > in_range_top]
    print(
        f"\nLEAVE-ONE-OUT   (pass fraction {fraction:g}; each sequence judged by a table built "
        f"without it; the statistics below are for the trusted range, up to {in_range_top} nt)"
    )
    per = defaultdict(lambda: [0, 0])
    for r, _f, ok in loo:
        key = (kind_of(r), r["ext_len"] + BACKBONE_LEN)
        per[key][0] += ok
        per[key][1] += 1
    for kind, want in (("neutral", "accepted"), ("nonbinder", "rejected")):
        good = bad = 0
        line = []
        for (k, length), (acc, n) in sorted(per.items(), key=lambda kv: kv[0][1]):
            if k != kind:
                continue
            right = acc if kind == "neutral" else n - acc
            good += right
            bad += n - right
            line.append(f"{length}:{right}/{n}")
        print(f"  {kind:9s} correctly {want}: {good}/{good + bad}   [{'  '.join(line)}]")
    fn = [f for r, f, _ in loo if kind_of(r) == "neutral"]
    fb = [f for r, f, _ in loo if kind_of(r) == "nonbinder"]
    if fn and fb:
        print(
            f"  position on the 0 (non-binder) .. 1 (neutral) scale: neutral mean "
            f"{st.fmean(fn):.2f} sd {_sd(fn):.2f} (lowest {min(fn):.2f});  "
            f"non-binder mean {st.fmean(fb):.2f} sd {_sd(fb):.2f} (highest {max(fb):.2f})"
        )

    if beyond:
        wrong = sum(ok != (kind_of(r) == "neutral") for r, _f, ok in beyond)
        print(
            f"  BEYOND the trusted range ({len(beyond)} sequences): {wrong} judged wrongly; "
            "shown for the record, not used for any statistic here"
        )
    print(
        "\nPASS FRACTION SWEEP  (leave-one-out; the fraction is where on that scale the cutoff sits)"
    )
    print("   fraction   neutral wrongly rejected   non-binder wrongly accepted")
    for x, rj, nn, ac, nb in fraction_sweep(loo, [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70]):
        mark = "  <-- in use" if abs(x - fraction) < 1e-9 else ""
        print(
            f"   {x:.2f}       {rj:3d}/{nn:<3d} ({100 * rj / max(nn, 1):4.1f}%)"
            f"           {ac:3d}/{nb:<3d} ({100 * ac / max(nb, 1):4.1f}%){mark}"
        )

    print(
        "\nBORDERLINE BAND  (leave-one-out; a verdict within +-margin of the cutoff, on the 0..1 "
        "scale, is flagged)"
    )
    print(
        "   margin   neutral flagged (still accepted)   non-binder flagged (still rejected)"
        "   clearly decided"
    )
    for m in (0.05, 0.10, 0.15, 0.20, 0.25):
        nf = sum(abs(f - fraction) < m for f in fn)
        bf = sum(abs(f - fraction) < m for f in fb)
        clear = 1 - (nf + bf) / max(len(fn) + len(fb), 1)
        mark = "  <-- in use" if abs(m - THRESHOLDS.docking_borderline_margin) < 1e-9 else ""
        print(
            f"   {m:.2f}      {nf:3d}/{len(fn):<3d} ({100 * nf / max(len(fn), 1):4.1f}%)"
            f"                  {bf:3d}/{len(fb):<3d} ({100 * bf / max(len(fb), 1):4.1f}%)"
            f"                {100 * clear:4.0f}%{mark}"
        )

    held = [r for r in records if r["source"] in HELD_OUT_SOURCES and kind_of(r)]
    if held:
        print(f"\nHELD OUT  (never in the table; sources {sorted(HELD_OUT_SOURCES)})")
        for kind, want in (("neutral", "accepted"), ("nonbinder", "rejected")):
            rows = [(r, *evaluate(r, table, fraction)) for r in held if kind_of(r) == kind]
            right = sum(ok if kind == "neutral" else not ok for _, _, ok in rows)
            detail = "  ".join(f"{r['len']}nt f={f:.2f}" for r, f, _ in rows)
            print(f"  {kind:9s} correctly {want}: {right}/{len(rows)}   [{detail}]")

    seq = [r for r in records if r["class"] == "sequestering"]
    if seq:
        print(
            "\nSEQUESTERING  (extension reverse-complementary to a scaffold segment; scaffold = 0..76)"
        )
        print(
            "   segment    n   mean f   rejected   scaffold nt paired to the extension (ViennaRNA MFE)"
        )
        by = defaultdict(list)
        for r in seq:
            by[(r.get("segment"), r["len"])].append(r)
        for (seg, length), rs in sorted(by.items(), key=lambda kv: (kv[0][1], str(kv[0][0]))):
            fs = [evaluate(r, table, fraction) for r in rs]
            paired = [r.get("scaffold_paired") for r in rs if r.get("scaffold_paired") is not None]
            print(
                f"   {str(seg):9s} {len(rs):2d}   {st.fmean(f for f, _ in fs):6.2f}   "
                f"{sum(not ok for _, ok in fs)}/{len(rs)}     "
                f"{('~' + str(round(st.fmean(paired), 1))) if paired else 'n/a'}   ({length} nt)"
            )

    designed = [r for r in records if r["class"] == "designed"]
    if designed:
        print("\nTHE PIPELINE'S OWN pegRNAs")
        for r in designed:
            f, ok = evaluate(r, table, fraction)
            print(
                f"   {r['len']} nt   mean ipTM {r['mean']:.3f}   f={f:.2f}   "
                f"{'accepted' if ok else 'REJECTED'}"
            )

    candidates = [r for r in records if r["class"] == "candidate"]
    if candidates:
        report_candidates(candidates, table, fraction, derived[1])

    controls = [r for r in records if r["class"] == "control-scaffold-only"]
    for r in controls:
        print(
            f"\n(scaffold alone, {r['len']} nt: mean {r['mean']:.3f}; outside the calibrated range)"
        )


OLD_FLOOR = 0.6  # the interface-ipTM floor the gate used before it was calibrated


def candidate_summary(candidates: list, table, fraction: float, top: int) -> dict:
    """What the gate does to the pegRNAs the pipeline itself produces."""
    inside = [r for r in candidates if r["len"] <= top]
    beyond = [r for r in candidates if r["len"] > top]
    judged = [(r, *evaluate(r, table, fraction)) for r in inside]
    margin = THRESHOLDS.docking_borderline_margin
    return {
        "n": len(candidates),
        "inside": len(inside),
        "beyond": len(beyond),
        "accepted": sum(ok for _, _f, ok in judged),
        "borderline": sum(abs(f - fraction) < margin for _, f, _ok in judged),
        "old_floor_pass": sum(r["mean"] >= OLD_FLOOR for r in candidates),
        "old_floor_pass_inside": sum(r["mean"] >= OLD_FLOOR for r in inside),
        "f_values": sorted(f for _, f, _ok in judged),
        "head": [(r, f, ok) for r, f, ok in judged if r.get("head")],
    }


def report_candidates(candidates: list, table, fraction: float, top: int) -> None:
    c = candidate_summary(candidates, table, fraction, top)
    print(
        f"\nREAL POOL CANDIDATES   (the pegRNAs the pipeline itself produces: {c['n']} measured; "
        "not in the table, so this is a test of it)"
    )
    print(
        f"  inside the trusted range (<= {top} nt): {c['inside']}   accepted by the gate: "
        f"{c['accepted']}/{c['inside']}   borderline: {c['borderline']}"
    )
    if c["f_values"]:
        fv = c["f_values"]
        print(
            f"     position on the 0 (non-binder) .. 1 (neutral) scale: lowest {fv[0]:.2f}, "
            f"median {st.median(fv):.2f}, highest {fv[-1]:.2f}"
        )
    print(
        f"  beyond it: {c['beyond']} (reported, not judged: they pass as advisory-only)   "
        f"| passing the OLD 0.6 floor: {c['old_floor_pass']}/{c['n']} overall, "
        f"{c['old_floor_pass_inside']}/{c['inside']} inside the range"
    )
    head = c["head"]
    if head:
        first = sum(ok for _, _f, ok in head)
        print(
            f"  the first candidates production would try (walk order): accepted {first}/{len(head)}"
        )
    by = defaultdict(list)
    for r in candidates:
        by[r.get("cargo", "?")].append(r)
    print("   cargo            n   length (nt)   mean ipTM range   position (inside the range)")
    for cargo, rs in sorted(by.items(), key=lambda kv: min(r["len"] for r in kv[1])):
        lens = [r["len"] for r in rs]
        ms = [r["mean"] for r in rs]
        fs = [evaluate(r, table, fraction)[0] for r in rs if r["len"] <= top]
        pos = f"{min(fs):.2f}..{max(fs):.2f}" if fs else "-"
        print(
            f"   {cargo:14s} {len(rs):3d}   {min(lens):3d}-{max(lens):3d}      "
            f"{min(ms):.2f}-{max(ms):.2f}          {pos}"
        )


def emit(records: list) -> str:
    ref = reference_records(records)
    c = cells(ref)
    rows = []
    for length, neutral, nonbinder in build_table(ref):
        e = length - BACKBONE_LEN
        rows.append(
            f"    ({length}, {neutral:.3f}, {nonbinder:.3f}),"
            f"  # n = {len(c[('neutral', e)])} neutral, {len(c[('nonbinder', e)])} non-binder"
        )
    lo, hi = calibrated_range(records)
    return (
        "DOCKING_REFERENCE = (\n" + "\n".join(rows) + "\n)\n\n"
        f"DOCKING_CALIBRATED_RANGE = ({lo}, {hi})"
    )


def stale(records: list) -> list:
    """Differences between the committed constant and the one the data implies."""
    want = build_table(reference_records(records))
    have = list(DOCKING_REFERENCE)
    problems = []
    if [w[0] for w in want] != [h[0] for h in have]:
        problems.append(
            f"lengths differ: data implies {[w[0] for w in want]}, "
            f"source has {[h[0] for h in have]}"
        )
        return problems
    for w, h in zip(want, have):
        for name, a, b in (("neutral", w[1], h[1]), ("non-binder", w[2], h[2])):
            if abs(a - b) > 0.0006:
                problems.append(f"{w[0]} nt {name}: data {a:.3f}, source {b:.3f}")
    want_range = calibrated_range(records)
    if tuple(DOCKING_CALIBRATED_RANGE) != want_range:
        problems.append(
            f"DOCKING_CALIBRATED_RANGE: data implies {want_range}, source has "
            f"{tuple(DOCKING_CALIBRATED_RANGE)}"
        )
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path, default=DATA)
    ap.add_argument("--fraction", type=float, default=THRESHOLDS.docking_pass_fraction)
    ap.add_argument("--emit", action="store_true", help="print the DOCKING_REFERENCE constant")
    ap.add_argument("--check", action="store_true", help="fail if gate_thresholds is stale")
    args = ap.parse_args()
    records = load(args.data)
    if args.emit:
        print(emit(records))
        return 0
    if args.check:
        problems = stale(records)
        for p in problems:
            print("STALE:", p)
        print("DOCKING_REFERENCE matches the data" if not problems else "")
        return 1 if problems else 0
    report(records, args.fraction)
    return 0


if __name__ == "__main__":
    sys.exit(main())
