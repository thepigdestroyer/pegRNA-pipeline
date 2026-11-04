"""Where does the SpliceAI cross-check's number come from for a large insertion, and is it evidence?

A real run of the CTCF example with the 501-bp GGGGS-FLAG-2A-resistance cassette passed the MaxEntScan gate
(the PRIMARY splice gate, ADR 0011) and was then rejected on every candidate by the SpliceAI cross-check:
"splice delta 0.336 >= 0.2". The same example with the 24-bp FLAG passes with 0.150. ADR 0011 made SpliceAI an
ADVISORY cross-check that binds only when a caller supplies real >= 10 kb genomic context, because a 500-bp
non-genomic insertion is out of its distribution. In ``api`` mode the window is N-padded to 10 kb instead, and
the number is applied as a hard gate. This probe measures what that number is.

It needs the real SpliceAI weights and TensorFlow, so it runs where the pipeline's real models are (the DSI box):

    dsi_run_py.sh probe_splice_gate.py /path/to/run/config.json [--flag-config PATH] [--skip-real-context] [--n-null 10]

Parts:

1. CAPTURE. Run the real ``design-recoded`` for the config with a recorder around the splice scorer, so the
   exact (reference, edited) pair the gate saw is kept. Confirms the failure is reproduced (0.336).
2. DECOMPOSE. ``SpliceAIScorer.delta`` is the max of several terms (flank perturbation, the strongest site in
   the replaced span = "loss", the strongest site in the inserted/recoded span = "gain"). Print every term, its
   position, whether the position is in the arm, the cargo or a flank, and the sequence around it.
3. SCALE. The "gain" term is the MAXIMUM of the raw site probability over every inserted base, not a
   difference of two predictions. A maximum over n positions grows with n. Score growing prefixes of the cargo.
4. NULL. Score synonymous-codon shuffles of the same cassette (same protein) and fully random sequences of the
   same length, in the same N-padded setting: how often does an arbitrary insert of this size reach 0.2?
5. REAL CONTEXT. Where Ensembl is reachable, put the real +/-5 kb flanks around the window instead of N and
   score again: does the number survive real genomic context (the case ADR 0011 says should bind)?
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import io
import json
import os
import random
import sys
import time
from pathlib import Path

import numpy as np

PAD = 5000  # SpliceAI's context on each side


# ---- tracks and the decomposition (a labelled copy of SpliceAIScorer.delta) ------------------------


class Tracks:
    def __init__(self, scorer, seq: str, left: str, right: str) -> None:
        padded = left + seq + right
        self.acc, self.don = scorer._tracks(padded)
        trim = (len(padded) - len(self.acc)) // 2
        core = padded[trim : len(padded) - trim] if trim else padded
        self.core = core


def _common_prefix(a: str, b: str) -> int:
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _common_suffix(a: str, b: str, cap: int) -> int:
    n = min(len(a), len(b), cap)
    i = 0
    while i < n and a[len(a) - 1 - i] == b[len(b) - 1 - i]:
        i += 1
    return i


def terms(R: Tracks, A: Tracks):
    """Every term of the delta, as (label, track, value, position, on) with ``on`` = "ref" or "alt"."""
    lr, la = len(R.acc), len(A.acc)
    p = _common_prefix(R.core, A.core)
    s = _common_suffix(R.core, A.core, cap=min(lr, la) - p)
    out = []
    for name, ta, tr in (("acceptor", A.acc, R.acc), ("donor", A.don, R.don)):
        if p > 0:
            d = np.abs(ta[:p] - tr[:p])
            i = int(d.argmax())
            out.append(("left flank, existing site changed", name, float(d[i]), i, "alt"))
        if s > 0:
            d = np.abs(ta[la - s :] - tr[lr - s :])
            i = int(d.argmax())
            out.append(("right flank, existing site changed", name, float(d[i]), la - s + i, "alt"))
        if lr - s > p:
            seg = tr[p : lr - s]
            i = int(seg.argmax())
            out.append(("LOSS: strongest native site in the replaced span", name, float(seg[i]), p + i, "ref"))
        if la - s > p:
            seg = ta[p : la - s]
            i = int(seg.argmax())
            out.append(("GAIN: strongest site anywhere in the inserted/recoded span", name, float(seg[i]), p + i, "alt"))
    return out, p, s


def delta_of(R: Tracks, A: Tracks) -> float:
    t, _, _ = terms(R, A)
    return max([0.0] + [x[2] for x in t])


def gain_of(R: Tracks, A: Tracks) -> float:
    t, _, _ = terms(R, A)
    return max([0.0] + [x[2] for x in t if x[0].startswith("GAIN")])


# ---- part 1: capture --------------------------------------------------------------------------


def capture(config_path: str, out_dir: Path):
    from pegrna_pipeline.agent import cli as cli_mod
    from pegrna_pipeline.hpc import recoded_driver

    pairs, deltas = [], []
    real_make = recoded_driver.make_scorers

    class Recorder:
        def __init__(self, inner) -> None:
            self._inner = inner

        def delta(self, ref_seq: str, alt_seq: str) -> float:
            d = self._inner.delta(ref_seq, alt_seq)
            pairs.append((ref_seq, alt_seq))
            deltas.append(d)
            return d

    def patched(**kw):
        bundle = real_make(**kw)
        return dataclasses.replace(bundle, splice_scorer=Recorder(bundle.splice_scorer))

    recoded_driver.make_scorers = patched
    buf = io.StringIO()
    code = 0
    try:
        with contextlib.redirect_stdout(buf):
            cli_mod.design_recoded_cmd.main(
                ["--config", config_path, "--mode", "api", "--output-dir", str(out_dir), "--progress"],
                standalone_mode=False,
            )
    except SystemExit as exc:
        code = int(exc.code or 0)
    finally:
        recoded_driver.make_scorers = real_make
    progress = []
    for line in buf.getvalue().splitlines():
        if line.startswith("@@P "):
            with contextlib.suppress(ValueError):
                progress.append(json.loads(line[4:]))
    return pairs, deltas, progress, code


# ---- helpers ----------------------------------------------------------------------------------


def where(pos: int, alt: str, cargo: str) -> str:
    c0 = alt.find(cargo)
    if c0 < 0:
        return f"alt[{pos}]"
    c1 = c0 + len(cargo)
    if pos < c0:
        return f"{c0 - pos} nt 5' of the cargo (arm / left flank)"
    if pos < c1:
        return f"cargo +{pos - c0} of {len(cargo)}"
    return f"{pos - c1} nt 3' of the cargo (right flank)"


def context(seq: str, pos: int, track: str) -> str:
    if track == "donor":  # exon | GT...
        return f"{seq[max(0, pos - 6) : pos]}|{seq[pos : pos + 6]}"
    return f"{seq[max(0, pos - 6) : pos + 1]}|{seq[pos + 1 : pos + 7]}"  # ...AG | exon


def synonymous_shuffle(cds: str, rng: random.Random) -> str:
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[1].forward_table
    by_aa: dict = {}
    for codon, aa in table.items():
        by_aa.setdefault(aa, []).append(codon)
    out = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i : i + 3]
        aa = table.get(codon)
        out.append(rng.choice(by_aa[aa]) if aa else codon)
    return "".join(out) + cds[len(cds) - len(cds) % 3 :]


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def find_flanks(transcript: str, window: str):
    """Real +/-5 kb flanks around ``window``, from Ensembl, or (None, reason)."""
    from pegrna_pipeline.hpc.genome_ingest import _ENSEMBL, _real_http_get

    lookup = _real_http_get(f"{_ENSEMBL}/lookup/id/{transcript}", expand=1)
    chrom, start, end, strand = (
        lookup["seq_region_name"],
        int(lookup["start"]),
        int(lookup["end"]),
        int(lookup["strand"]),
    )
    wide = _real_http_get(
        f"{_ENSEMBL}/sequence/region/human/{chrom}:{start - 6000}..{end + 6000}:{strand}"
    ).get("seq", "").upper()
    idx = wide.find(window.upper())
    if idx < 0:
        return None, f"the window is not inside the fetched region (chr{chrom} {start}-{end} strand {strand})"
    left, right = wide[max(0, idx - PAD) : idx], wide[idx + len(window) : idx + len(window) + PAD]
    left, right = "N" * (PAD - len(left)) + left, right + "N" * (PAD - len(right))
    return (left, right), f"chr{chrom}:{start}-{end} strand {strand}; window at +{idx} of the fetched region"


# ---- main ---------------------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--skip-real-context", action="store_true")
    ap.add_argument("--n-null", type=int, default=10)
    ap.add_argument(
        "--flag-config",
        default=str(Path.home() / "pegrna" / "prime-editing" / "examples" / "ctcf_cterm" / "config_autospacer_flag.json"),
        help="the shipped 24-bp FLAG example: its cargo is scored in the same arm and site as a control",
    )
    ap.add_argument("--out", default=str(Path.home() / "pegrna" / "dev-runs" / "splice_probe"))
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    cargo = cfg["cargo"].upper()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print(f"config: cargo {len(cargo)} nt, window {len(cfg['genome'])} nt, S = {cfg['site_s']}, "
          f"arm_side {cfg.get('arm_side')}, transcript {cfg.get('transcript')}")

    # 1. capture ------------------------------------------------------------------------------
    print("\n[1] CAPTURE: the real pipeline path, with the splice scorer recorded")
    pairs, deltas, progress, code = capture(args.config, out_dir)
    cryptic = [p for p in progress if p.get("gate") == "cryptic"]
    for p in cryptic[:1]:
        print(f"    MaxEntScan (primary gate), candidate {p['cand']}: passed={p['passed']}  {p['info']}  {p.get('metrics') or ''}")
    print(f"    the SpliceAI cross-check was called {len(pairs)}x; deltas = {sorted(set(round(d, 3) for d in deltas))}; exit code {code}")
    if not pairs:
        print("    the splice gate was never reached; nothing to decompose")
        return
    ref, alt = pairs[0]
    print(f"    reference {len(ref)} nt, edited {len(alt)} nt (+{len(alt) - len(ref)})")

    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

    scorer = SpliceAIScorer()
    NL = "N" * PAD
    R = Tracks(scorer, ref, NL, NL)
    A = Tracks(scorer, alt, NL, NL)

    # 2. decompose ----------------------------------------------------------------------------
    print("\n[2] DECOMPOSE: every term of the delta (N-padded, as in api mode)")
    t, p, s = terms(R, A)
    print(f"    aligned left flank {p} nt, right flank {s} nt; variant region: ref[{p}:{len(ref) - s}] -> alt[{p}:{len(alt) - s}]")
    print(f"    recomputed delta = {max(x[2] for x in t):.3f}   (recorded {deltas[0]:.3f})")
    for label, track, value, pos, on in sorted(t, key=lambda x: -x[2]):
        seq = alt if on == "alt" else ref
        loc = where(pos, alt, cargo) if on == "alt" else f"ref[{pos}], {pos - cfg['site_s']:+d} from S"
        print(f"    {value:6.3f}  {track:8s} {label}")
        print(f"            at {loc}; sequence {context(seq, pos, track)}")

    # 3. scale --------------------------------------------------------------------------------
    print("\n[3] SCALE: the gain term for growing prefixes of the cargo (same arm, same site)")
    c0 = alt.find(cargo)
    left_part, right_part = alt[:c0], alt[c0 + len(cargo) :]
    print("    insert length   gain   overall delta")
    for k in sorted({24, 60, 120, 240, 360, len(cargo)}):
        if k > len(cargo):
            continue
        Ak = Tracks(scorer, left_part + cargo[:k] + right_part, NL, NL)
        print(f"    {k:11d}   {gain_of(R, Ak):.3f}   {delta_of(R, Ak):.3f}")
    if Path(args.flag_config).exists():
        flag = json.loads(Path(args.flag_config).read_text())["cargo"].upper()
        Af = Tracks(scorer, left_part + flag + right_part, NL, NL)
        print(f"    {len(flag):11d}   {gain_of(R, Af):.3f}   {delta_of(R, Af):.3f}   <- the shipped FLAG example (the gate reported 0.150)")

    # 4. null ---------------------------------------------------------------------------------
    rng = random.Random(20261104)
    print(f"\n[4] NULL: {args.n_null} synonymous-codon shuffles of this cassette, and {args.n_null} random sequences of the same length")
    syn, rnd = [], []
    for _ in range(args.n_null):
        a = Tracks(scorer, left_part + synonymous_shuffle(cargo, rng) + right_part, NL, NL)
        syn.append(delta_of(R, a))
        b = Tracks(scorer, left_part + random_dna(len(cargo), rng) + right_part, NL, NL)
        rnd.append(delta_of(R, b))

    def summary(x):
        x = sorted(x)
        return f"min {x[0]:.3f}  median {x[len(x) // 2]:.3f}  max {x[-1]:.3f}  >= 0.2: {sum(v >= 0.2 for v in x)}/{len(x)}"

    print(f"    same protein, random synonymous codons : {summary(syn)}")
    print(f"    random {len(cargo)}-nt sequences             : {summary(rnd)}")
    print(f"    this cassette                          : {deltas[0]:.3f}")

    # 5. real context -------------------------------------------------------------------------
    if args.skip_real_context or not cfg.get("transcript"):
        print("\n[5] REAL CONTEXT: skipped")
    else:
        print("\n[5] REAL CONTEXT: real +/-5 kb flanks from Ensembl instead of N")
        os.environ.setdefault("PEGRNA_ENSEMBL_RETRIES", "2")
        os.environ.setdefault("PEGRNA_ENSEMBL_TIMEOUT", "60")
        try:
            flanks, note = find_flanks(cfg["transcript"], cfg["genome"])
        except Exception as exc:  # noqa: BLE001 - a probe reports, it does not crash
            flanks, note = None, f"{type(exc).__name__}: {str(exc)[:200]}"
        print(f"    {note}")
        if flanks:
            L, Rr = flanks
            Rt = Tracks(scorer, ref, L, Rr)
            At = Tracks(scorer, alt, L, Rr)
            tt, _, _ = terms(Rt, At)
            print(f"    delta with real flanks: {max(x[2] for x in tt):.3f}   (N-padded: {deltas[0]:.3f})")
            for label, track, value, pos, on in sorted(tt, key=lambda x: -x[2])[:3]:
                seq = alt if on == "alt" else ref
                loc = where(pos, alt, cargo) if on == "alt" else f"ref[{pos}]"
                print(f"      {value:6.3f}  {track:8s} {label}: {loc}; {context(seq, pos, track)}")
            syn_real = []
            for _ in range(min(args.n_null, 6)):
                a = Tracks(scorer, left_part + synonymous_shuffle(cargo, rng) + right_part, L, Rr)
                syn_real.append(delta_of(Rt, a))
            print(f"    synonymous shuffles, real flanks       : {summary(syn_real)}")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
