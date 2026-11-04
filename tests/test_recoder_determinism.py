"""The recoder must give the same arm for the same input, in every process (issue #30).

DnaChisel resolves constraints in the order it is given them. The recoder built that order from a ``set`` of
k-mer strings, and a ``set`` of strings iterates in a per-process order (string hashing is randomised), so the
same 91-nt CTCF arm came out as 17, 14 or 15 edits, and as several different arms, depending only on
``PYTHONHASHSEED``. Every result was valid; none could be regenerated.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from pegrna_pipeline.hpc.recode_arm import _avoid_kmers

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json"

# the PYTHONHASHSEED values that disagreed in the bug report, with ones that did not
HASH_SEEDS = (0, 1, 2, 5, 9, 11)

_PROBE = """
import hashlib, json, sys
from pegrna_pipeline.hpc.recode_arm import recode_arm
from pegrna_pipeline.hpc.recoded_tjpe_workflow import _frame_at
cfg = json.load(open(sys.argv[1]))
g, s, mask = cfg["genome"].upper(), cfg["site_s"], cfg["cds_mask"]
out = []
for lo in (101, 70, 56, 55):
    seg = g[lo:s]
    r = recode_arm(seg, genome=seg, frame=_frame_at(mask, lo), max_run=8)
    out.append(f"{lo}:{r.n_edits}:{hashlib.md5(r.recoded.encode()).hexdigest()[:10]}")
print(" ".join(out))
"""


def _run_with_hashseed(seed: int) -> str:
    proc = subprocess.run(
        [sys.executable, "-c", _PROBE, str(EXAMPLE)],
        capture_output=True,
        text=True,
        cwd=str(REPO),
        env={"PYTHONHASHSEED": str(seed), "PATH": "/usr/bin:/bin", "PYTHONPATH": str(REPO / "src")},
        timeout=300,
    )
    assert proc.returncode == 0, proc.stderr[-800:]
    return proc.stdout.strip().splitlines()[-1]


# ---- the ordering, pure Python (runs in CI, which has no DnaChisel) ------------------------


def test_the_avoided_kmers_come_in_sorted_order():
    kmers = _avoid_kmers("ACGTACGTTTGACCA", 3)
    assert kmers == sorted(kmers)


def test_the_avoided_kmers_are_unique_and_the_right_length():
    kmers = _avoid_kmers("AAAAAAAACGT", 3)
    assert len(kmers) == len(set(kmers))
    assert all(len(k) == 4 for k in kmers)


def test_every_window_of_the_genome_is_represented():
    genome = "ACGTTGCAAC"
    kmers = _avoid_kmers(genome, 3)
    assert {genome[i : i + 4] for i in range(len(genome) - 3)} == set(kmers)


def test_the_order_does_not_depend_on_the_input_order_of_equal_content():
    a = _avoid_kmers("GATTACAGATTACA", 4)
    b = _avoid_kmers("GATTACAGATTACA", 4)
    assert a == b


# ---- the recoder itself (needs DnaChisel; the CI extras do not install it) ------------------


@pytest.fixture(scope="module")
def dnachisel():
    return pytest.importorskip("dnachisel")


def test_the_long_ctcf_arms_are_identical_in_every_process(dnachisel):
    results = {seed: _run_with_hashseed(seed) for seed in HASH_SEEDS}
    assert len(set(results.values())) == 1, results


def test_the_recode_is_no_worse_than_the_best_the_bug_ever_found(dnachisel):
    """The 91-nt arm gave 17 edits typically and 14 at best under the old ordering. A multi-start over
    fixed orderings must reach the better answer every time, not only by luck."""
    line = _run_with_hashseed(0)
    edits = {int(p.split(":")[0]): int(p.split(":")[1]) for p in line.split()}
    assert edits[56] <= 14 and edits[55] <= 14, line
    assert edits[101] <= 9, line  # the 46-nt demo arm


def test_the_result_still_satisfies_the_contract(dnachisel):
    """Protein-identical and no run longer than max_run: determinism must not cost either."""
    from Bio.Seq import Seq

    from pegrna_pipeline.hpc.microhomology import longest_common_run
    from pegrna_pipeline.hpc.recode_arm import recode_arm
    from pegrna_pipeline.hpc.recoded_tjpe_workflow import _frame_at

    cfg = json.loads(EXAMPLE.read_text())
    genome, site_s, mask = cfg["genome"].upper(), cfg["site_s"], cfg["cds_mask"]
    for lo in (101, 70, 56):
        seg = genome[lo:site_s]
        frame = _frame_at(mask, lo)
        rec = recode_arm(seg, genome=seg, frame=frame, max_run=8)
        assert rec is not None and len(rec.recoded) == len(seg)
        pad = (3 - frame) % 3
        a, b = seg[pad:], rec.recoded[pad:]
        a, b = a[: len(a) // 3 * 3], b[: len(b) // 3 * 3]
        assert str(Seq(a).translate()) == str(Seq(b).translate())
        assert longest_common_run(rec.recoded, seg) <= 8
