"""Task 4: DnaChisel constrained-optimisation recode (skips if dnachisel absent)."""

import pytest
from Bio.Seq import Seq

pytest.importorskip("dnachisel")

from pegrna_pipeline.hpc.microhomology import longest_common_run  # noqa: E402
from pegrna_pipeline.hpc.recode_arm import (  # noqa: E402
    recode_arm,
    recode_min_edits_dnachisel,
)


def test_dnachisel_recode_preserves_protein_and_breaks_run():
    cds_span = "GAAGAAGAAGAAGAAGAA"
    res = recode_min_edits_dnachisel(cds_span, genome=cds_span, frame=0, max_run=6)
    assert res is not None
    assert str(Seq(res.recoded).translate()) == str(Seq(cds_span).translate())
    assert longest_common_run(res.recoded, cds_span) <= 6


def test_dnachisel_noop_when_already_clean():
    cds_span = "GAAGAAGAA"
    res = recode_min_edits_dnachisel(cds_span, genome="TTTTTTTTT", frame=0, max_run=6)
    assert res is not None
    assert res.n_edits == 0


def test_dnachisel_rejects_invalid_frame():
    with pytest.raises(ValueError, match="frame"):
        recode_min_edits_dnachisel("GAAGAA", genome="GAAGAA", frame=3, max_run=3)


def test_recode_arm_dispatch_prefers_dnachisel_when_available():
    # recode_arm() should return a valid RecodeResult via whichever backend
    cds_span = "GAAGAAGAAGAAGAAGAA"
    res = recode_arm(cds_span, genome=cds_span, frame=0, max_run=6)
    assert res is not None
    assert str(Seq(res.recoded).translate()) == "EEEEEE"
    assert longest_common_run(res.recoded, cds_span) <= 6
