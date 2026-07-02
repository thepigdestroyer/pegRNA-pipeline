"""Offline tests for the real-scorer fixes (SpliceAI indel delta, Cas-OFFinder
aggregation, RiNALMo->RNAformer structure fallback). The heavy models themselves
are CX3-only; these tests exercise the pure logic with stubbed model calls.
"""

from __future__ import annotations

import numpy as np

from pegrna_pipeline.hpc.full_structure import AutoStructureScorer
from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer
from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer


# --------------------------------------------------------------- SpliceAI indel
class _StubSpliceAI(SpliceAIScorer):
    """SpliceAIScorer with a deterministic, length-addressable _tracks stub so the
    indel-alignment math is testable without TensorFlow/weights. The stub encodes a
    strong donor 'site' at any position whose base is 'G' (acceptor track all-zero),
    so an inserted G-run manifests as a splice GAIN and a deleted one as a LOSS."""

    def _tracks(self, seq):  # type: ignore[override]
        donor = np.array([1.0 if b == "G" else 0.0 for b in seq], dtype=np.float32)
        acceptor = np.zeros(len(seq), dtype=np.float32)
        return acceptor, donor


def test_splice_delta_equal_length_positionwise():
    s = _StubSpliceAI(context=2)
    ref = "ATATATAT"  # no G anywhere -> donor track all 0
    alt = "ATGTATAT"  # one G introduced -> a site appears at that position
    assert s.delta(ref, alt) == 1.0  # SNV: position-wise |Δ| = 1.0 at the new G


def test_splice_delta_insertion_detects_gain():
    """The recoded-TJPE case: alt is LONGER (cargo insert). A G-bearing insert must
    be caught as a splice GAIN, not silently pass and not raise on length mismatch."""
    s = _StubSpliceAI(context=2)
    ref = "ATATATAT"
    alt = "ATAT" + "GGG" + "ATAT"  # 3-nt insert carrying donor sites
    d = s.delta(ref, alt)  # must NOT raise despite len(ref)=8 != len(alt)=11
    assert d == 1.0


def test_splice_delta_deletion_detects_loss():
    s = _StubSpliceAI(context=2)
    ref = "ATAT" + "GGG" + "ATAT"  # ref has donor sites in the middle
    alt = "ATATATAT"  # the G-region is deleted -> LOSS
    assert s.delta(ref, alt) == 1.0


def test_splice_delta_length_mismatch_no_change_is_low():
    """A length-changing edit that introduces no splice signal stays ~0 (not a
    spurious high delta from mis-alignment)."""
    s = _StubSpliceAI(context=2)
    ref = "ATATATAT"
    alt = "ATAT" + "ATA" + "ATAT"  # insert with no G -> no site
    assert s.delta(ref, alt) == 0.0


# --------------------------------------------------------------- Cas-OFFinder
class _StubCasOffinder(CasOffinderScorer):
    """CasOffinderScorer whose _run writes a canned cas-offinder output instead of
    invoking the binary, so aggregation is testable offline."""

    def __init__(self, canned_rows, **kw):
        super().__init__(**kw)
        self._canned = canned_rows  # list of (crRNA, chrom, pos, seq, strand, mism)

    def _run(self, input_path, output_path):  # type: ignore[override]
        with open(output_path, "w") as fh:
            for row in self._canned:
                fh.write("\t".join(str(x) for x in row) + "\n")


def test_casoffinder_unique_guide_is_clean():
    # One exact hit only -> ontarget 1.0, 0 off-targets, has_off False.
    g = "ACGTACGTACGTACGTACGT"  # 20 nt
    rows = [(g + "NGG", "chr1", 100, g + "TGG", "+", 0)]
    sc = _StubCasOffinder(rows, genome="/x.fa")
    res = sc([g])
    assert res[g] == (1.0, 0, False)


def test_casoffinder_counts_offtargets():
    g = "ACGTACGTACGTACGTACGT"
    rows = [
        (g + "NGG", "chr1", 100, g + "TGG", "+", 0),  # intended on-target
        (g + "NGG", "chr5", 900, g + "AGG", "+", 2),  # a 2-mismatch (near-tier) off-target
        (g + "NGG", "chr9", 55, g + "CGG", "-", 3),  # a 3-mismatch (distant) hit -> NOT counted
    ]
    sc = _StubCasOffinder(rows, genome="/x.fa")
    on, off, has = sc([g])[g]
    # Only the 2-mismatch near-tier hit disqualifies; the distant 3-mismatch hit is
    # searched but not counted (default count_mismatch_threshold=2).
    assert on == 1.0 and off == 1 and has is True


def test_casoffinder_guide_absent_is_flagged():
    g = "ACGTACGTACGTACGTACGT"
    sc = _StubCasOffinder([], genome="/x.fa")  # no hits at all
    on, off, has = sc([g])[g]
    assert on == 0.0  # not found exactly -> ontarget proxy 0 -> screen drops it


def test_casoffinder_empty_input():
    sc = _StubCasOffinder([], genome="/x.fa")
    assert sc([]) == {}


# --------------------------------------------------------------- structure fallback
def test_auto_structure_falls_back_to_rnaformer_on_importerror():
    scorer = AutoStructureScorer()

    class _RiNALMoMissing:
        def pair_probs(self, seq):
            raise ImportError("no flash_attn / no GPU")

    class _RNAformerOK:
        def pair_probs(self, seq):
            return {(0, 5): 0.9}

    scorer._rinalmo = _RiNALMoMissing()
    scorer._rnaformer = _RNAformerOK()
    out = scorer.pair_probs("ACGUACGU")
    assert out == {(0, 5): 0.9}
    assert scorer.backend_name == "rnaformer"


def test_auto_structure_prefers_rinalmo_when_available():
    scorer = AutoStructureScorer()

    class _RiNALMoOK:
        def pair_probs(self, seq):
            return {(1, 9): 0.7}

    scorer._rinalmo = _RiNALMoOK()
    out = scorer.pair_probs("ACGUACGUAC")
    assert out == {(1, 9): 0.7}
    assert scorer.backend_name == "rinalmo"
