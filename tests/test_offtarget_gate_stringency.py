"""Off-target gate stringency + 0-nick diagnosability (root-cause fix for job 3158586).

The real-mode CTCF run ended in `NoValidPegRNA after 0 tried` with EMPTY diagnostics:
`screen_nick_candidates` dropped every nick because `CasOffinderScorer` counted all
<=3-mismatch genomic hits (17-448 per guide vs 3.1 Gb GRCh38) as off-targets, and the
gate requires 0. A 3-mismatch hit is not a functional off-target. The scorer now counts
only near-tier hits (mismatches <= count_mismatch_threshold, default 2) plus extra exact
duplicates; distant hits are still searched but not disqualifying. And run_recoded_tjpe
now emits an explanatory diagnostic when the nick pool is empty.
"""

from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer
from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA, run_recoded_tjpe
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics

G = "ACGTACGTACGTACGTACGT"


class _StubCasOffinder(CasOffinderScorer):
    """CasOffinderScorer whose _run writes canned cas-offinder rows (offline)."""

    def __init__(self, canned_rows, **kw):
        super().__init__(**kw)
        self._canned = canned_rows

    def _run(self, input_path, output_path):  # type: ignore[override]
        with open(output_path, "w") as fh:
            for row in self._canned:
                fh.write("\t".join(str(x) for x in row) + "\n")


def _rows(guide, tiers):
    """tiers: {mismatch_count: n_hits} -> canned cas-offinder rows."""
    out = []
    for mm, n in tiers.items():
        for _ in range(n):
            out.append((guide + "NGG", "chr1", 1, guide + "TGG", "+", mm))
    return out


# --- scorer counting ---------------------------------------------------------


def test_default_count_threshold_is_two():
    assert CasOffinderScorer(genome="/x.fa").count_mismatch_threshold == 2


def test_distant_3mm_hits_are_not_offtargets():
    # Unique exact hit + 17 distant 3mm hits -> a clean guide (mm3 ignored).
    sc = _StubCasOffinder(_rows(G, {0: 1, 3: 17}), genome="/x.fa")
    assert sc([G])[G] == (1.0, 0, False)


def test_near_tier_hits_are_counted():
    # 2 hits at 1mm + 5 at 2mm = 7 near-tier off-targets; 30 at 3mm ignored.
    sc = _StubCasOffinder(_rows(G, {0: 1, 1: 2, 2: 5, 3: 30}), genome="/x.fa")
    assert sc([G])[G] == (1.0, 7, True)


def test_exact_duplicate_counts_as_offtarget():
    # Two exact (0mm) hits: one on-target + one exact-duplicate off-target.
    sc = _StubCasOffinder(_rows(G, {0: 2}), genome="/x.fa")
    assert sc([G])[G] == (1.0, 1, True)


def test_count_threshold_is_configurable():
    # With threshold 3, the 3mm hits DO count.
    sc = _StubCasOffinder(_rows(G, {0: 1, 3: 17}), genome="/x.fa", count_mismatch_threshold=3)
    assert sc([G])[G] == (1.0, 17, True)


# --- 0-nick diagnosability ---------------------------------------------------


def _importance():
    class _S:
        def metrics(self, i):
            return ResidueMetrics(
                plddt=50.0,
                rsa=0.5,
                contact_number=4,
                is_catalytic=False,
                is_disulfide=False,
                is_positive_phi_gly=False,
                is_structural_pro=False,
            )

    return _S()


def _inputs(**over):
    base = dict(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome="ACGATGCATGGCATGACTGA" * 16,
        site_s=150,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        ngrna_distance_min=40,
        ngrna_distance_max=70,
        pegrna_strand="Rv",
        max_run=8,
    )
    base.update(over)
    return base


def test_empty_nick_pool_gives_explanatory_diagnostic():
    def reject_all(protospacers):
        return {p: (0.0, 999, True) for p in protospacers}

    try:
        run_recoded_tjpe(
            **_inputs(),
            splice_scorer=SpliceFake(delta=0.0),
            full_structure_scorer=FullStructureFake(pairs={}),
            importance_scorer=_importance(),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
            offtarget_scorer=reject_all,
        )
        raise AssertionError("expected NoValidPegRNA")
    except NoValidPegRNA as exc:
        msg = str(exc).lower()
        assert "0 nicks after off-target screen" in msg
        # the diagnostic must not be blank (the old bug)
        assert not msg.rstrip().endswith("diagnostics:")
