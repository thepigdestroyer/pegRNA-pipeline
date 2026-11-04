"""The candidate walk evaluates each DISTINCT design once, and ``rtt_aa`` does not change a design.

Found while sampling the pipeline's real candidates for the docking calibration: the CTCF pool of
18 candidates (6 synonymous-codon "RTT variants" x 3 PBS lengths) is really 3 pegRNAs, because the
walk only ever uses a variant's PBS length. Five different anchor residues produce byte-identical
pegRNAs. The repeats were not free: a design that failed re-ran the docking gate for every
duplicate, and each run is a ~40 s paid Boltz-2 request. And the UI told users the pipeline
"re-encodes" the anchor residue with a synonymous codon, which it never did.
"""

import pytest
import test_pridict_ngrna_wiring as W

from pegrna_pipeline.hpc import recoded_driver
from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    enumerate_pbs_rtt_variants,
    run_recoded_tjpe,
)

PBS = (13, 15, 17)


class RejectingDock(DockFake):
    """Rejects every design and records each pegRNA it is shown."""

    def __init__(self):
        super().__init__(apo_iptm=0.85, candidate_iptm=0.10)
        self.shown = []

    def interface_iptm(self, pegrna, cas9):
        self.shown.append(pegrna)
        return super().interface_iptm(pegrna, cas9)


def _walk(**over):
    dock = RejectingDock()
    inputs = W._wf_inputs(dock_scorer=dock, pbs_lengths=PBS, rtt_aa="L")
    inputs.update(over)
    with pytest.raises(NoValidPegRNA) as exc:
        run_recoded_tjpe(**inputs)
    return dock, exc.value


def test_each_distinct_design_reaches_the_docking_gate_exactly_once():
    dock, _ = _walk()
    assert dock.shown, "the walk must reach docking for this to mean anything"
    assert len(dock.shown) == len(set(dock.shown)), "a design was docked more than once"


def test_the_number_of_docking_calls_is_designs_not_variants():
    """Leucine has 6 codons, so 6 RTT variants x 3 PBS lengths = 18 candidates per nick."""
    n_variants = len(enumerate_pbs_rtt_variants(spacer="A" * 20, rtt_aa="L", pbs_lengths=PBS))
    assert n_variants == 18
    dock, _ = _walk()
    nicks = len({p[96 + 0 : 96 + 17] for p in dock.shown})  # rcPBS2 identifies the nick
    assert len(dock.shown) <= nicks * len(PBS), "at most one design per (nick, PBS length)"
    assert len(dock.shown) < n_variants * nicks


def test_candidates_tried_counts_designs():
    dock, err = _walk()
    assert f"after {len(dock.shown)} tried" in str(err)


def test_a_residue_with_one_codon_gives_the_same_count_as_one_with_six():
    """Methionine has one codon, leucine six; the designs and the docking calls are identical."""
    a, _ = _walk(rtt_aa="M")
    b, _ = _walk(rtt_aa="L")
    assert a.shown == b.shown


# ---- rtt_aa does not change the design --------------------------------------------------


@pytest.mark.parametrize("aa", ["L", "A", "W", "M", "S", "MKP"])
def test_the_anchor_residue_does_not_change_the_pegrna(aa):
    base = run_recoded_tjpe(**W._wf_inputs(pbs_lengths=PBS, rtt_aa="L")).pegRNA_sequence
    other = run_recoded_tjpe(**W._wf_inputs(pbs_lengths=PBS, rtt_aa=aa)).pegRNA_sequence
    assert other == base


def test_the_anchor_residue_is_optional():
    inputs = W._wf_inputs(pbs_lengths=PBS)
    inputs.pop("rtt_aa")
    with_aa = run_recoded_tjpe(**W._wf_inputs(pbs_lengths=PBS, rtt_aa="L")).pegRNA_sequence
    assert run_recoded_tjpe(**inputs).pegRNA_sequence == with_aa


def test_without_an_anchor_there_is_one_variant_per_pbs_length():
    assert len(enumerate_pbs_rtt_variants(spacer="A" * 20, pbs_lengths=PBS)) == 3


def test_the_driver_accepts_a_config_without_it():
    """A blank field in the web form is dropped from the config; that must not be an error."""
    inputs = W._planted_inputs()
    inputs.pop("rtt_aa", None)
    runner = W._NickRunner(spacer="TGCATGCATGCATGCATGCA", nicking=[W._row(1)])
    out = recoded_driver.design_recoded(**inputs, mode="fake", pridict_runner=runner)
    assert out.result.pegRNA_sequence
