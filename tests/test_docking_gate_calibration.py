"""The docking gate, pinned to the Boltz-2 measurements it was calibrated on.

Every number below was measured on the live NVIDIA Boltz-2 NIM (SpCas9 + RNA, 5 diffusion
samples per call). They are fixtures, not guesses: the point of this file is that the gate
accepts what was measured to be a normal pegRNA and rejects what was measured to be a
non-binder, at the lengths where that is a non-trivial thing to do.

Background. The gate used an absolute floor of interface ipTM >= 0.6. But ipTM falls with
RNA length regardless of design quality, so the floor was unreachable for tjPE-sized
pegRNAs: neutral extensions of 194/232 nt average 0.48/0.43, and the designed 194-nt
pegRNA scored 0.544. Meanwhile the canonical SpCas9:sgRNA complex scored 0.69 and a
scrambled non-binder reached 0.634 in a single draw, so single draws were neither
reachable-for-valid nor safe-against-invalid.
"""

import pytest

from pegrna_pipeline.hpc.complex_dock import (
    DockFake,
    DockingServiceError,
    docking_reference,
    passes_calibrated_docking_gate,
    passes_docking_gate,
)
from pegrna_pipeline.hpc.complex_dock_nim import Boltz2NIMDockScorer


class _Samples:
    """Returns fixed per-sample ipTMs, as the NIM scorer does."""

    def __init__(self, samples):
        self._s = samples

    def interface_iptm_samples(self, pegrna, cas9):  # noqa: ARG002
        return self._s


def _rna(n):
    return "A" * n


def gate(length, samples):
    return passes_calibrated_docking_gate(_rna(length), "nCas9", scorer=_Samples(samples))


# ---- designs and controls that were measured -------------------------------------------

# 5 samples each, total RNA length in nt.
DESIGNED_194 = (194, [0.518, 0.541, 0.628, 0.502, 0.529])  # the pipeline's own design
DESIGNED_232 = (232, [0.451, 0.409, 0.400, 0.439, 0.459])
NEUTRAL_RANDOM_194 = [
    [0.452, 0.439, 0.507, 0.518, 0.483],
    [0.523, 0.618, 0.482, 0.631, 0.548],
    [0.394, 0.467, 0.433, 0.504, 0.455],
]
NEUTRAL_RANDOM_232 = [
    [0.470, 0.496, 0.485, 0.489, 0.433],
    [0.385, 0.413, 0.429, 0.428, 0.383],
    [0.477, 0.382, 0.370, 0.339, 0.336],
]
NEUTRAL_RANDOM_156 = [
    [0.538, 0.537, 0.649, 0.602, 0.614],
    [0.645, 0.688, 0.633, 0.634, 0.602],
    [0.549, 0.649, 0.643, 0.683, 0.560],
]
SEQUESTERING_156 = (156, [0.324, 0.315, 0.303, 0.284, 0.398])  # positive control
SEQUESTERING_194 = (194, [0.279, 0.286, 0.250, 0.286, 0.276])
SCRAMBLED_PLUS_EXT_194 = [[0.265, 0.264, 0.245, 0.249, 0.243], [0.252, 0.346, 0.257, 0.302, 0.296]]
SCRAMBLED_PLUS_EXT_232 = [[0.272, 0.240, 0.223, 0.302, 0.239], [0.301, 0.239, 0.274, 0.226, 0.303]]
SCRAMBLED_96 = [0.477, 0.331, 0.448, 0.347, 0.634]  # includes a single draw above 0.6
CANONICAL_96 = [0.763, 0.676, 0.752, 0.625, 0.648]


def test_the_designed_pegrnas_that_the_old_floor_rejected_now_pass():
    for length, samples in (DESIGNED_194, DESIGNED_232):
        ok, why = gate(length, samples)
        assert ok, why


def test_the_old_floor_really_did_reject_them():
    """Guard the premise: under the absolute floor these same designs fail."""
    for _, samples in (DESIGNED_194, DESIGNED_232):
        assert sum(samples) / len(samples) < 0.6


@pytest.mark.parametrize("samples", NEUTRAL_RANDOM_194)
def test_neutral_extensions_pass_at_194(samples):
    ok, why = gate(194, samples)
    assert ok, why


@pytest.mark.parametrize("samples", NEUTRAL_RANDOM_232)
def test_neutral_extensions_pass_at_232(samples):
    ok, why = gate(232, samples)
    assert ok, why


@pytest.mark.parametrize("samples", NEUTRAL_RANDOM_156)
def test_neutral_extensions_pass_at_156(samples):
    ok, why = gate(156, samples)
    assert ok, why


def test_the_scaffold_sequestering_positive_control_is_rejected():
    """The failure mode this gate exists for: an extension that ties up the scaffold.
    The folding gate ignores every pair touching the scaffold, so only this catches it."""
    for length, samples in (SEQUESTERING_156, SEQUESTERING_194):
        ok, why = gate(length, samples)
        assert not ok, why
        assert "non-binder" in why


@pytest.mark.parametrize("samples", SCRAMBLED_PLUS_EXT_194)
def test_non_binders_are_rejected_at_194(samples):
    ok, why = gate(194, samples)
    assert not ok, why


@pytest.mark.parametrize("samples", SCRAMBLED_PLUS_EXT_232)
def test_non_binders_are_rejected_at_232(samples):
    ok, why = gate(232, samples)
    assert not ok, why


def test_one_lucky_draw_does_not_rescue_a_non_binder():
    """The scrambled RNA drew 0.634 once, above the old 0.6 floor. Its mean is 0.448."""
    assert max(SCRAMBLED_96) > 0.6
    ok, _ = gate(96, SCRAMBLED_96)
    assert not ok


def test_the_canonical_functional_complex_passes():
    ok, why = gate(96, CANONICAL_96)
    assert ok, why


# ---- reference behaviour --------------------------------------------------------------


def test_reference_decreases_with_length():
    lengths = [96, 126, 156, 194, 232]
    neutral = [docking_reference(n).neutral for n in lengths]
    assert neutral == sorted(neutral, reverse=True) or neutral[1] <= neutral[2]
    assert docking_reference(232).neutral < docking_reference(96).neutral


def test_non_binder_is_always_below_neutral():
    for n in range(96, 241, 4):
        r = docking_reference(n)
        assert r.nonbinder < r.neutral and r.span > 0.1


def test_the_cutoff_is_a_fraction_of_the_way_up():
    r = docking_reference(194)
    assert r.cutoff(0.0) == pytest.approx(r.nonbinder)
    assert r.cutoff(1.0) == pytest.approx(r.neutral)
    assert r.cutoff(0.5) == pytest.approx((r.neutral + r.nonbinder) / 2)


def test_interpolation_is_continuous_across_table_rows():
    a, b = docking_reference(193), docking_reference(195)
    assert abs(a.neutral - b.neutral) < 0.01


# ---- outside the calibrated range: report, do not judge --------------------------------


def test_a_very_long_pegrna_is_advisory_not_rejected_and_not_silently_accepted():
    metrics = {}
    ok, why = passes_calibrated_docking_gate(
        _rna(673), "nCas9", scorer=_Samples([0.2, 0.2, 0.2, 0.2, 0.2]), metrics=metrics
    )
    assert ok is True
    assert why.startswith("ADVISORY ONLY")
    assert metrics["advisory"] is True and "note" in metrics


def test_the_calibrated_range_boundary():
    from pegrna_pipeline.hpc.gate_thresholds import DOCKING_CALIBRATED_RANGE

    lo, hi = DOCKING_CALIBRATED_RANGE
    assert docking_reference(hi).calibrated and not docking_reference(hi + 1).calibrated
    assert docking_reference(lo).calibrated and not docking_reference(lo - 1).calibrated


def test_metrics_report_value_cutoff_and_both_references():
    metrics = {}
    passes_calibrated_docking_gate(
        _rna(194), "nCas9", scorer=_Samples(DESIGNED_194[1]), metrics=metrics
    )
    assert metrics["cmp"] == ">="
    assert metrics["value"] == pytest.approx(0.544, abs=1e-3)
    assert metrics["cutoff"] == pytest.approx(docking_reference(194).cutoff(0.5), abs=1e-3)
    assert "neutral" in metrics["extra"] and "non-binder" in metrics["extra"]


def test_no_samples_is_a_service_error_not_a_verdict():
    with pytest.raises(DockingServiceError):
        gate(194, [])


# ---- dispatch ---------------------------------------------------------------------------


def test_sampling_scorers_are_routed_to_the_calibrated_gate():
    ok, why = passes_docking_gate(_rna(194), "nCas9", scorer=_Samples(DESIGNED_194[1]))
    assert ok and "neutral extension" in why


def test_scorers_without_samples_keep_the_legacy_gate():
    """DockFake has no per-sample API, so existing behaviour is unchanged."""
    ok, why = passes_docking_gate(
        _rna(194), "nCas9", scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.83)
    )
    assert ok and "interface preserved" in why
    ok, why = passes_docking_gate(
        _rna(194), "nCas9", scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.55)
    )
    assert not ok and "floor" in why


# ---- NIM response handling --------------------------------------------------------------


def _pcs(*values):
    return [{"0": {"0": 0.0, "1": v}, "1": {"0": v, "1": 0.0}} for v in values]


def test_nim_extracts_every_samples_interface_value():
    s = Boltz2NIMDockScorer(api_key="x")
    vals = s._extract_samples({"pair_chains_iptm_scores": _pcs(0.5, 0.6, 0.7, 0.55, 0.45)})
    assert vals == [0.5, 0.6, 0.7, 0.55, 0.45]
    assert s.last_metric_source == "pair_chains_iptm_scores"


def test_nim_single_number_view_is_the_first_sample():
    s = Boltz2NIMDockScorer(api_key="x")
    assert s._extract_iptm({"pair_chains_iptm_scores": _pcs(0.5, 0.9)}) == 0.5


def test_calibrated_gate_refuses_a_metric_it_was_not_calibrated_on(monkeypatch):
    """Falling back to the global ipTM or an aggregate proxy would compare the wrong
    quantity against the references, so it must fail loudly instead."""
    s = Boltz2NIMDockScorer(api_key="x")
    monkeypatch.setattr(s, "_resolve_cas9", lambda c: "M" * 50)
    monkeypatch.setattr(
        s,
        "_predict_samples",
        lambda rna, prot: (setattr(s, "last_metric_source", "confidence_scores(proxy)"), [0.5])[1],
    )
    with pytest.raises(DockingServiceError) as exc:
        s.interface_iptm_samples("ACGU" * 20, "nCas9")
    assert "calibration" in str(exc.value)


def test_nim_requests_the_configured_number_of_samples(monkeypatch):
    from pegrna_pipeline.hpc import complex_dock_nim

    seen = {}

    class _R:
        def read(self):
            import json

            return json.dumps({"pair_chains_iptm_scores": _pcs(0.5, 0.5, 0.5)}).encode()

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        import json

        seen["payload"] = json.loads(req.data)
        return _R()

    monkeypatch.setattr(complex_dock_nim.urllib.request, "urlopen", fake_urlopen)
    Boltz2NIMDockScorer(api_key="x", samples=3)._predict_samples("ACGU" * 10, "M" * 30)
    assert seen["payload"]["diffusion_samples"] == 3


# ---- borderline verdicts ----------------------------------------------------------------
#
# In the calibration data the closest calls on either side (the lowest neutral design and the
# highest non-binder) both sit within 15 points of the cutoff, so a verdict there is not robust.
# The flag reports that; it must never change what the rule decides.


def _mean_at(length: int, f: float) -> float:
    """The mean ipTM that sits ``f`` of the way from non-binder (0) to neutral (1)."""
    ref = docking_reference(length)
    return ref.nonbinder + f * ref.span


def _gate_at(length: int, f: float, **kw):
    metrics = {}
    ok, why = passes_calibrated_docking_gate(
        _rna(length), "nCas9", scorer=_Samples([_mean_at(length, f)] * 5), metrics=metrics, **kw
    )
    return ok, why, metrics


@pytest.mark.parametrize("length", [126, 156, 194, 232])
def test_a_verdict_just_inside_the_band_is_flagged_either_way(length):
    ok_lo, why_lo, m_lo = _gate_at(length, 0.42)
    ok_hi, why_hi, m_hi = _gate_at(length, 0.58)
    assert (ok_lo, ok_hi) == (False, True), "the flag must not move the verdict"
    assert m_lo["borderline"] and m_hi["borderline"]
    assert "BORDERLINE" in why_lo and "BORDERLINE" in why_hi
    assert "not robust" in m_lo["note"]


@pytest.mark.parametrize("f", [0.0, 0.2, 0.3, 0.7, 0.9, 1.1])
def test_a_clear_verdict_is_not_flagged(f):
    ok, why, metrics = _gate_at(194, f)
    assert ok is (f >= 0.5)
    assert "borderline" not in metrics and "BORDERLINE" not in why


def test_the_band_is_fifteen_points_either_side_of_the_cutoff():
    assert _gate_at(194, 0.5 - 0.14)[2].get("borderline")
    assert _gate_at(194, 0.5 + 0.14)[2].get("borderline")
    assert not _gate_at(194, 0.5 - 0.16)[2].get("borderline")
    assert not _gate_at(194, 0.5 + 0.16)[2].get("borderline")


def test_the_flag_never_changes_the_verdict():
    for i in range(0, 121, 3):
        f = i / 100
        flagged, _, _ = _gate_at(194, f)
        unflagged, _, _ = _gate_at(194, f, borderline_margin=0.0)
        assert flagged == unflagged == (f >= 0.5 - 1e-9), f


def test_the_margin_is_configurable():
    assert not _gate_at(194, 0.5, borderline_margin=0.0)[2].get("borderline")
    assert _gate_at(194, 0.05, borderline_margin=1.0)[2].get("borderline")


def test_the_advisory_range_is_never_called_borderline():
    ok, why, metrics = _gate_at(300, 0.5)
    assert ok and why.startswith("ADVISORY ONLY") and "borderline" not in metrics


def test_the_margin_comes_from_the_documented_threshold():
    from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS

    assert THRESHOLDS.docking_borderline_margin == 0.15
    assert "docking_borderline_margin" in THRESHOLDS.citations


def test_a_borderline_pass_reaches_the_designs_warnings():
    """A reader ordering DNA needs to see it without opening the full report."""
    import test_pridict_ngrna_wiring as W

    from pegrna_pipeline.hpc.recoded_tjpe_workflow import run_recoded_tjpe

    class OnTheLine:
        def interface_iptm_samples(self, pegrna, cas9):  # noqa: ARG002
            return [_mean_at(len(pegrna), 0.55)] * 5

    res = run_recoded_tjpe(**W._wf_inputs(dock_scorer=OnTheLine()))
    assert any(w.startswith("docking gate: borderline") for w in res.warnings), res.warnings


def test_a_clear_pass_adds_no_docking_warning():
    import test_pridict_ngrna_wiring as W

    from pegrna_pipeline.hpc.recoded_tjpe_workflow import run_recoded_tjpe

    class Comfortable:
        def interface_iptm_samples(self, pegrna, cas9):  # noqa: ARG002
            return [_mean_at(len(pegrna), 0.95)] * 5

    res = run_recoded_tjpe(**W._wf_inputs(dock_scorer=Comfortable()))
    assert not any("docking" in w for w in res.warnings)
