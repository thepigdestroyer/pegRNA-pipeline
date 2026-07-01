"""Real-mode scorer-resource wiring for ``design_recoded`` (ADR 0008).

Regression tests for CX3 job 3158373: the 5/5 gate self-check passed but the
design aborted with ``real off-target scoring requires genome=`` because
``design_recoded`` never threaded the real-scorer resource paths into
``make_scorers``. These tests pin the contract:

- In ``real`` mode, ``design_recoded`` resolves the off-target genome FASTA, the
  Boltz cache, and the RiNALMo/RNAformer weights from an explicit parameter OR
  the same environment variable the smoke harness / ``run_recoded_pipeline.pbs``
  set (``PEGRNA_GENOME`` / ``BOLTZ_CACHE`` / ``RINALMO_SS_WEIGHTS`` /
  ``RNAFORMER_CKPT``), with the explicit parameter winning.
- Fake mode passes no scorer resources.
- Wiring never fabricates a path: with neither param nor env, ``genome`` stays
  ``None`` so the off-target scorer still fails loud.

The heavy real models are never imported: ``make_scorers`` is spied and returns
the deterministic fake bundle so the rest of ``design_recoded`` runs offline.
"""

from pegrna_pipeline.hpc import recoded_driver, scorer_factory

import pytest

_RESOURCE_ENV = ("PEGRNA_GENOME", "BOLTZ_CACHE", "RINALMO_SS_WEIGHTS", "RNAFORMER_CKPT")


def _demo_inputs():
    # ``genome`` here is the LOCAL ~240 bp design-context sequence consumed by
    # run_recoded_tjpe -- NOT a FASTA path. The off-target FASTA is genome_fasta.
    return dict(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome="ACGATGCATGGCATGACTGA" * 16,
        site_s=150,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
    )


@pytest.fixture
def spy(monkeypatch):
    """Record the kwargs make_scorers is called with; return the fake bundle so
    the rest of design_recoded runs deterministically without heavy models."""
    captured = {}

    def _spy(*, mode, **kwargs):
        captured["mode"] = mode
        captured["kwargs"] = kwargs
        return scorer_factory._fake_bundle()

    monkeypatch.setattr(recoded_driver, "make_scorers", _spy)
    for var in _RESOURCE_ENV:
        monkeypatch.delenv(var, raising=False)
    return captured


def test_fake_mode_passes_no_scorer_resources(spy):
    recoded_driver.design_recoded(**_demo_inputs(), mode="fake")
    assert spy["mode"] == "fake"
    assert spy["kwargs"] == {}


def test_real_mode_threads_explicit_genome_fasta(spy):
    recoded_driver.design_recoded(
        **_demo_inputs(), mode="real", genome_fasta="/staged/GRCh38.fa", uniprot_acc="P49711"
    )
    assert spy["mode"] == "real"
    assert spy["kwargs"]["genome"] == "/staged/GRCh38.fa"
    assert spy["kwargs"]["uniprot_acc"] == "P49711"


def test_real_mode_reads_pegrna_genome_env(spy, monkeypatch):
    monkeypatch.setenv("PEGRNA_GENOME", "/env/GRCh38.fa")
    recoded_driver.design_recoded(**_demo_inputs(), mode="real")
    assert spy["kwargs"]["genome"] == "/env/GRCh38.fa"


def test_real_mode_explicit_genome_beats_env(spy, monkeypatch):
    monkeypatch.setenv("PEGRNA_GENOME", "/env/GRCh38.fa")
    recoded_driver.design_recoded(**_demo_inputs(), mode="real", genome_fasta="/explicit/GRCh38.fa")
    assert spy["kwargs"]["genome"] == "/explicit/GRCh38.fa"


def test_real_mode_no_genome_leaves_none(spy):
    """Fail-loud preserved: the wiring must not invent a genome path."""
    recoded_driver.design_recoded(**_demo_inputs(), mode="real")
    assert spy["kwargs"]["genome"] is None


def test_real_mode_threads_boltz_and_structure_resources_from_env(spy, monkeypatch):
    monkeypatch.setenv("BOLTZ_CACHE", "/data/boltz_cache")
    monkeypatch.setenv("RINALMO_SS_WEIGHTS", "/data/rinalmo/w.pt")
    monkeypatch.setenv(
        "RNAFORMER_CKPT", "/data/rnaformer_weights/RNAformer_32M_state_dict_biophysical.pth"
    )
    recoded_driver.design_recoded(**_demo_inputs(), mode="real")
    kw = spy["kwargs"]
    assert kw["boltz_cache"] == "/data/boltz_cache"
    assert kw["rinalmo_weights"] == "/data/rinalmo/w.pt"
    # RNAFORMER_CKPT resolves to its containing directory (mirrors smoke_models.py).
    assert kw["rnaformer_model_dir"] == "/data/rnaformer_weights"


def test_real_mode_explicit_resource_params_beat_env(spy, monkeypatch):
    monkeypatch.setenv("BOLTZ_CACHE", "/env/boltz")
    monkeypatch.setenv("RINALMO_SS_WEIGHTS", "/env/rinalmo.pt")
    monkeypatch.setenv("RNAFORMER_CKPT", "/env/rnaformer/ckpt.pth")
    recoded_driver.design_recoded(
        **_demo_inputs(),
        mode="real",
        boltz_cache="/explicit/boltz",
        rinalmo_weights="/explicit/rinalmo.pt",
        rnaformer_model_dir="/explicit/rnaformer",
    )
    kw = spy["kwargs"]
    assert kw["boltz_cache"] == "/explicit/boltz"
    assert kw["rinalmo_weights"] == "/explicit/rinalmo.pt"
    assert kw["rnaformer_model_dir"] == "/explicit/rnaformer"


def test_real_mode_threads_pam_and_max_mismatches_defaults(spy):
    recoded_driver.design_recoded(**_demo_inputs(), mode="real")
    kw = spy["kwargs"]
    assert kw["pam"] == "NGG"
    assert kw["max_mismatches"] == 3


def test_real_mode_threads_custom_pam_and_max_mismatches(spy):
    recoded_driver.design_recoded(**_demo_inputs(), mode="real", pam="NNGRRT", max_mismatches=4)
    kw = spy["kwargs"]
    assert kw["pam"] == "NNGRRT"
    assert kw["max_mismatches"] == 4


def test_real_mode_bare_rnaformer_ckpt_resolves_to_none(spy, monkeypatch):
    # A ckpt env value with no directory component must degrade to None, not "".
    monkeypatch.setenv("RNAFORMER_CKPT", "ckpt.pth")
    recoded_driver.design_recoded(**_demo_inputs(), mode="real")
    assert spy["kwargs"]["rnaformer_model_dir"] is None


def test_real_mode_empty_env_normalized_to_none(spy, monkeypatch):
    # Hand-exported empty env vars must not reach a scorer as an empty-string path.
    monkeypatch.setenv("PEGRNA_GENOME", "")
    monkeypatch.setenv("BOLTZ_CACHE", "")
    monkeypatch.setenv("RINALMO_SS_WEIGHTS", "")
    recoded_driver.design_recoded(**_demo_inputs(), mode="real")
    kw = spy["kwargs"]
    assert kw["genome"] is None
    assert kw["boltz_cache"] is None
    assert kw["rinalmo_weights"] is None
