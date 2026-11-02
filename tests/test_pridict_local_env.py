"""PRIDICT2's bundled DeepCas9 needs Keras 2, and says so very badly when it
doesn't get it.

`trained_models/DeepCas9_TestCode.py` is TF1 graph-mode code that passes a
symbolic `tf.bool` placeholder as `training=` to `tf.keras.layers.Dropout`. Under
Keras 3 (the default `tf.keras` from TensorFlow >= 2.16) that raises "Using a
symbolic `tf.Tensor` as a Python `bool` is not allowed". PRIDICT catches it per
row, logs it, and exits 0 having written no pegRNA CSV -- which the runner then
reports as "the edit is likely outside its design domain (insertions must be
<= 40 bp)". The real cause is nowhere in that message, so the runner requests the
legacy-Keras shim itself rather than leaving it to a deployment note.
"""

import os

from pegrna_pipeline.hpc.pridict_local import PridictLocalRunner


def _runner(tmp_path):
    return PridictLocalRunner(
        pridict_dir=str(tmp_path),
        pridict_python="/usr/bin/python3",
    )


def test_requests_legacy_keras_for_the_subprocess(tmp_path, monkeypatch):
    monkeypatch.delenv("TF_USE_LEGACY_KERAS", raising=False)
    env = _runner(tmp_path)._subprocess_env()
    assert env["TF_USE_LEGACY_KERAS"] == "1"


def test_does_not_override_an_explicit_caller_setting(tmp_path, monkeypatch):
    """An operator who has deliberately pinned this keeps control."""
    monkeypatch.setenv("TF_USE_LEGACY_KERAS", "0")
    env = _runner(tmp_path)._subprocess_env()
    assert env["TF_USE_LEGACY_KERAS"] == "0"


def test_preserves_the_rest_of_the_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("PEGRNA_TEST_SENTINEL", "kept")
    env = _runner(tmp_path)._subprocess_env()
    assert env["PEGRNA_TEST_SENTINEL"] == "kept"
    assert "PATH" in env and env["PATH"] == os.environ["PATH"]
