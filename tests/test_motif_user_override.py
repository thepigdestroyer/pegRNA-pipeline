"""Tests for user motif overrides at ~/.pegrna/motifs.yaml."""

import pytest
import yaml

from pegrna_pipeline.agent.motif_library import MotifLibrary


@pytest.fixture
def tmp_home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    return tmp_path


def test_load_user_default_returns_empty_when_missing(tmp_home):
    lib = MotifLibrary.load_user_default()
    assert isinstance(lib, MotifLibrary)
    assert len(lib) == 0


def test_load_user_default_reads_home_motifs_yaml(tmp_home):
    user_dir = tmp_home / ".pegrna"
    user_dir.mkdir()
    user_yaml = user_dir / "motifs.yaml"
    user_yaml.write_text(
        yaml.safe_dump(
            {
                "motifs": [
                    {
                        "name": "MyTag",
                        "category": "custom",
                        "dna": "ATGAAG",
                        "aa": "MK",
                        "description": "lab-specific tag",
                    }
                ]
            }
        )
    )
    lib = MotifLibrary.load_user_default()
    assert "MyTag" in {m.name for m in lib}
    assert lib.lookup("mytag").aa == "MK"


def test_load_full_merges_user_over_bundled(tmp_home):
    user_dir = tmp_home / ".pegrna"
    user_dir.mkdir()
    (user_dir / "motifs.yaml").write_text(
        yaml.safe_dump(
            {
                "motifs": [
                    # Override bundled FLAG with different DNA.
                    {
                        "name": "FLAG",
                        "category": "epitope_tag",
                        "dna": "GATTATAAGGATGATGATGATAAG",
                        "aa": "DYKDDDDK",
                        "description": "lab flavor FLAG",
                    },
                    {
                        "name": "LocalTag",
                        "category": "custom",
                        "dna": "ATG",
                        "aa": "M",
                        "description": "",
                    },
                ]
            }
        )
    )
    lib = MotifLibrary.load_full()
    # User FLAG takes precedence.
    flag = lib.lookup("FLAG")
    assert flag.dna == "GATTATAAGGATGATGATGATAAG"
    # LocalTag is available.
    assert "LocalTag" in {m.name for m in lib}
    # Bundled motifs still resolvable.
    assert lib.lookup("GGGGS").aa == "GGGGS"


def test_load_full_works_without_user_file(tmp_home):
    lib = MotifLibrary.load_full()
    # Bundled FLAG is pristine.
    assert lib.lookup("FLAG").dna == "GACTACAAGGACGACGATGACAAG"
