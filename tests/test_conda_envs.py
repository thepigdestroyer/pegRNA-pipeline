"""Verify bundled conda environment YAML files exist and are well-formed."""

from pathlib import Path

import yaml

ENVS_DIR = Path(__file__).resolve().parent.parent / "envs"
REQUIRED_ENVS = ("agent.yml", "pipeline.yml", "offtarget.yml")


def test_all_required_env_files_exist():
    for name in REQUIRED_ENVS:
        assert (ENVS_DIR / name).is_file(), f"missing envs/{name}"


def test_env_files_are_valid_yaml_with_name_and_deps():
    for name in REQUIRED_ENVS:
        data = yaml.safe_load((ENVS_DIR / name).read_text())
        assert isinstance(data, dict), f"{name} is not a mapping"
        assert "name" in data, f"{name} missing 'name' field"
        assert "dependencies" in data, f"{name} missing 'dependencies' field"
        assert isinstance(data["dependencies"], list) and data["dependencies"], (
            f"{name} has empty dependencies"
        )
