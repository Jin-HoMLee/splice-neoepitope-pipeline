from pathlib import Path

import pytest

from app.config_contract import ConfigValidationError, LocalRunConfig


def make_config(tmp_path: Path, **overrides):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("patient_id\tsample_type\tfastq1\nP1\tPrimary Tumor\ttumor.fastq.gz\n")
    values = {"manifest": manifest, "output_root": tmp_path / "run"}
    values.update(overrides)
    return LocalRunConfig(**values)


def test_default_is_ten_reads_hisat2_and_safe(tmp_path):
    config = make_config(tmp_path).validate()
    overlay = config.overlay()
    assert overlay["filtering"]["min_junction_reads"] == 10
    assert overlay["alignment"] == {"aligner": "hisat2", "uniqueness_filter": {"enabled": False}}
    assert overlay["tcrdock"]["enabled"] is False


def test_five_reads_requires_acknowledgement(tmp_path):
    with pytest.raises(ConfigValidationError, match="requires acknowledgement"):
        make_config(tmp_path, depth=5).validate()
    make_config(tmp_path, depth=5, sensitivity_acknowledged=True).validate()


@pytest.mark.parametrize("depth", [0, 4, 6, 11])
def test_unlisted_depth_is_rejected(tmp_path, depth):
    with pytest.raises(ConfigValidationError, match="depth must be"):
        make_config(tmp_path, depth=depth).validate()


def test_remote_and_unsafe_controls_are_rejected(tmp_path):
    with pytest.raises(ConfigValidationError, match="local path"):
        LocalRunConfig(Path("gs://bucket/samples.tsv"), tmp_path / "run").validate(require_paths=False)
    with pytest.raises(ConfigValidationError, match="locked off"):
        make_config(tmp_path, uniqueness_filter=True).validate()
    with pytest.raises(ConfigValidationError, match="outside"):
        make_config(tmp_path, tcrdock=True).validate()


def test_star_is_deferred_from_phase_zero(tmp_path):
    with pytest.raises(ConfigValidationError, match="HISAT2 only"):
        make_config(tmp_path, aligner="star").validate()


def test_manifest_requires_one_patient_and_primary_tumor(tmp_path):
    config = make_config(tmp_path)
    config.manifest.write_text(
        "patient_id\tsample_type\tfastq1\nP1\tSolid Tissue Normal\tn.fastq.gz\n"
    )
    with pytest.raises(ConfigValidationError, match="Primary Tumor"):
        config.validate()
    config.manifest.write_text(
        "patient_id\tsample_type\tfastq1\nP1\tPrimary Tumor\ta.fastq.gz\nP2\tPrimary Tumor\tb.fastq.gz\n"
    )
    with pytest.raises(ConfigValidationError, match="exactly one patient"):
        config.validate()
