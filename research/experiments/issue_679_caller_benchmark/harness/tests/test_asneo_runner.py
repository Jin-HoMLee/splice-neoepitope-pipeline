"""Tests for the ASNEO runner (#966, AC-1).

Only the pure output-contract helpers are unit-tested here; the actual
subprocess run against the asneo conda env is the recorded live-integration
smoke (plan Task 6), not a unit test.
"""

import pytest

from runners.asneo import peptide_output_path, validate_output


def test_peptide_output_path_is_under_workdir(tmp_path):
    assert peptide_output_path(str(tmp_path)).endswith("putative_peptide.txt")


def test_validate_output_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_output(str(tmp_path / "nope.txt"))


def test_validate_output_zero_peptides_raises(tmp_path):
    f = tmp_path / "putative_peptide.txt"
    f.write_text("# only a comment, no peptides\n\n")
    with pytest.raises(ValueError):
        validate_output(str(f))


def test_validate_output_with_peptides_passes(tmp_path):
    f = tmp_path / "putative_peptide.txt"
    f.write_text("# header\nLEQGTHPKFQ\n")
    validate_output(str(f))  # must not raise
