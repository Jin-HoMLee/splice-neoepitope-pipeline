"""Unit tests for the HISAT2 read-input arg helper used by hisat2_align.

Guards the single-end (`-U`) vs paired-end (`-1/-2`) command selection. That
branch was previously inline shell in `alignment.smk` (unreachable by a test);
it is exercised in production but had no reproducible automated coverage after
the GCP patient_002 runs were retired - see Issue #962.
"""
import pytest

from hisat2_command import build_read_args


class TestBuildReadArgs:
    def test_single_end_uses_U(self):
        assert build_read_args("r1.fq.gz", "") == "-U r1.fq.gz"

    def test_single_end_none_fastq2(self):
        assert build_read_args("r1.fq.gz", None) == "-U r1.fq.gz"

    def test_whitespace_fastq2_is_single_end(self):
        # Matches strandness.get_strandness_from_row: whitespace-only fastq2 is SE.
        assert build_read_args("r1.fq.gz", "   ") == "-U r1.fq.gz"

    def test_paired_end_uses_1_2(self):
        assert build_read_args("r1.fq.gz", "r2.fq.gz") == "-1 r1.fq.gz -2 r2.fq.gz"

    def test_paired_end_with_directory_paths(self):
        # Realistic staged paths (mirrors the patient_002 BG003082_T0 PE case).
        r1 = "data/patient_002/BG003082_T0_R1.fastq.gz"
        r2 = "data/patient_002/BG003082_T0_R2.fastq.gz"
        assert build_read_args(r1, r2) == f"-1 {r1} -2 {r2}"

    def test_strips_surrounding_whitespace(self):
        assert build_read_args("  r1.fq.gz  ", "  r2.fq.gz  ") == "-1 r1.fq.gz -2 r2.fq.gz"

    def test_empty_fastq1_raises(self):
        with pytest.raises(ValueError):
            build_read_args("", "r2.fq.gz")

    def test_none_fastq1_raises(self):
        with pytest.raises(ValueError):
            build_read_args(None, "")
