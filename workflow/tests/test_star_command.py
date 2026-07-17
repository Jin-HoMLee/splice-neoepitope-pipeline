"""Unit tests for the STAR --readFilesIn arg helper used by star_align.

Guards the single-end vs paired-end file-list selection. That branch was
previously inline shell in `alignment.smk` (unreachable by a test); this is the
STAR sibling of `test_hisat2_command.py`, extracted as the Issue #1081 symmetry
follow-up to #962. STAR takes a bare space-separated file list, not `-U`/`-1/-2`.
"""
import pytest

from star_command import build_read_files_in


class TestBuildReadFilesIn:
    def test_single_end_is_r1_only(self):
        assert build_read_files_in("r1.fq.gz", "") == "r1.fq.gz"

    def test_single_end_none_fastq2(self):
        assert build_read_files_in("r1.fq.gz", None) == "r1.fq.gz"

    def test_whitespace_fastq2_is_single_end(self):
        # Mirrors build_read_args: whitespace-only fastq2 is single-end.
        assert build_read_files_in("r1.fq.gz", "   ") == "r1.fq.gz"

    def test_paired_end_is_space_separated(self):
        assert build_read_files_in("r1.fq.gz", "r2.fq.gz") == "r1.fq.gz r2.fq.gz"

    def test_paired_end_with_directory_paths(self):
        r1 = "data/patient_002/BG003082_T0_R1.fastq.gz"
        r2 = "data/patient_002/BG003082_T0_R2.fastq.gz"
        assert build_read_files_in(r1, r2) == f"{r1} {r2}"

    def test_strips_surrounding_whitespace(self):
        assert build_read_files_in("  r1.fq.gz  ", "  r2.fq.gz  ") == "r1.fq.gz r2.fq.gz"

    def test_empty_fastq1_raises(self):
        with pytest.raises(ValueError):
            build_read_files_in("", "r2.fq.gz")

    def test_none_fastq1_raises(self):
        with pytest.raises(ValueError):
            build_read_files_in(None, "")
