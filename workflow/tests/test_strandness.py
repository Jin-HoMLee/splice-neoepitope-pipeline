"""Unit tests for the strandness helper used by hisat2_align."""
import pytest

from strandness import get_strandness_flag, get_strandness_from_row


class TestGetStrandnessFlag:
    @pytest.mark.parametrize(
        "strandness,is_paired_end,expected",
        [
            ("forward", False, "F"),
            ("forward", True, "FR"),
            ("reverse", False, "R"),
            ("reverse", True, "RF"),
            ("unstranded", False, ""),
            ("unstranded", True, ""),
        ],
    )
    def test_recognized_values(self, strandness, is_paired_end, expected):
        assert get_strandness_flag(strandness, is_paired_end) == expected

    @pytest.mark.parametrize("missing_or_empty", [None, "", "  "])
    def test_missing_or_empty_means_unstranded(self, missing_or_empty):
        # Backward compat: samples.tsv without the column or with blank cells must
        # default to unstranded (no flag passed) — preserves current pipeline behavior.
        assert get_strandness_flag(missing_or_empty, False) == ""
        assert get_strandness_flag(missing_or_empty, True) == ""

    def test_case_insensitive(self):
        assert get_strandness_flag("Forward", False) == "F"
        assert get_strandness_flag("REVERSE", True) == "RF"

    @pytest.mark.parametrize("bad", ["garbage", "forwrd", "F", "RF", "fwd"])
    def test_unrecognized_value_raises(self, bad):
        # Typos like `forwrd` would silently produce unstranded alignment, masking
        # a misconfigured sample. HISAT2-flag literals (`F`, `RF`) are also rejected —
        # the column takes biological direction, not the HISAT2 flag value.
        with pytest.raises(ValueError, match="Unrecognized strandness value"):
            get_strandness_flag(bad, False)


class TestGetStrandnessFromRow:
    def test_forward_se_row(self):
        # PBMC scRNA R2-only row from patient_002.tsv — SE forward-stranded.
        row = {"sample_id": "PBMC_scRNA_Pool1_L002", "fastq1": "x.fastq.gz", "fastq2": "", "strandness": "forward"}
        assert get_strandness_from_row(row) == "F"

    def test_unstranded_pe_row(self):
        # BG003082 tumor row from patient_002.tsv — PE, no strandness column set.
        row = {"sample_id": "BG003082_T0", "fastq1": "r1.fq.gz", "fastq2": "r2.fq.gz", "strandness": ""}
        assert get_strandness_from_row(row) == ""

    def test_missing_strandness_key(self):
        # samples.tsv predating the column — `strandness` key absent entirely.
        row = {"sample_id": "SRR9143066", "fastq1": "x.fq.gz", "fastq2": ""}
        assert get_strandness_from_row(row) == ""

    def test_whitespace_fastq2_is_se(self):
        # Whitespace-only fastq2 must be treated as SE, not PE.
        row = {"sample_id": "s", "fastq1": "x.fq.gz", "fastq2": "   ", "strandness": "reverse"}
        assert get_strandness_from_row(row) == "R"

    def test_reverse_pe_row(self):
        # Hypothetical PE reverse-stranded sample (TruSeq dUTP-style library).
        row = {"sample_id": "s", "fastq1": "r1.fq.gz", "fastq2": "r2.fq.gz", "strandness": "reverse"}
        assert get_strandness_from_row(row) == "RF"
