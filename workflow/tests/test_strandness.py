"""Unit tests for the strandness helper used by hisat2_align."""
import pytest

from strandness import get_strandness_flag


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

    def test_unrecognized_value_means_unstranded(self):
        # Defensive: unknown strings fall through to "" rather than raising. Failure
        # mode is graceful unstranded fallback, not a crashed alignment rule.
        assert get_strandness_flag("garbage", False) == ""
        # HISAT2-flag literal (not a biological direction) is also rejected.
        assert get_strandness_flag("F", False) == ""
