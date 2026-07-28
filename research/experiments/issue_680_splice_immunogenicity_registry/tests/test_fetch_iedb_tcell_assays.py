"""Integrity tests for the Issue #1290 IEDB reverse-discovery puller.

The deliverable of Issue #1290 is a *count*: how many assayed epitopes reverse-map to
candidate neojunctions, split positive / negative. Every failure mode below corrupts
that count silently rather than loudly, which is why each one is asserted rather than
assumed:

  - a truncated page walk under-reports yield and looks like a clean run,
  - an overlapping page walk double-counts,
  - an unrecognised IEDB outcome value bucketed by a permissive default would move
    rows between the positive and negative classes, and the negative class is the
    whole reason this Issue exists.

Every check here is paired: one case that must pass and one that must fail, so a
check that could only ever confirm cannot survive.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "issue_1290_reverse_discovery"))

import fetch_iedb_tcell_assays as fetch  # noqa: E402


class TestNormaliseOutcome:
    """IEDB's `qualitative_measure` vocabulary -> the registry's binary label."""

    @pytest.mark.parametrize(
        "measure",
        ["Positive", "Positive-High", "Positive-Intermediate", "Positive-Low"],
    )
    def test_every_positive_variant_maps_to_positive(self, measure):
        assert fetch.normalise_outcome(measure) == "positive"

    def test_negative_maps_to_negative(self):
        assert fetch.normalise_outcome("Negative") == "negative"

    def test_unknown_measure_raises_rather_than_defaulting(self):
        # The falsifier for the whole vocabulary: if IEDB adds a category, this must
        # break the run instead of quietly folding it into one of the two classes.
        with pytest.raises(fetch.IntegrityError, match="unrecognised"):
            fetch.normalise_outcome("Equivocal")

    def test_empty_measure_raises(self):
        with pytest.raises(fetch.IntegrityError):
            fetch.normalise_outcome("")


class TestPlanPages:
    def test_covers_a_partial_final_page(self):
        assert fetch.plan_pages(total=25, page_size=10) == [0, 10, 20]

    def test_exact_multiple_has_no_empty_trailing_page(self):
        assert fetch.plan_pages(total=20, page_size=10) == [0, 10]

    def test_total_smaller_than_page_size_is_one_page(self):
        assert fetch.plan_pages(total=3, page_size=10) == [0]

    def test_zero_total_plans_no_pages(self):
        assert fetch.plan_pages(total=0, page_size=10) == []


class TestAssertNoTruncation:
    def test_passes_when_counts_agree(self):
        fetch.assert_no_truncation(fetched=66358, expected_total=66358)

    def test_raises_when_short(self):
        with pytest.raises(fetch.IntegrityError, match="66357"):
            fetch.assert_no_truncation(fetched=66357, expected_total=66358)

    def test_raises_when_over(self):
        # Over-fetching means the page walk overlapped; just as wrong as truncation.
        with pytest.raises(fetch.IntegrityError):
            fetch.assert_no_truncation(fetched=66359, expected_total=66358)


class TestAssertUnique:
    def test_passes_on_distinct_ids(self):
        rows = [{"tcell_id": 1}, {"tcell_id": 2}, {"tcell_id": 3}]
        fetch.assert_unique(rows, key="tcell_id")

    def test_raises_and_names_the_duplicate(self):
        rows = [{"tcell_id": 1}, {"tcell_id": 2}, {"tcell_id": 1}]
        with pytest.raises(fetch.IntegrityError, match="1"):
            fetch.assert_unique(rows, key="tcell_id")
