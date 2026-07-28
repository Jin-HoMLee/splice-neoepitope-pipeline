"""Tests for the #737 sparsity tool (Issue #1230).

The point of the refactor is that no headline number is hand-maintained. So the two
tests that matter are the two that can actually catch drift:

  - `test_compute_stats_reproduces_committed_json` - the regenerator must reproduce
    the committed artifact exactly. Committing a distiller that no longer produces
    the committed bytes is how a stats file quietly becomes fiction.
  - `TestProseCanary` - the writeup's prose numbers must agree with stats.json, and
    the canary is itself proved able to fail by feeding it a deliberately stale text.

Issue #1069 is the case these exist for: the same headline numbers were hand-updated
across four files and a bot review caught a stale one the manual sweep missed.
"""

import json
import math
from pathlib import Path

import pytest

import sparsity

EXPERIMENT_DIR = Path(__file__).resolve().parent.parent
REGISTRY = EXPERIMENT_DIR.parent / "registry.tsv"
COMMITTED_STATS = EXPERIMENT_DIR / "outputs" / "sparsity_stats.json"
WRITEUP = EXPERIMENT_DIR / "sparsity_writeup.md"


class TestConcentration:
    """Herfindahl-Hirschman index and the inverse-Simpson effective count."""

    def test_even_four_way_split_has_effective_n_four(self):
        result = sparsity.concentration({"a": 5, "b": 5, "c": 5, "d": 5})
        assert result["hhi"] == pytest.approx(0.25)
        assert result["effective_n"] == pytest.approx(4.0)

    def test_single_category_is_maximally_concentrated(self):
        result = sparsity.concentration({"only": 7})
        assert result["hhi"] == pytest.approx(1.0)
        assert result["effective_n"] == pytest.approx(1.0)
        assert result["top_share"] == pytest.approx(1.0)
        # With one category there is no second, so top2 cannot exceed top.
        assert result["top2_share"] == pytest.approx(1.0)

    def test_skewed_split_shares_and_index(self):
        result = sparsity.concentration({"big": 3, "small": 1})
        assert result["n"] == 4
        assert result["distinct"] == 2
        assert result["top"] == "big"
        assert result["top_share"] == pytest.approx(0.75)
        assert result["top2_share"] == pytest.approx(1.0)
        assert result["hhi"] == pytest.approx(0.625)
        assert result["effective_n"] == pytest.approx(1.6)

    def test_effective_n_is_the_reciprocal_of_hhi(self):
        # Reported values are rounded per the stats.json contract (hhi 3dp,
        # effective_n 2dp), and effective_n is derived from the UNROUNDED hhi. So the
        # identity is asserted at the contract's precision, against the exact
        # unrounded reciprocal rather than against the rounded hhi.
        counts = {"a": 9, "b": 4, "c": 2, "d": 1}
        total = sum(counts.values())
        exact_hhi = sum((c / total) ** 2 for c in counts.values())
        result = sparsity.concentration(counts)
        assert result["hhi"] == round(exact_hhi, 3)
        assert result["effective_n"] == round(1.0 / exact_hhi, 2)

    def test_ties_pick_a_deterministic_top(self):
        # Set/dict iteration order must not decide the reported top contributor.
        first = sparsity.concentration({"b": 4, "a": 4})
        second = sparsity.concentration({"a": 4, "b": 4})
        assert first["top"] == second["top"]

    def test_empty_counts_raises_rather_than_returning_nan(self):
        with pytest.raises(ValueError):
            sparsity.concentration({})


class TestComputeStats:
    def test_compute_stats_reproduces_committed_json(self):
        """The regenerator must reproduce the committed artifact exactly.

        This is the known-answer control for the whole tool. If the registry grows and
        this fails, the fix is to re-run the tool and commit the new stats, never to
        relax the assertion.
        """
        committed = json.loads(COMMITTED_STATS.read_text())
        computed = sparsity.compute_stats(REGISTRY)
        assert computed == committed

    def test_scorable_positives_are_the_stats_denominator(self):
        computed = sparsity.compute_stats(REGISTRY)
        n = computed["scorable_positives"]
        for facet in ("by_study", "by_allele", "by_mechanism"):
            assert computed[facet]["n"] == n
            assert sum(computed[facet]["counts"].values()) == n

    def test_dual_gate_rows_do_not_exceed_total_rows(self):
        computed = sparsity.compute_stats(REGISTRY)
        assert computed["scorable_positives"] <= computed["dual_gate_rows"]
        assert computed["dual_gate_rows"] <= computed["registry_rows_total"]


class TestPowerCurve:
    def test_power_requirement_decreases_as_effect_size_grows(self):
        curve = sparsity.power_n_per_arm([0.70, 0.75, 0.80, 0.85])
        values = [curve[k] for k in sorted(curve, key=float)]
        assert values == sorted(values, reverse=True), (
            "a larger true AUC must need fewer samples per arm"
        )

    def test_matches_the_published_writeup_figures(self):
        curve = sparsity.power_n_per_arm([0.70, 0.75])
        assert curve["0.7"] == 31
        assert curve["0.75"] == 19


class TestProseCanary:
    """AC3: a stale number in a prose reader must be caught."""

    def test_the_committed_writeup_agrees_with_the_committed_stats(self):
        stats = json.loads(COMMITTED_STATS.read_text())
        mismatches = sparsity.check_prose_consistency(WRITEUP.read_text(), stats)
        assert mismatches == []

    def test_canary_catches_a_stale_scorable_count(self):
        # The falsifier: if this passes, the canary cannot detect drift and is theatre.
        stats = json.loads(COMMITTED_STATS.read_text())
        stale = WRITEUP.read_text().replace("82 scorable", "80 scorable")
        mismatches = sparsity.check_prose_consistency(stale, stats)
        assert mismatches, "canary failed to notice a changed scorable-positive count"

    def test_canary_catches_a_stale_effective_allele_count(self):
        stats = json.loads(COMMITTED_STATS.read_text())
        stale = WRITEUP.read_text().replace("1.26", "1.99")
        mismatches = sparsity.check_prose_consistency(stale, stats)
        assert mismatches, "canary failed to notice a changed effective allele count"

    def test_canary_reports_which_claim_drifted(self):
        stats = json.loads(COMMITTED_STATS.read_text())
        stale = WRITEUP.read_text().replace("1.26", "1.99")
        mismatches = sparsity.check_prose_consistency(stale, stats)
        assert any("allele" in m.lower() for m in mismatches), mismatches


class TestAllReaders:
    """Every reader that quotes a headline number is held to stats.json."""

    def test_all_three_committed_readers_are_consistent(self):
        stats = json.loads(COMMITTED_STATS.read_text())
        results = sparsity.check_all_readers(stats)
        assert set(results) == {"sparsity_writeup.md", "slides.qmd", "registry README.md"}
        assert all(v == [] for v in results.values()), results

    def test_every_occurrence_is_checked_not_just_the_first(self):
        """A number stated twice, stale in one place only, must still be caught.

        This is the Issue #1069 shape: the same headline repeated across readers and
        updated in some of them. A checker that stops at the first match reports clean.
        """
        stats = json.loads(COMMITTED_STATS.read_text())
        text = (
            "Of 97 registry rows, 95 pass both gates; **82 are scorable**. "
            "Later on, we restate it wrongly as 80 scorable positives."
        )
        mismatches = sparsity.check_prose_consistency(
            text, stats, ["registry rows total", "dual-gate rows", "scorable positives"]
        )
        assert any("scorable" in m for m in mismatches), mismatches

    def test_a_missing_claim_is_reported_not_skipped(self):
        # Deleting the sentence must not be a way to make the check pass.
        stats = json.loads(COMMITTED_STATS.read_text())
        mismatches = sparsity.check_prose_consistency(
            "prose with no numbers at all", stats, ["scorable positives"]
        )
        assert any("not found" in m for m in mismatches), mismatches

    def test_slides_stale_number_is_caught(self):
        stats = json.loads(COMMITTED_STATS.read_text())
        stale = (EXPERIMENT_DIR / "slides.qmd").read_text().replace("1.26", "1.99")
        mismatches = sparsity.check_prose_consistency(
            stale, stats, sparsity.READER_CLAIMS["slides.qmd"]
        )
        assert any("allele" in m.lower() for m in mismatches), mismatches


class TestNumberWords:
    """The slides spell small counts out; the canary must still read them."""

    def test_number_word_matches_its_digit(self):
        stats = json.loads(COMMITTED_STATS.read_text())
        text = "**One** hard true-negative field-wide."
        assert sparsity.check_prose_consistency(
            text, stats, ["hard true-negative count"]
        ) == []

    def test_wrong_number_word_is_still_caught(self):
        stats = json.loads(COMMITTED_STATS.read_text())
        text = "**Three** hard true-negative rows field-wide."
        assert sparsity.check_prose_consistency(
            text, stats, ["hard true-negative count"]
        ) != []


def test_hhi_of_a_uniform_split_scales_as_one_over_k():
    # Property check across sizes; guards the formula rather than one hand-picked case.
    for k in (2, 3, 5, 10):
        result = sparsity.concentration({str(i): 1 for i in range(k)})
        # Compared at the contract's rounding precision; 1/3 reports as 0.333.
        assert result["hhi"] == round(1.0 / k, 3)
        assert result["effective_n"] == round(float(k), 2)
        assert not math.isnan(result["hhi"])
