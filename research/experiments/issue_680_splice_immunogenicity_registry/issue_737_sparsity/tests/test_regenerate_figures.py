"""Tests for the #737 deck-figure regenerator (Issue #1230, option C).

The figures render numbers visually, so they are a reader the prose canary structurally
cannot check - it reads text, not pixels. What substitutes for that check is this: the
regenerator derives every figure from `sparsity_stats.json`, the same single source of
truth the prose is held to, and never recomputes from `registry.tsv`. These tests pin
that property, plus the pure Pareto arithmetic the first figure needs.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "figures"))

import _regenerate_figures as figs  # noqa: E402

EXPERIMENT_DIR = Path(__file__).resolve().parent.parent
STATS = EXPERIMENT_DIR / "outputs" / "sparsity_stats.json"


class TestParetoSeries:
    def test_cumulative_share_ends_at_one(self):
        labels, values, cumulative = figs.pareto_series({"a": 3, "b": 1})
        assert labels == ["a", "b"]
        assert values == [3, 1]
        assert cumulative[-1] == pytest.approx(1.0)

    def test_cumulative_is_monotonic_non_decreasing(self):
        _, _, cumulative = figs.pareto_series({"a": 5, "b": 3, "c": 2, "d": 1})
        assert cumulative == sorted(cumulative)

    def test_second_point_matches_the_top_two_share(self):
        # The figure annotates the top-two share at index 1, so this must agree with
        # what `concentration` reports, or the annotation contradicts the caption.
        counts = {"a": 35, "b": 18, "c": 10, "d": 5}
        _, _, cumulative = figs.pareto_series(counts)
        expected = (35 + 18) / sum(counts.values())
        assert cumulative[1] == pytest.approx(expected)

    def test_order_is_descending_by_count(self):
        labels, values, _ = figs.pareto_series({"small": 1, "big": 9, "mid": 4})
        assert values == [9, 4, 1]
        assert labels == ["big", "mid", "small"]

    def test_empty_counts_raises(self):
        with pytest.raises(ValueError):
            figs.pareto_series({})


class TestPowerKeyLookup:
    """The notebook used float keys; JSON gives strings. Regression guard."""

    def test_power_lookup_accepts_the_json_string_keys(self):
        stats = json.loads(STATS.read_text())
        assert figs.power_at(stats, 0.75) == stats["power_n_per_arm"]["0.75"]
        assert figs.power_at(stats, 0.70) == stats["power_n_per_arm"]["0.7"]

    def test_power_lookup_raises_on_an_absent_auc(self):
        stats = json.loads(STATS.read_text())
        with pytest.raises(KeyError):
            figs.power_at(stats, 0.99)


class TestRegenerate:
    def test_writes_all_three_deck_figures(self, tmp_path):
        stats = json.loads(STATS.read_text())
        written = figs.regenerate(stats, tmp_path)
        assert {p.name for p in written} == {
            "fig_study_pareto.png",
            "fig_allele_mechanism.png",
            "fig_negative_scarcity.png",
        }
        for path in written:
            assert path.exists()
            assert path.stat().st_size > 5_000, f"{path.name} looks empty"

    def test_regenerate_reads_only_the_stats_it_is_given(self, tmp_path):
        """Figures must derive from stats.json, not recompute from the registry.

        Handing it altered stats must change the output. If a figure silently reread
        registry.tsv, this would produce identical bytes and pass vacuously.
        """
        stats = json.loads(STATS.read_text())
        baseline = figs.regenerate(stats, tmp_path / "base")
        baseline_bytes = {p.name: p.read_bytes() for p in baseline}

        altered = json.loads(STATS.read_text())
        altered["negatives"]["hard_true_negative"] = 42
        altered["negatives"]["tested_total"] = 50
        changed = figs.regenerate(altered, tmp_path / "altered")
        changed_bytes = {p.name: p.read_bytes() for p in changed}

        assert (
            changed_bytes["fig_negative_scarcity.png"]
            != baseline_bytes["fig_negative_scarcity.png"]
        ), "negative-scarcity figure ignored the stats it was handed"
