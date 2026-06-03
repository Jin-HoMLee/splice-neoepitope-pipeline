"""Offline unit tests for the morning-news poller (Issue #639).

All network / gh calls are dependency-injected; nothing here touches the wire.
The headline test is ``test_repetition_dies_on_second_run`` — the property the
whole feature exists to guarantee.
"""
import pytest

import poll_releases as pr


# --------------------------------------------------------------------------- #
# version parsing / comparison
# --------------------------------------------------------------------------- #
def test_version_tuple_strips_local_and_prefix():
    assert pr.version_tuple("2.12.0+cu126") == (2, 12, 0)
    assert pr.version_tuple("v9.22.0") == (9, 22, 0)
    assert pr.version_tuple("1.3.5") == (1, 3, 5)


def test_version_gt_orders_numerically():
    assert pr.version_gt("0.7.0", "0.6.1") is True
    assert pr.version_gt("0.6.1", "0.6.1") is False
    assert pr.version_gt("9.22.0", "9.9.0") is True  # not lexicographic


# --------------------------------------------------------------------------- #
# guards
# --------------------------------------------------------------------------- #
def test_guard_frozen_suppresses():
    tool = {"tool": "jax", "frozen": True}
    assert pr.apply_guards(tool, "0.4.30") == (None, "frozen")


def test_guard_watch_only_is_a_watch_not_a_bump():
    tool = {"tool": "gcp_dlvm_image", "watch_only": True}
    assert pr.apply_guards(tool, "v20260601") == ("v20260601", "watch")


def test_guard_max_version_caps_imgtgenedl():
    tool = {"tool": "IMGTgeneDL", "max_version": "0.6.1"}
    # PyPI publishing 0.7.0 must NOT flag (conda solve breaks past the cap).
    assert pr.apply_guards(tool, "0.7.0") == (None, "above-cap")
    # At the cap is fine.
    assert pr.apply_guards(tool, "0.6.1") == ("0.6.1", "bump")


def test_guard_fetch_failure():
    assert pr.apply_guards({"tool": "snakemake"}, None) == (None, "fetch-failed")


def test_guard_normal_bump():
    assert pr.apply_guards({"tool": "snakemake"}, "9.22.0") == ("9.22.0", "bump")


# --------------------------------------------------------------------------- #
# compute_delta
# --------------------------------------------------------------------------- #
def test_first_sighting_seeds_baseline_silently():
    rec = pr.compute_delta({"tool": "snakemake"}, "9.22.0", None)
    assert rec["kind"] == "baseline"
    assert rec["surface"] is False
    assert rec["watermark"] == "9.22.0"


def test_genuine_delta_surfaces_and_advances_watermark():
    rec = pr.compute_delta({"tool": "snakemake"}, "9.22.0", "9.21.0")
    assert rec["surface"] is True
    assert rec["kind"] == "bump"
    assert (rec["from"], rec["to"]) == ("9.21.0", "9.22.0")
    assert rec["watermark"] == "9.22.0"


def test_unchanged_is_muted():
    rec = pr.compute_delta({"tool": "snakemake"}, "9.22.0", "9.22.0")
    assert rec["surface"] is False
    assert rec["watermark"] == "9.22.0"


def test_frozen_delta_never_surfaces_and_leaves_watermark_untouched():
    rec = pr.compute_delta({"tool": "jax", "frozen": True}, "0.4.30", "0.3.25")
    assert rec["surface"] is False
    assert rec["suppressed"] is True
    assert rec["watermark"] is None  # do not advance a frozen pin


def test_capped_delta_suppressed():
    rec = pr.compute_delta({"tool": "IMGTgeneDL", "max_version": "0.6.1"}, "0.7.0", "0.6.1")
    assert rec["surface"] is False
    assert rec["suppressed"] is True


# --------------------------------------------------------------------------- #
# torch cu126 channel parsing
# --------------------------------------------------------------------------- #
def test_parse_cu126_index_keeps_tag_and_picks_highest():
    html = (
        '<a href="torch-2.11.0+cu126-cp311.whl">torch-2.11.0+cu126</a>'
        '<a href="torch-2.12.0+cu126-cp311.whl">torch-2.12.0+cu126</a>'
    )
    assert pr.parse_cu126_index(html) == "2.12.0+cu126"


def test_parse_cu126_index_empty():
    assert pr.parse_cu126_index("<html>nothing here</html>") is None


# --------------------------------------------------------------------------- #
# run() — end-to-end with injected fetcher + tracker
# --------------------------------------------------------------------------- #
def _fetcher(versions):
    def _f(tool):
        return versions.get(tool["tool"])
    return _f


def test_run_surfaces_only_genuine_deltas_and_updates_watermark():
    basket = {"software": [
        {"tool": "snakemake", "feed_type": "pypi", "feed": "snakemake"},
        {"tool": "bedtools", "feed_type": "github", "feed": "arq5x/bedtools2"},
        {"tool": "jax", "feed_type": "pypi", "feed": "jax", "frozen": True},
        {"tool": "IMGTgeneDL", "feed_type": "pypi", "feed": "IMGTgeneDL", "max_version": "0.6.1"},
    ]}
    versions = {
        "snakemake": "9.22.0",   # delta vs 9.21.0
        "bedtools": "2.31.1",    # unchanged
        "jax": "0.4.30",         # frozen -> suppressed
        "IMGTgeneDL": "0.7.0",   # above cap -> suppressed
    }
    watermark = {"tools": {"snakemake": "9.21.0", "bedtools": "2.31.1"}}

    surfaced, wm = pr.run(
        basket, watermark, check_tracked=False, fetcher=_fetcher(versions),
    )

    names = [r["tool"] for r in surfaced]
    assert names == ["snakemake"]
    assert wm["tools"]["snakemake"] == "9.22.0"
    assert wm["tools"]["bedtools"] == "2.31.1"
    # frozen / capped tools never enter the watermark
    assert "jax" not in wm["tools"]
    assert "IMGTgeneDL" not in wm["tools"]


def test_pure_delta_tracked_item_still_surfaces_with_annotation():
    basket = {"software": [{"tool": "snakemake", "feed_type": "pypi", "feed": "snakemake"}]}
    watermark = {"tools": {"snakemake": "9.21.0"}}

    surfaced, _ = pr.run(
        basket, watermark, check_tracked=True,
        fetcher=_fetcher({"snakemake": "9.22.0"}),
        tracker=lambda name: 200,
    )
    assert len(surfaced) == 1
    # Pure-delta: a genuine bump surfaces even though Issue #200 tracks it,
    # but it is annotated with the tracking issue.
    assert surfaced[0]["tracking_issue"] == 200


def test_watch_only_change_surfaces_as_watch_without_tracking_lookup():
    basket = {"software": [
        {"tool": "gcp_dlvm_image", "feed_type": "pypi", "feed": "x", "watch_only": True},
    ]}
    watermark = {"tools": {"gcp_dlvm_image": "img-A"}}
    called = []

    surfaced, _ = pr.run(
        basket, watermark, check_tracked=True,
        fetcher=_fetcher({"gcp_dlvm_image": "img-B"}),
        tracker=lambda name: called.append(name) or 999,
    )
    assert surfaced[0]["kind"] == "watch"
    assert "tracking_issue" not in surfaced[0]  # watch items skip the tracker
    assert called == []


def test_repetition_dies_on_second_run():
    """The headline property: a tool surfaces once, never again until it moves."""
    basket = {"software": [{"tool": "snakemake", "feed_type": "pypi", "feed": "snakemake"}]}
    fetch = _fetcher({"snakemake": "9.22.0"})
    watermark = {"tools": {"snakemake": "9.21.0"}}

    first, watermark = pr.run(basket, watermark, check_tracked=False, fetcher=fetch)
    assert [r["tool"] for r in first] == ["snakemake"]

    # Same upstream version, watermark carried forward -> silent.
    second, watermark = pr.run(basket, watermark, check_tracked=False, fetcher=fetch)
    assert second == []


def test_gap_yields_single_clean_delta_not_a_repeat():
    """Skipping runs while upstream moves 9.21 -> 9.25 gives ONE delta, not stale repeats."""
    basket = {"software": [{"tool": "snakemake", "feed_type": "pypi", "feed": "snakemake"}]}
    watermark = {"tools": {"snakemake": "9.21.0"}}

    surfaced, watermark = pr.run(
        basket, watermark, check_tracked=False, fetcher=_fetcher({"snakemake": "9.25.0"}),
    )
    assert len(surfaced) == 1
    assert (surfaced[0]["from"], surfaced[0]["to"]) == ("9.21.0", "9.25.0")


def test_reference_data_section_is_polled_too():
    basket = {
        "software": [],
        "reference_data": [{"tool": "vdjdb", "feed_type": "github", "feed": "antigenomics/vdjdb-db"}],
    }
    watermark = {"tools": {"vdjdb": "2026-05-16"}}
    surfaced, _ = pr.run(
        basket, watermark, check_tracked=False, fetcher=_fetcher({"vdjdb": "2026-06-01"}),
    )
    assert [r["tool"] for r in surfaced] == ["vdjdb"]


# --------------------------------------------------------------------------- #
# briefing formatting
# --------------------------------------------------------------------------- #
def test_format_briefing_empty_says_none():
    out = pr.format_briefing([], {"watch": []})
    assert "none" in out.lower()


def test_format_briefing_lists_delta_and_watchlist():
    surfaced = [{"tool": "snakemake", "kind": "bump", "from": "9.21.0", "to": "9.22.0", "tracking_issue": 200}]
    basket = {"watch": [{"tool": "uv"}, {"tool": "ruff"}]}
    out = pr.format_briefing(surfaced, basket)
    assert "snakemake" in out and "9.21.0" in out and "9.22.0" in out
    assert "Issue #200" in out
    assert "uv, ruff" in out


# --------------------------------------------------------------------------- #
# yaml IO round-trip (skipped if PyYAML is unavailable)
# --------------------------------------------------------------------------- #
def test_watermark_roundtrip(tmp_path):
    pytest.importorskip("yaml")
    path = tmp_path / "watermark_developer.yaml"
    data = {"tools": {"snakemake": "9.22.0"}, "last_briefing": "2026-06-03"}
    pr.save_watermark(path, data)
    loaded = pr.load_watermark(path)
    assert loaded["tools"]["snakemake"] == "9.22.0"
    assert loaded["last_briefing"] == "2026-06-03"


def test_load_watermark_missing_file_returns_empty(tmp_path):
    pytest.importorskip("yaml")
    loaded = pr.load_watermark(tmp_path / "does_not_exist.yaml")
    assert loaded == {"tools": {}, "last_briefing": None}


def test_basket_file_parses_and_has_known_guards():
    pytest.importorskip("yaml")
    basket = pr.load_basket()
    names = {t["tool"] for t in basket["software"]}
    assert {"snakemake", "torch", "IMGTgeneDL", "jax", "gcp_dlvm_image"} <= names
    by_name = {t["tool"]: t for t in basket["software"]}
    assert by_name["jax"].get("frozen") is True
    assert by_name["IMGTgeneDL"].get("max_version") == "0.6.1"
    # gcp_dlvm_image is suppressed by feed_type: manual (no machine feed), NOT by
    # a watch_only guard — manual feeds return None before apply_guards sees them.
    assert by_name["gcp_dlvm_image"]["feed_type"] == "manual"
    assert "watch_only" not in by_name["gcp_dlvm_image"]
    assert by_name["torch"]["feed_type"] == "pytorch_cu126"
