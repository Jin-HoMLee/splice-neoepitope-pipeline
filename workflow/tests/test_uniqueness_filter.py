"""Unit tests for the pure NH-uniqueness prefilter helpers (Issue #919).

Guards the helper logic. The sibling test_alignment_uniqueness_filter.py guards
that the `hisat2_align` rule actually calls these helpers, so an edit dropping
the block from the rule fails there rather than passing silently here.
"""
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from uniqueness_filter import (  # noqa: E402
    MIN_SAMTOOLS_VERSION,
    NH_UNIQUE_EXPR,
    build_prefilter_block,
    build_prefilter_gate,
    build_preflight_block,
    build_preflight_gate,
    filtered_bam_outputs,
    filtered_bam_path,
    is_enabled,
    parse_samtools_version,
    regtools_input_bam,
    supports_filter_expressions,
)

_BAM = "results/{patient_id}/alignment/{sample}/{sample}.bam"


# ── is_enabled: default off ──────────────────────────────────────────────────

@pytest.mark.parametrize("config", [
    {},
    {"alignment": {}},
    {"alignment": None},
    {"alignment": {"uniqueness_filter": {}}},
    {"alignment": {"uniqueness_filter": None}},
    {"alignment": {"uniqueness_filter": {"enabled": False}}},
])
def test_is_enabled_defaults_off(config):
    """An untouched (or partially-specified) config must not enable the filter."""
    assert is_enabled(config) is False


def test_is_enabled_true_when_set():
    assert is_enabled({"alignment": {"uniqueness_filter": {"enabled": True}}}) is True


# ── filtered_bam_path ────────────────────────────────────────────────────────

def test_filtered_bam_path_replaces_bam_suffix():
    assert filtered_bam_path("a/b/s.bam") == "a/b/s.nh_unique.bam"


def test_filtered_bam_path_preserves_wildcards():
    """Snakemake wildcards must survive untouched - they expand at job time."""
    assert filtered_bam_path(_BAM) == (
        "results/{patient_id}/alignment/{sample}/{sample}.nh_unique.bam"
    )


def test_filtered_bam_path_appends_when_no_bam_suffix():
    assert filtered_bam_path("a/b/s") == "a/b/s.nh_unique.bam"


# ── samtools version gate ────────────────────────────────────────────────────

@pytest.mark.parametrize("text,expected", [
    ("samtools 1.13-4\nUsing htslib 1.13", (1, 13)),   # Ubuntu 22.04 apt build
    ("samtools 1.23.1\nUsing htslib 1.23.1", (1, 23)),  # local homebrew build
    ("samtools 1.12", (1, 12)),
    ("samtools 1.9", (1, 9)),
    ("samtools 2", (2, 0)),
    ("", None),
    ("bcftools 1.20", None),
    ("samtools", None),
    ("samtools notaversion", None),
])
def test_parse_samtools_version(text, expected):
    assert parse_samtools_version(text) == expected


@pytest.mark.parametrize("version,expected", [
    ((1, 11), False),   # one release before -e landed
    ((1, 12), True),    # -e introduced here
    ((1, 13), True),    # apt / production
    ((1, 23), True),    # local
    ((2, 0), True),
    (None, False),      # unparseable version must not be assumed capable
])
def test_supports_filter_expressions(version, expected):
    assert supports_filter_expressions(version) is expected


def test_min_version_matches_upstream_introduction():
    """`-e` landed in samtools 1.12; the constant must not drift below that."""
    assert MIN_SAMTOOLS_VERSION == (1, 12)


# ── block composition: off ───────────────────────────────────────────────────

def test_both_gates_empty_when_disabled():
    """Off must contribute no shell at all, so the command is a provable no-op."""
    assert build_preflight_gate(False) == ""
    assert build_prefilter_gate(False) == ""


def test_regtools_reads_unfiltered_bam_when_disabled():
    assert regtools_input_bam(False) == "{output.bam}"


def test_no_filtered_outputs_when_disabled():
    assert filtered_bam_outputs(False, _BAM) == {}


# ── block composition: on ────────────────────────────────────────────────────

def test_regtools_reads_filtered_bam_when_enabled():
    assert regtools_input_bam(True) == "{output.filtered_bam}"


def test_filtered_outputs_when_enabled():
    out = filtered_bam_outputs(True, _BAM)
    assert set(out) == {"filtered_bam", "filtered_bai"}
    assert out["filtered_bam"].endswith("{sample}.nh_unique.bam")
    assert out["filtered_bai"] == out["filtered_bam"] + ".bai"


def test_enabled_gates_carry_preflight_and_prefilter_separately():
    """The two gates are emitted at different points in the shell: the preflight
    before the aligner (so a bad samtools fails fast) and the prefilter after it."""
    assert build_preflight_gate(True) == build_preflight_block()
    assert build_prefilter_gate(True) == build_prefilter_block()


def test_preflight_gate_carries_no_prefilter():
    """Regression guard on the hoist: if the preflight gate also emitted the
    prefilter, the filter would run before the BAM it reads exists."""
    gate = build_preflight_gate(True)
    assert "--expr" in gate
    assert "{output.filtered_bam}" not in gate


def test_enabled_gate_uses_nh_expression_not_mapq_floor():
    """NH is the primary lever; a `-q` MAPQ floor would drop unique low-MAPQ reads."""
    gate = build_prefilter_gate(True)
    assert "-e '[NH]==1'" in gate
    assert "-q 2" not in gate


def test_preflight_probes_capability_and_exits_nonzero():
    """The preflight must abort, not warn - a silent fallback emits unfiltered BED12."""
    block = build_preflight_block()
    assert "--expr" in block
    assert "exit 1" in block


def test_prefilter_writes_and_indexes_the_filtered_bam():
    block = build_prefilter_block()
    assert "-o {output.filtered_bam}" in block
    assert "samtools index {output.filtered_bam}" in block


# ── brace safety: these strings are Snakemake templates ──────────────────────

def test_filter_expression_is_brace_free():
    """A `{` in the expression would be read as a Snakemake placeholder."""
    assert "{" not in NH_UNIQUE_EXPR and "}" not in NH_UNIQUE_EXPR


def test_every_placeholder_in_enabled_block_is_resolvable():
    """Guard against a typo'd placeholder that Snakemake would fail to expand.

    Snakemake expands `{...}` against the job namespace, so any brace group in
    the block must name something the `hisat2_align` rule actually provides.
    """
    allowed = {"log", "threads", "output.bam", "output.filtered_bam"}
    shell = build_preflight_gate(True) + build_prefilter_gate(True)
    found = set(re.findall(r"\{([^{}]*)\}", shell))
    assert found <= allowed, f"unresolvable placeholders: {sorted(found - allowed)}"
