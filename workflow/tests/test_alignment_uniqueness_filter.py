"""Snapshot test of the rendered hisat2_align uniqueness-filter block (Issue #919).

Complements the pure-helper test (test_uniqueness_filter.py) at the rule layer:
it asserts the `hisat2_align` rule actually wires the prefilter, so a future
edit that drops the block from the shell (or stops feeding regtools the
filtered BAM) fails here even though the pure tests still pass. Sibling of
test_alignment_hisat2_command.py.

The load-bearing assertion is the *off* case: the knob defaults off, and an
opt-in filter that silently perturbs the default DAG would be a correctness
bug, not a feature. So the off case asserts no filtered BAM is declared, no
`samtools view` prefilter runs, and regtools still reads the unfiltered BAM.

Runs `snakemake --dry-run --printshellcmds` programmatically with an
aligner=hisat2 override; skips where snakemake is unavailable.
"""
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

_TARGET = "results/tp/alignment/se_sample/raw_junctions.tsv"
_BAM = "results/tp/alignment/se_sample/se_sample.bam"
_FILTERED_BAM = "results/tp/alignment/se_sample/se_sample.nh_unique.bam"


def _render(tmp_path, enabled):
    """Dry-run the hisat2_align rule with the knob on or off; return stdout."""
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not on PATH - activate the snakemake conda env")

    r1 = tmp_path / "se_R1.fq.gz"
    r1.touch()
    index_dir = tmp_path / "hisat2_index"
    index_dir.mkdir()
    (index_dir / "index.done").touch()
    (index_dir / "genome.1.ht2").touch()

    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "patient_id\tsample_id\tsample_type\tfastq1\tfastq2\tstrandness\n"
        f"tp\tse_sample\tPrimary Tumor\t{r1}\t\t\n"
    )

    override = tmp_path / "override.yaml"
    override.write_text(textwrap.dedent(f"""\
        samples_tsv: "{samples}"
        alignment:
          aligner: "hisat2"
          hisat2_index_dir: "{index_dir}"
          hisat2_prebuilt_url: ""
          uniqueness_filter:
            enabled: {str(bool(enabled)).lower()}
    """))

    result = subprocess.run(
        ["snakemake", "-n", "--printshellcmds",
         "--configfile", "config/config.yaml", str(override), "--", _TARGET],
        cwd=REPO_ROOT, capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        pytest.fail(
            f"snakemake dry-run failed (enabled={enabled}, exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )
    return result.stdout


def regtools_input_bam(stdout):
    """Return the positional BAM `regtools junctions extract` is invoked on.

    Asserting on a bare `"<path> \\"` substring is not enough: the prefilter's
    own `-o <filtered_bam> \\` line matches it too, so such a test passes even
    when regtools is wired to the wrong BAM (caught by mutating the helper).
    This isolates the regtools invocation and returns its positional argument.
    """
    lines = stdout.splitlines()
    for i, line in enumerate(lines):
        if "regtools junctions extract" not in line:
            continue
        for cont in lines[i + 1:]:
            token = cont.strip().rstrip("\\").strip()
            if not token:
                continue
            if token.startswith("-") or token.startswith("2>>"):
                continue          # flag, flag value, or the log redirect
            if token.endswith(".bam"):
                return token
            break
    raise AssertionError("no regtools junctions extract invocation found")


@pytest.fixture(scope="module")
def rendered_off(tmp_path_factory):
    return _render(tmp_path_factory.mktemp("uniq_off"), enabled=False)


@pytest.fixture(scope="module")
def rendered_on(tmp_path_factory):
    return _render(tmp_path_factory.mktemp("uniq_on"), enabled=True)


# ── default off: the filter must be a no-op ──────────────────────────────────

def test_off_declares_no_filtered_bam(rendered_off):
    assert "nh_unique" not in rendered_off


def test_off_runs_no_prefilter(rendered_off):
    assert "samtools view" not in rendered_off
    assert "--expr" not in rendered_off


def test_off_regtools_reads_unfiltered_bam(rendered_off):
    assert regtools_input_bam(rendered_off) == _BAM


# ── opt-in on: the filter must be fully wired ────────────────────────────────

def test_on_declares_filtered_bam_and_index(rendered_on):
    assert _FILTERED_BAM in rendered_on
    assert f"{_FILTERED_BAM}.bai" in rendered_on


def test_on_preflights_samtools_expression_support(rendered_on):
    assert "--expr" in rendered_on
    assert "exit 1" in rendered_on


def test_on_prefilters_on_nh_tag(rendered_on):
    assert "-e '[NH]==1'" in rendered_on


def test_on_regtools_reads_the_filtered_bam(rendered_on):
    """The whole point: junctions come from the filtered alignments, not the raw BAM."""
    assert regtools_input_bam(rendered_on) == _FILTERED_BAM


def test_on_still_produces_the_unfiltered_bam(rendered_on):
    """The raw BAM stays a declared output; downstream rules and QC still want it."""
    assert f"-o {_BAM} -" in rendered_on or f"{_BAM}," in rendered_on
