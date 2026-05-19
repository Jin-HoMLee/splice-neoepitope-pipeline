"""Snapshot test of the rendered star_align shell command.

Asserts the production STAR alignment command contains the flags required
for novel-junction sensitivity (--twopassMode Basic, --limitSjdbInsertNsj)
and does not silently throttle multi-mappers
(--outSJfilterReads Unique, --outSJfilterCountUniqueMin, --outSJfilterCountTotalMin).

The test runs `snakemake --dry-run --printshellcmds` programmatically with an
aligner=star config override and stub inputs in a tmp_path. The captured
shell output is asserted against substring presence/absence.

This is the first shell-level rule snapshot test in the codebase. Existing
tests cover pure Python helpers (test_strandness.py, test_bed12_to_junctions.py,
test_star_sj_to_junctions.py, etc.); this complements them at the rule layer.
"""
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def star_dry_run_output(tmp_path):
    """Render the star_align shell command and return captured stdout."""
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not on PATH — activate the snakemake conda env")

    # Stub FASTQs (must exist for `ancient(_get_fastq*)` to skip producer-rule lookup)
    fq1 = tmp_path / "test_R1.fq.gz"
    fq2 = tmp_path / "test_R2.fq.gz"
    fq1.touch()
    fq2.touch()

    # Stub STAR index (sentinel + dir must exist for star_align inputs to resolve)
    star_index_dir = tmp_path / "star_index"
    star_index_dir.mkdir()
    (star_index_dir / "index.done").touch()

    # Stub samples.tsv
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "patient_id\tsample_id\tsample_type\tfastq1\tfastq2\n"
        f"testpatient\ttestsample\tPrimary Tumor\t{fq1}\t{fq2}\n"
    )

    # Config override: aligner=star, point at stub paths
    override = tmp_path / "override.yaml"
    override.write_text(textwrap.dedent(f"""\
        samples_tsv: "{samples}"
        alignment:
          aligner: "star"
          star_index_dir: "{star_index_dir}"
    """))

    target = "results/testpatient/alignment/testsample/junctions.tsv"

    result = subprocess.run(
        [
            "snakemake",
            "-n",
            "--printshellcmds",
            "--configfile", "config/config.yaml", str(override),
            "--",
            target,
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.fail(
            f"snakemake dry-run failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )
    return result.stdout


def test_star_align_uses_twopass_mode(star_dry_run_output):
    """2-pass mode is the flag that gives STAR its novel-junction edge over HISAT2."""
    assert "--twopassMode Basic" in star_dry_run_output


def test_star_align_raises_sjdb_insert_limit(star_dry_run_output):
    """Default 1M cap can abort runs when 2-pass + multi-mappers discover more."""
    assert "--limitSjdbInsertNsj 2000000" in star_dry_run_output


def test_star_align_does_not_filter_to_unique_reads_only(star_dry_run_output):
    """Paper keeps multi-mappers; --outSJfilterReads Unique silently drops them."""
    assert "--outSJfilterReads Unique" not in star_dry_run_output


def test_star_align_does_not_throttle_min_unique_count(star_dry_run_output):
    """--outSJfilterCountUniqueMin throttles output below paper baseline."""
    assert "--outSJfilterCountUniqueMin" not in star_dry_run_output


def test_star_align_does_not_throttle_min_total_count(star_dry_run_output):
    """--outSJfilterCountTotalMin throttles output below paper baseline."""
    assert "--outSJfilterCountTotalMin" not in star_dry_run_output
