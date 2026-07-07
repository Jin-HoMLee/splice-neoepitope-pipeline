"""Snapshot test of the rendered hisat2_align read-input args.

Complements the pure-helper test (test_hisat2_command.py) at the rule layer:
it asserts the `hisat2_align` shell command actually wires `{params.read_args}`
so a single-end sample renders `-U <r1>` and a paired-end sample renders
`-1 <r1> -2 <r2>`. The pure test guards the helper's logic; this guards that
the rule still calls it (a future edit dropping `{params.read_args}` from the
shell would pass the unit test but fail here). Sibling of
test_alignment_star_command.py. See Issue #962.

Runs `snakemake --dry-run --printshellcmds` programmatically with an
aligner=hisat2 override and stub inputs; skips where snakemake is unavailable.
"""
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def hisat2_dry_run_output(tmp_path_factory):
    """Render hisat2_align for one SE and one PE sample; return captured stdout.

    Both targets are rendered in a single dry-run so the SE (`-U`) and PE
    (`-1/-2`) commands appear in one output; assertions key on the distinct
    stub paths so the two samples can't be confused.
    """
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not on PATH - activate the snakemake conda env")

    tmp_path = tmp_path_factory.mktemp("hisat2_stub")
    se_r1 = tmp_path / "se_R1.fq.gz"
    pe_r1 = tmp_path / "pe_R1.fq.gz"
    pe_r2 = tmp_path / "pe_R2.fq.gz"
    for f in (se_r1, pe_r1, pe_r2):
        f.touch()

    index_dir = tmp_path / "hisat2_index"
    index_dir.mkdir()
    (index_dir / "index.done").touch()
    (index_dir / "genome.1.ht2").touch()

    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "patient_id\tsample_id\tsample_type\tfastq1\tfastq2\tstrandness\n"
        f"tp\tse_sample\tSolid Tissue Normal\t{se_r1}\t\t\n"
        f"tp\tpe_sample\tPrimary Tumor\t{pe_r1}\t{pe_r2}\treverse\n"
    )

    override = tmp_path / "override.yaml"
    override.write_text(textwrap.dedent(f"""\
        samples_tsv: "{samples}"
        alignment:
          aligner: "hisat2"
          hisat2_index_dir: "{index_dir}"
          hisat2_prebuilt_url: ""
    """))

    result = subprocess.run(
        [
            "snakemake", "-n", "--printshellcmds",
            "--configfile", "config/config.yaml", str(override),
            "--",
            "results/tp/alignment/se_sample/raw_junctions.tsv",
            "results/tp/alignment/pe_sample/raw_junctions.tsv",
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
    return {"out": result.stdout, "se_r1": str(se_r1),
            "pe_r1": str(pe_r1), "pe_r2": str(pe_r2)}


def test_single_end_sample_uses_U(hisat2_dry_run_output):
    o = hisat2_dry_run_output
    assert f"-U {o['se_r1']}" in o["out"]


def test_single_end_sample_not_rendered_paired(hisat2_dry_run_output):
    o = hisat2_dry_run_output
    assert f"-1 {o['se_r1']}" not in o["out"]


def test_paired_end_sample_uses_1_2(hisat2_dry_run_output):
    o = hisat2_dry_run_output
    assert f"-1 {o['pe_r1']} -2 {o['pe_r2']}" in o["out"]


def test_paired_end_sample_not_rendered_single(hisat2_dry_run_output):
    o = hisat2_dry_run_output
    assert f"-U {o['pe_r1']}" not in o["out"]
