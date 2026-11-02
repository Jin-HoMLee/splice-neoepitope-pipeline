from pathlib import Path

import pytest

from app.command_builder import build_snakemake_argv, display_command


def test_configfiles_share_one_flag_and_target_follows_separator():
    argv = build_snakemake_argv(
        configfiles=[Path("config/config.yaml"), Path("overlay.yaml")],
        target=Path("results/P1/reports/report.html"),
        cores=4,
    )
    assert argv.count("--configfile") == 1
    config_index = argv.index("--configfile")
    separator = argv.index("--")
    assert argv[config_index + 1:separator] == ["config/config.yaml", "overlay.yaml"]
    assert argv[separator + 1:] == ["results/P1/reports/report.html"]
    assert "--rerun-incomplete" in argv


def test_dry_run_uses_identical_safe_shape():
    argv = build_snakemake_argv(
        configfiles=[Path("base.yaml")], target=Path("report.html"), cores=2, dry_run=True
    )
    assert "--dry-run" in argv
    assert "--printshellcmds" in argv
    assert display_command(argv).startswith("snakemake --cores 2")


def test_missing_config_or_unbounded_cores_fail():
    with pytest.raises(ValueError, match="config"):
        build_snakemake_argv(configfiles=[], target=Path("x"), cores=1)
    with pytest.raises(ValueError, match="cores"):
        build_snakemake_argv(configfiles=[Path("x")], target=Path("y"), cores=9)
