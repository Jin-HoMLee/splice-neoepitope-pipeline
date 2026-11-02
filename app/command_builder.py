"""Safe Snakemake argv construction. No shell strings are executed here."""

import shlex
from pathlib import Path
from typing import Iterable, List, Sequence


def build_snakemake_argv(
    *,
    configfiles: Sequence[Path],
    target: Path,
    cores: int,
    snakemake_bin: str = "snakemake",
    dry_run: bool = False,
) -> List[str]:
    if not configfiles:
        raise ValueError("at least one config file is required")
    if cores < 1 or cores > 8:
        raise ValueError("cores must be between 1 and 8")
    argv = [
        snakemake_bin,
        "--cores",
        str(cores),
        "--use-conda",
        "--rerun-triggers",
        "mtime",
        "--rerun-incomplete",
    ]
    if dry_run:
        argv.extend(["--dry-run", "--printshellcmds"])
    argv.append("--configfile")
    argv.extend(str(path) for path in configfiles)
    argv.extend(["--", str(target)])
    return argv


def display_command(argv: Iterable[str]) -> str:
    return shlex.join(list(argv))
