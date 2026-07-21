"""ASNEO runner (#966, AC-1).

Invokes the open-only (option-B patched) ASNEO caller on a chr22-scale
SJ.out.tab in the ``asneo`` conda env, returning the path to the produced
putative_peptide.txt. The clone / patch / genome / run sequence lives in the
proven ``run_asneo.sh`` (adapted from the #965 smoke); this module shells out
to it and validates the result. Command orchestration is bash; the output
contract is checked here in pure, unit-tested Python.
"""

import subprocess
from pathlib import Path

_SCRIPT = Path(__file__).with_name("run_asneo.sh")


def peptide_output_path(workdir: str) -> str:
    """Where run_asneo.sh writes the caller's candidate peptides (pure)."""
    return str(Path(workdir) / "out" / "putative_peptide.txt")


def validate_output(path: str) -> None:
    """Fail loud if ASNEO produced no candidate peptides (pure)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"ASNEO produced no output at {path}")
    peptides = [ln.strip() for ln in p.read_text().splitlines()
                if ln.strip() and not ln.strip().startswith("#")]
    if not peptides:
        raise ValueError(f"ASNEO output {path} has 0 candidate peptides")


def run_asneo(sj_tab: str, workdir: str) -> str:
    """Run patched ASNEO on ``sj_tab``; return the putative_peptide.txt path."""
    subprocess.run(["bash", str(_SCRIPT), sj_tab, workdir], check=True)
    out = peptide_output_path(workdir)
    validate_output(out)
    return out
