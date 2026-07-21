"""ASNEO runner (#966, AC-1).

Invokes the open-only (option-B patched) ASNEO caller on a chr22-scale
SJ.out.tab in the ``asneo`` conda env. The clone / patch / genome / run
sequence lives in the proven ``run_asneo.sh`` (adapted from the #965 smoke);
this module shells out to it and validates the result.

Runner contract: ``run_asneo(sj_tab, workdir) -> (output_path, provenance)``.
The provenance dict is the only handle tying ASNEO's peptide-level (null
junction) records back to their origin, so it carries the input SJ.out.tab
name + a content digest + the thresholds used. ``DEFAULT_THRESHOLDS`` is the
single source for those thresholds - passed to the bash script AND recorded in
provenance, so the two cannot drift.
"""

import hashlib
import subprocess
from pathlib import Path

_SCRIPT = Path(__file__).with_name("run_asneo.sh")

DEFAULT_THRESHOLDS = {"reads": 2, "psi": 0.05, "lengths": "8,9,10,11"}


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


def _sha1(path: str) -> str:
    return hashlib.sha1(Path(path).read_bytes()).hexdigest()


def run_asneo(sj_tab: str, workdir: str, thresholds: dict = None):
    """Run patched ASNEO on ``sj_tab``; return ``(peptide_path, provenance)``.

    ``provenance`` records the input basename + content SHA-1 + the thresholds
    used, so the (peptide-level, null-junction) records stay traceable to their
    origin - the recoverability the design spec asks for pending #1258.
    """
    t = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    subprocess.run(
        ["bash", str(_SCRIPT), sj_tab, workdir,
         str(t["reads"]), str(t["psi"]), str(t["lengths"])],
        check=True,
    )
    out = peptide_output_path(workdir)
    validate_output(out)
    provenance = {
        "sj_tab": Path(sj_tab).name,
        "sj_tab_sha1": _sha1(sj_tab),
        "thresholds": dict(t),
    }
    return out, provenance
