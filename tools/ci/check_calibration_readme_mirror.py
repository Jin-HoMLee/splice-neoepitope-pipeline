#!/usr/bin/env python3
"""Canary: the calibration README's LOCO table must mirror its notebook (Issue #1184).

`research/experiments/issue_547_immunogenicity_calibration/README.md` quotes a
per-cohort LOCO discrimination table (AUPRC / lift / AUPRG). Those numbers are
computed in `notebook.ipynb` and the README is a hand-copied mirror of them, so
it drifts silently whenever the notebook is refreshed and the README is not -
exactly what happened in Issue #1065 (the README kept quoting pre-#805 values for
weeks). PR #1183 fixed the instance and named the notebook as the source of
truth; this promotes the check PR #1183 ran by hand into a mechanism.

The notebook's committed output cell wins. This module parses the source-of-truth
table out of that output cell and the mirror table out of the README, and asserts
they agree: AUPRC and AUPRG exactly, lift within rounding (the README rounds
146.7 to ~147). A mismatch exits non-zero and names every disagreeing cell.

Pure stdlib; the pure functions are unit-tested in
`test_check_calibration_readme_mirror.py`, and a path-filtered CI job
(`.github/workflows/calibration-readme-canary.yml`) runs this end-to-end.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
EXPERIMENT_DIR = REPO_ROOT / "research" / "experiments" / "issue_547_immunogenicity_calibration"
DEFAULT_NOTEBOOK = EXPERIMENT_DIR / "notebook.ipynb"
DEFAULT_README = EXPERIMENT_DIR / "README.md"

# Lift is rounded for display (146.7 -> ~147), so it is compared within this
# absolute tolerance; AUPRC and AUPRG are asserted to the README's full precision.
LIFT_TOL = 0.5


@dataclass(frozen=True)
class LocoRow:
    auprc: float
    lift: float
    auprg: float


# --- notebook (source of truth) ----------------------------------------------


def _code_cell_output_text(cell: dict) -> str:
    """Concatenate the human-visible text of a code cell's outputs."""
    text = ""
    for out in cell.get("outputs", []):
        if out.get("output_type") == "stream":
            text += "".join(out.get("text", []))
        elif "data" in out:
            plain = out["data"].get("text/plain")
            if plain:
                text += "".join(plain)
    return text


def _is_header_row(line: str) -> bool:
    """True for the table's column-header row: cohort + AUPRC + lift + AUPRG as
    exact whitespace tokens.

    Token-exact, not substring: the cell's *title* line ("Per-cohort discrimination
    - AUPRC/prevalence lift vs AUPRG ...") contains all four words as substrings
    but none as standalone columns, so a substring test would mistake it for the
    header and parse zero rows.
    """
    toks = {t.lower() for t in line.split()}
    return {"cohort", "auprc", "lift", "auprg"} <= toks


def _is_loco_source_table(text: str) -> bool:
    """A discrimination table whose header carries cohort + AUPRC + lift + AUPRG.

    This distinguishes the source table (the `Per-cohort discrimination` cell) from
    the neighbouring `LOCO results` cell, which has auprc/auprg but no `lift`
    column - so a header row with all four is the unambiguous discriminator.
    """
    return any(_is_header_row(line) for line in text.splitlines())


def parse_notebook_loco(notebook_path: Path) -> dict[str, LocoRow]:
    """Parse the per-cohort AUPRC/lift/AUPRG table from the notebook's output cell.

    Locates the single code cell whose output holds the discrimination table
    (identified by content, not index, so a re-ordered notebook still resolves)
    and parses its whitespace-delimited rows by header-column name. Raises if zero
    or more than one such table is found - a notebook restructure should fail this
    loudly, never silently pick the wrong cell.
    """
    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    matches = [
        _code_cell_output_text(c)
        for c in nb.get("cells", [])
        if c.get("cell_type") == "code" and _is_loco_source_table(_code_cell_output_text(c))
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one LOCO discrimination table in {notebook_path.name}, "
            f"found {len(matches)} (notebook restructured? update the canary)"
        )
    return _parse_whitespace_table(matches[0])


def _parse_whitespace_table(text: str) -> dict[str, LocoRow]:
    """Parse a pandas-`to_string` table into {cohort -> LocoRow} by column name."""
    lines = text.splitlines()
    header_idx = next(i for i, ln in enumerate(lines) if _is_header_row(ln))
    header = lines[header_idx].split()
    col = {name.lower(): j for j, name in enumerate(header)}
    out: dict[str, LocoRow] = {}
    for ln in lines[header_idx + 1 :]:
        toks = ln.split()
        if len(toks) != len(header):
            continue  # blank line or trailing prose after the table body
        try:
            out[toks[col["cohort"]]] = LocoRow(
                auprc=float(toks[col["auprc"]]),
                lift=float(toks[col["lift"]]),
                auprg=float(toks[col["auprg"]]),
            )
        except (KeyError, ValueError):
            continue
    return out


# --- README (mirror) ---------------------------------------------------------

_NUM_RE = re.compile(r"[-+]?\d*\.?\d+")


def _first_float(cell: str) -> Optional[float]:
    m = _NUM_RE.search(cell)
    return float(m.group()) if m else None


def parse_readme_loco(readme_path: Path) -> dict[str, LocoRow]:
    """Parse the LOCO mirror table out of the README markdown.

    Finds the pipe table whose header row names AUPRC, lift and AUPRG, then reads
    each data row, extracting the leading float from each cell (so `~147×` and
    `~1.2× (near baseline)` both parse). Raises if the table is not found.
    """
    lines = readme_path.read_text(encoding="utf-8").splitlines()
    header_idx = None
    cols: dict[str, int] = {}
    for i, ln in enumerate(lines):
        if not ln.lstrip().startswith("|"):
            continue
        low = ln.lower()
        if "auprc" in low and "auprg" in low and "lift" in low and "cohort" in low:
            cells = [c.strip().lower() for c in ln.strip().strip("|").split("|")]
            for j, name in enumerate(cells):
                if "cohort" in name:
                    cols["cohort"] = j
                elif name == "auprc":
                    cols["auprc"] = j
                elif "lift" in name:
                    cols["lift"] = j
                elif name == "auprg":
                    cols["auprg"] = j
            header_idx = i
            break
    if header_idx is None or not {"cohort", "auprc", "lift", "auprg"} <= set(cols):
        raise ValueError(f"LOCO mirror table not found in {readme_path.name}")

    out: dict[str, LocoRow] = {}
    # skip the header separator row (|---|---|), then read until the table ends
    for ln in lines[header_idx + 2 :]:
        if not ln.lstrip().startswith("|"):
            break
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if len(cells) <= max(cols.values()):
            continue
        cohort = cells[cols["cohort"]].strip()
        auprc = _first_float(cells[cols["auprc"]])
        lift = _first_float(cells[cols["lift"]])
        auprg = _first_float(cells[cols["auprg"]])
        if cohort and None not in (auprc, lift, auprg):
            out[cohort] = LocoRow(auprc=auprc, lift=lift, auprg=auprg)
    return out


# --- comparison --------------------------------------------------------------


def compare(notebook: dict[str, LocoRow], readme: dict[str, LocoRow]) -> list[str]:
    """Return a list of human-readable mismatches; empty means the mirror is clean.

    AUPRC and AUPRG must match the notebook's value exactly (the README quotes them
    to full precision, e.g. `0.0322`, so a lower-precision `0.032` is genuine drift,
    not a faithful rounding - this is the exact class PR #1183 fixed on IMPROVE).
    Lift alone is compared within `LIFT_TOL`, because the README deliberately rounds
    it for display (146.7 -> ~147).
    """
    problems: list[str] = []
    nb_cohorts, rm_cohorts = set(notebook), set(readme)
    if nb_cohorts != rm_cohorts:
        missing = nb_cohorts - rm_cohorts
        extra = rm_cohorts - nb_cohorts
        if missing:
            problems.append(f"cohorts in notebook but not README: {sorted(missing)}")
        if extra:
            problems.append(f"cohorts in README but not notebook: {sorted(extra)}")

    for cohort in sorted(nb_cohorts & rm_cohorts):
        nb, rm = notebook[cohort], readme[cohort]
        for field in ("auprc", "auprg"):
            nb_v, rm_v = getattr(nb, field), getattr(rm, field)
            if nb_v != rm_v:
                problems.append(
                    f"{cohort} {field.upper()}: README {rm_v} != notebook {nb_v}"
                )
        if abs(nb.lift - rm.lift) > LIFT_TOL:
            problems.append(
                f"{cohort} lift: README {rm.lift} vs notebook {nb.lift} "
                f"(exceeds rounding tolerance {LIFT_TOL})"
            )
    return problems


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--notebook", type=Path, default=DEFAULT_NOTEBOOK)
    ap.add_argument("--readme", type=Path, default=DEFAULT_README)
    args = ap.parse_args(argv)

    notebook = parse_notebook_loco(args.notebook)
    readme = parse_readme_loco(args.readme)

    # Coverage floor (mirrors the annotate-flag-canary's "agreement AND coverage"
    # convention): compare() returns [] on two empty dicts, so a run where BOTH
    # sides parse to zero cohorts - a header that resolves but every data row is
    # dropped (a leading pandas index column, a cohort label with a space) on both
    # notebook and README at once - would pass vacuously. Refuse a zero-cohort
    # parse outright so the canary can only be green on a real, populated match.
    if not notebook or not readme:
        empty = "notebook" if not notebook else "README"
        print(
            f"Coverage floor: parsed 0 cohorts from the {empty} - the mirror cannot "
            "be verified. A zero-cohort parse must never pass vacuously (a table "
            "header resolved but no data rows parsed)."
        )
        return 1

    problems = compare(notebook, readme)
    if problems:
        print("Calibration README LOCO table drifted from the notebook (source of truth):")
        for p in problems:
            print(f"  - {p}")
        print("\nFix: refresh the README table from notebook.ipynb (the notebook wins).")
        return 1
    print(f"Calibration README mirrors the notebook: {len(readme)} cohorts agree.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
