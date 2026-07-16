"""Unit tests for the calibration README<->notebook mirror canary (Issue #1184).

The load-bearing test is `test_stale_pre_1183_readme_is_red`: a check that cannot
fail on the actual historical drift (Issue #1065, the pre-#805 numbers surviving
in the README) is not a check. It reconstructs that stale table and asserts the
canary flags every drifted cell.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_calibration_readme_mirror as m  # noqa: E402

# The pre-#1183 (pre-#805) values that silently survived in the README while the
# notebook already carried the refreshed ones - straight from the PR #1183 body.
_STALE_README = """# Calibration

| cohort left out | prevalence | AUPRC | LOCO lift | AUPRG |
|---|---|---|---|---|
| NCI | 0.000245 | 0.027 | ~111x | 0.998 |
| TESLA | 0.046 | 0.201 | ~4.4x | 0.859 |
| HiTIDE | 0.026 | 0.071 | ~2.7x | 0.815 |
| IMPROVE | 0.027 | 0.032 | ~1.2x | 0.268 |
"""

# A synthetic notebook whose ONLY discrimination-table cell leads with the exact
# prose title that broke the first parser (its substrings look like a header).
_TITLE_TRAP_NOTEBOOK = {
    "cells": [
        {
            "cell_type": "code",
            "outputs": [
                {
                    "output_type": "stream",
                    "text": [
                        "Per-cohort discrimination - AUPRC/prevalence lift vs AUPRG (comparable):\n",
                        " cohort  true_pos  prevalence  AUPRC  lift  AUPRG\n",
                        "    NCI       103    0.000245 0.0359 146.7 0.9999\n",
                        "  TESLA        34    0.046196 0.2482   5.4 0.8779\n",
                    ],
                }
            ],
        }
    ],
}


def test_parse_notebook_loco_matches_source_of_truth():
    nb = m.parse_notebook_loco(m.DEFAULT_NOTEBOOK)
    assert set(nb) == {"NCI", "TESLA", "HiTIDE", "IMPROVE"}
    assert nb["NCI"] == m.LocoRow(auprc=0.0359, lift=146.7, auprg=0.9999)
    assert nb["IMPROVE"] == m.LocoRow(auprc=0.0322, lift=1.2, auprg=0.2678)


def test_parse_readme_loco_reads_the_mirror_table():
    rm = m.parse_readme_loco(m.DEFAULT_README)
    assert set(rm) == {"NCI", "TESLA", "HiTIDE", "IMPROVE"}
    # `~147x` and `~1.2x (near baseline)` both parse to their leading float
    assert rm["NCI"].lift == 147.0
    assert rm["IMPROVE"].auprg == 0.2678


def test_current_readme_mirrors_the_notebook_green():
    nb = m.parse_notebook_loco(m.DEFAULT_NOTEBOOK)
    rm = m.parse_readme_loco(m.DEFAULT_README)
    assert m.compare(nb, rm) == []


def test_stale_pre_1183_readme_is_red(tmp_path):
    """The canary MUST flag the real historical drift - the AC2 falsifier."""
    stale = tmp_path / "README.md"
    stale.write_text(_STALE_README, encoding="utf-8")
    nb = m.parse_notebook_loco(m.DEFAULT_NOTEBOOK)
    rm = m.parse_readme_loco(stale)
    problems = m.compare(nb, rm)
    assert problems, "the canary did not flag the pre-#1183 drift - it is a hollow check"
    # every cohort's AUPRC and AUPRG drifted, so each is named
    joined = "\n".join(problems)
    for cohort in ("NCI", "TESLA", "HiTIDE", "IMPROVE"):
        assert f"{cohort} AUPRC" in joined, f"missed the {cohort} AUPRC drift"


def test_main_exit_codes(tmp_path, capsys):
    # green on the real pair
    assert m.main([]) == 0
    # red on a stale README
    stale = tmp_path / "README.md"
    stale.write_text(_STALE_README, encoding="utf-8")
    assert m.main(["--readme", str(stale)]) == 1


def test_lift_compared_within_rounding_but_auprc_exact():
    nb = {"X": m.LocoRow(auprc=0.0359, lift=146.7, auprg=0.9999)}
    # README rounds 146.7 -> 147: within tolerance, AUPRC/AUPRG exact -> clean
    assert m.compare(nb, {"X": m.LocoRow(auprc=0.0359, lift=147.0, auprg=0.9999)}) == []
    # a lift off by more than the rounding tolerance is caught
    assert m.compare(nb, {"X": m.LocoRow(auprc=0.0359, lift=150.0, auprg=0.9999)})
    # an AUPRC that differs in the last printed digit is caught (no lift excuse)
    assert m.compare(nb, {"X": m.LocoRow(auprc=0.0360, lift=146.7, auprg=0.9999)})


def test_title_line_is_not_mistaken_for_the_header():
    """Regression: the cell's prose title has cohort/AUPRC/lift/AUPRG as substrings;
    token-exact header matching must skip it and parse the real rows."""
    import json

    nb_path = None
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".ipynb", delete=False) as f:
        json.dump(_TITLE_TRAP_NOTEBOOK, f)
        nb_path = Path(f.name)
    try:
        rows = m.parse_notebook_loco(nb_path)
        assert rows["NCI"] == m.LocoRow(auprc=0.0359, lift=146.7, auprg=0.9999)
        assert rows["TESLA"].auprg == 0.8779
    finally:
        nb_path.unlink()


def test_missing_table_raises():
    import json
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".ipynb", delete=False) as f:
        json.dump({"cells": []}, f)
        p = Path(f.name)
    try:
        with pytest.raises(ValueError):
            m.parse_notebook_loco(p)
    finally:
        p.unlink()
