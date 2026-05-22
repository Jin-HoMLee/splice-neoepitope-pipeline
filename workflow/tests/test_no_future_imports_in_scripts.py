"""Regression: workflow/scripts/*.py must not use `from __future__`.

Snakemake's `PythonScript.write_script` prepends an executable preamble
(`snakemake = pickle.loads(...)`) before the script source — see
`snakemake/script/__init__.py` (snakemake 8.30, line ~807). Any `__future__`
import in the source then lands *after* executable code, raising:

    SyntaxError: from __future__ imports must occur at the beginning of the file

The `PY_PREAMBLE_RE` regex at line 44 is defined but not wired up
(see the TODO comment in upstream Snakemake source). Until that's fixed
upstream, scripts invoked via `script:` directive must not use `from
__future__`.

This test simulates Snakemake's wrapping and asserts every workflow script
parses cleanly. Caught latently in Issue #461 — PR #428 introduced future
imports across 6 scripts to support lazy pandas imports; chr22 integration
runs targeted only `panel.tsv`, so the breakage never executed.
"""
import pathlib
import py_compile
import tempfile

import pytest

SCRIPTS_DIR = pathlib.Path(__file__).resolve().parent.parent / "scripts"

# Minimal stand-in for the Snakemake preamble. Only property that matters
# for this test is "executable code precedes the source" — the upstream
# preamble is one long minified line; the trigger is identical.
SNAKEMAKE_PREAMBLE = (
    "######## snakemake preamble start (automatically inserted, do not edit) ########\n"
    "import sys, pickle; snakemake = None  # stubbed for parse-check\n"
    "######## snakemake preamble end #########\n"
)


@pytest.mark.parametrize(
    "script_path",
    sorted(SCRIPTS_DIR.glob("*.py")),
    ids=lambda p: p.name,
)
def test_script_parses_under_snakemake_wrapper(script_path: pathlib.Path) -> None:
    source = script_path.read_text()
    wrapped = SNAKEMAKE_PREAMBLE + source
    with tempfile.NamedTemporaryFile(
        suffix=".py", mode="w", delete=False
    ) as f:
        f.write(wrapped)
        tmp_path = f.name
    try:
        py_compile.compile(tmp_path, doraise=True)
    except py_compile.PyCompileError as e:
        pytest.fail(
            f"{script_path.name} fails to parse under Snakemake's preamble "
            f"injection. Likely cause: `from __future__ import annotations` "
            f"at the top of the file (Snakemake's preamble runs before it, "
            f"making the future import non-leading). See Issue #461.\n"
            f"\nOriginal compile error:\n{e}"
        )
