# Pipeline Tests

Unit + snapshot tests for the pipeline's Python helpers and Snakemake rule shell commands.

## Test venv setup (one-time per clone)

The test runner is a pyenv-managed Python venv separate from the workflow's `snakemake` conda env. The two envs cooperate at run time: the venv provides `pytest` + test deps, and the activated conda env puts the `snakemake` binary on PATH for tests that invoke it via subprocess (e.g. [test_alignment_star_command.py](test_alignment_star_command.py)).

Prerequisite: [`uv`](https://docs.astral.sh/uv/) (Astral's fast installer/resolver) - `brew install uv` or the standalone installer `curl -LsSf https://astral.sh/uv/install.sh | sh`.

```bash
# 1. Pin the Python version used by this clone (writes .python-version, gitignored).
#    Any stable Python >=3.10 works; 3.13.5 is recommended (recent, widely
#    available via pyenv on fresh clones). Test deps in requirements-test.txt
#    have no upper Python bound.
pyenv local 3.13.5

# 2. Create the test venv inside workflow/tests/ (also gitignored) with uv.
uv venv --python 3.13.5 workflow/tests/.venv

# 3. Install test dependencies with uv (10-100x faster than pip; no separate
#    pip-upgrade step needed).
uv pip install --python workflow/tests/.venv/bin/python -r workflow/tests/requirements-test.txt
```

Only the two per-clone pyenv venvs move to `uv`. The `snakemake` conda env and the per-rule `--use-conda` envs (`workflow/envs/*.yaml`) are unchanged - `uv` operates on the pip/PyPI side and does not touch conda's binary solver.

## Running tests

```bash
# Tests that don't invoke snakemake: venv alone is enough
workflow/tests/.venv/bin/python -m pytest workflow/tests/ -v

# Tests that invoke snakemake as a subprocess (e.g. the star_align snapshot
# test) need snakemake on PATH — activate the conda env first
conda activate snakemake
workflow/tests/.venv/bin/python -m pytest workflow/tests/ -v
```

The fixture in [test_alignment_star_command.py](test_alignment_star_command.py) calls `pytest.skip` if `snakemake` is not on PATH, so the test suite still runs cleanly without the conda env — those specific tests just skip.

## Running on memory-tight machines (e.g. M1 8 GB)

When the OS page cache is under pressure, `pytest --collect-only` can stall for minutes — fresh `import pandas` and pytest's assertion-rewriter pass both turn into disk-bound page faults. Two workarounds while the suite stays pandas-heavy at module scope (see [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364)):

```bash
# 1. Run only the touched test file(s); let CI cover the full suite.
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_filter_junctions.py -q

# 2. Disable assertion rewriting (faster collection, less helpful assert errors).
workflow/tests/.venv/bin/python -m pytest workflow/tests/ -p no:assertion -q
```

To re-measure the static import cost surface (per-test-file and per-script cold module-load):

```bash
workflow/tests/.venv/bin/python workflow/tests/scripts/profile_imports.py
```

Last measurement (non-RAM-pressured run): full collection of 287 tests in ~20s; per-test cold module-load 0.03–0.62s, dominated by `pandas` (0.26s cold) imported at module scope in 6 test files and 6 project scripts.

## Why this split

- **`snakemake` conda env = workflow runner.** Snakemake + solver. Keep pristine; adding test deps risks dep drift in the workhorse env.
- **`.venv` = test runner.** Pinned via [requirements-test.txt](requirements-test.txt). Fast to recreate. Isolated.
- **CI parity.** CI runs tests via pip in a clean venv, not a conda env — local mirrors that.

The migration from git worktrees to separate clones (2026-05-14) means each clone needs its own one-time `.venv` setup; no sharing across clones.
