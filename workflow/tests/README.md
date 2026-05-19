# Pipeline Tests

Unit + snapshot tests for the pipeline's Python helpers and Snakemake rule shell commands.

## Test venv setup (one-time per clone)

The test runner is a pyenv-managed Python venv separate from the workflow's `snakemake` conda env. The two envs cooperate at run time: the venv provides `pytest` + test deps, and the activated conda env puts the `snakemake` binary on PATH for tests that invoke it via subprocess (e.g. [test_alignment_star_command.py](test_alignment_star_command.py)).

```bash
# 1. Pin the Python version used by this clone (writes .python-version, gitignored)
pyenv local 3.14.4

# 2. Create the test venv inside workflow/tests/ (also gitignored)
python -m venv workflow/tests/.venv

# 3. Install test dependencies
workflow/tests/.venv/bin/pip install --upgrade pip
workflow/tests/.venv/bin/pip install -r workflow/tests/requirements-test.txt
```

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

## Why this split

- **`snakemake` conda env = workflow runner.** Snakemake + solver. Keep pristine; adding test deps risks dep drift in the workhorse env.
- **`.venv` = test runner.** Pinned via [requirements-test.txt](requirements-test.txt). Fast to recreate. Isolated.
- **CI parity.** CI runs tests via pip in a clean venv, not a conda env — local mirrors that.

The migration from git worktrees to separate clones (2026-05-14) means each clone needs its own one-time `.venv` setup; no sharing across clones.
