"""Shared skip guard for live ``gh`` Projects-v2 integration tests in ``tools/ci``.

The default repo-scoped ``GITHUB_TOKEN`` cannot read Projects v2, so live tests
that resolve project items (or hit any Projects-v2 API) must *skip*, not
*error*, when that capability is absent. That happens on fork PRs where the
``GH_PROJECT_TOKEN`` secret isn't exposed, and in local environments without a
``read:project``-scoped ``gh`` login.

``test_recheck_dispatch.py`` and ``test_recheck_milestone.py`` both apply the
``REQUIRES_LIVE_GH`` marker exported here, so the probe runs once and the two
live suites share a single graceful-skip contract.

Not a ``test_*`` module, so pytest does not collect it as tests.
"""
import shutil
import subprocess

import pytest


def _gh_has_project_read_scope() -> bool:
    """Probe whether the current ``gh`` auth can read Projects v2 (one GraphQL call).

    Returns ``False`` (→ skip) when ``gh`` is missing, the call times out, or the
    query fails (e.g. the token lacks ``read:project``). Never raises, so a
    missing capability skips the dependent test instead of erroring it.
    """
    if not shutil.which("gh"):
        return False
    try:
        result = subprocess.run(
            ["gh", "api", "graphql", "-f",
             "query=query { viewer { projectsV2(first: 1) { totalCount } } }"],
            capture_output=True, timeout=10, check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


REQUIRES_LIVE_GH = pytest.mark.skipif(
    not _gh_has_project_read_scope(),
    reason="requires gh auth with project read scope (set GH_TOKEN to a PAT with `read:project`)",
)
