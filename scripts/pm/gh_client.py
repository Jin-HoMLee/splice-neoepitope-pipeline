#!/usr/bin/env python3
"""Shared hardened ``gh`` subprocess wrapper for the board/hygiene tooling (Issue #1017).

One home for the ``gh()`` helper that every board script re-rolled by hand. It bakes
in the retry+backoff logic first proven in ``recheck_milestone.py`` (Issue #711 D4) so
the same three-mode failure family is handled once, not re-patched per script.

The three failure modes every ``gh`` call site faces (Issue #1017 body):

- **(a) transient transport failure** (5xx / secondary-rate-limit 403 / replication
  lag / network) -> a non-zero exit unrelated to the request. Retried here with
  exponential backoff, honoring a ``Retry-After`` hint; a terminal failure raises
  :class:`GhError`.
- **(b) data-shape error routed through ``--jq``** -> jq exits non-zero and is
  conflated with (a), so an over-broad ``except`` silently drops good data (a
  deleted-account null ``.user`` once dropped *every* comment on an issue, #1011).
  **House rule: keep ``jq`` out of ``gh``.** Do NOT pass ``--jq`` to this wrapper;
  fetch raw JSON (``parse_json=True``, the default) and filter in Python, or pipe to
  a separate ``jq`` process. Then a data-shape error can never masquerade as a
  transport failure and be retried/swallowed.
- **(c) wrong-but-successful response** -> no exception at all, a silent wrong answer.
  The known instance: ``gh search "is:open is:blocked"`` returns 0 for genuinely
  blocked issues (#745). **Caveat: prefer a bare ``is:blocked`` search or the
  GraphQL ``blockedBy`` connection**; this wrapper cannot detect a 0-exit wrong body.

Callers isolate a per-item failure by catching :class:`GhError` (or the base
``subprocess.CalledProcessError`` it subclasses) around the per-item ``gh`` call, so
one item's terminal failure skips that item instead of aborting the whole sweep
(the #989/#1012 shape).

Import from a sibling ``scripts/pm`` script via the established pattern::

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from gh_client import gh, GhError
"""
from __future__ import annotations

import json
import re
import subprocess
import time

GH_MAX_ATTEMPTS = 4
assert GH_MAX_ATTEMPTS >= 1, "gh() runs the loop at least once; a terminal raise needs a result"
GH_BACKOFF_BASE_SECONDS = 2.0
# Clamp a Retry-After hint: GitHub can legally emit a large value (e.g. 3600s)
# during sustained degradation, which would otherwise stall a live job for that
# whole duration before raising (Issue #711 review).
GH_RETRY_AFTER_CAP_SECONDS = 120.0

# Deterministic client errors that won't fix themselves on retry - a bad request
# stays bad. Everything else (5xx, 403/secondary-rate-limit, network/timeout, an
# empty-stderr crash) is treated as transient and retried. We fail TOWARD retrying:
# over-retrying a genuine bad arg only wastes a bounded few seconds, whereas
# under-retrying a transient blip reds the whole job (Issue #711 D4).
_DETERMINISTIC_HTTP_RE = re.compile(r"HTTP (400|401|404|410|422)\b")
_RETRY_AFTER_RE = re.compile(r"retry[- ]after[:\s]+(\d+)", re.IGNORECASE)


class GhError(subprocess.CalledProcessError):
    """A terminal ``gh`` failure after retries are exhausted (Issue #1017).

    Subclasses ``subprocess.CalledProcessError`` so a caller can catch the specific
    type to isolate a per-item failure while every pre-existing
    ``except subprocess.CalledProcessError`` site keeps catching it unchanged.
    """


def _is_transient_gh_error(stderr: str) -> bool:
    return not _DETERMINISTIC_HTTP_RE.search(stderr or "")


def _retry_after_seconds(stderr: str) -> float | None:
    m = _RETRY_AFTER_RE.search(stderr or "")
    return float(m.group(1)) if m else None


def gh(*args: str, parse_json: bool = True, _runner=subprocess.run, _sleep=time.sleep) -> object:
    """Run ``gh`` with retry + exponential backoff on transient failures.

    A transient non-zero exit (mode (a)) is retried up to ``GH_MAX_ATTEMPTS`` with
    exponential backoff, honoring a ``Retry-After`` hint when GitHub provides one
    (clamped to ``GH_RETRY_AFTER_CAP_SECONDS``). A deterministic 4xx is not retried,
    and a terminal failure raises :class:`GhError` - preserving the ``check=True``
    contract callers depend on. ``_runner`` / ``_sleep`` are injection seams for tests.

    ``--jq`` (mode (b) house rule above) is **rejected with ``ValueError``**: fetch
    raw JSON and filter in Python so a data-shape error cannot be misread as a
    transport failure. This enforces the rule at the single chokepoint rather than
    leaving it a convention each call site could slip (mechanism over memory).

    ``parse_json`` defaults to True (every board read expects JSON). For a mutation
    call with empty stdout (``gh issue edit`` / ``gh pr comment``) pass
    ``parse_json=False``, else a *successful* call raises ``json.JSONDecodeError`` on
    the empty body - a confusing success-that-looks-like-failure.
    """
    # Enforce the mode-(b) house rule at the chokepoint: gh's -q/--jq routes a jq
    # data-shape error through the same non-zero exit as a transport failure, so it
    # would be retried/swallowed. Reject it before any subprocess runs.
    if any(a == "--jq" or a == "-q" or a.startswith("--jq=") for a in args):
        raise ValueError(
            "gh(): pass raw JSON and filter in Python; --jq/-q conflates a data-shape "
            "error with a transport failure (mode b, #1011)"
        )
    cmd = ["gh", *args]
    result = None
    for attempt in range(GH_MAX_ATTEMPTS):
        result = _runner(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return json.loads(result.stdout) if parse_json else result.stdout
        if not _is_transient_gh_error(result.stderr):
            break
        if attempt < GH_MAX_ATTEMPTS - 1:
            delay = _retry_after_seconds(result.stderr)
            if delay is None:
                delay = GH_BACKOFF_BASE_SECONDS * (2 ** attempt)
            _sleep(min(delay, GH_RETRY_AFTER_CAP_SECONDS))
    raise GhError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
