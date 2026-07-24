#!/usr/bin/env python3
"""Board-item lookup that refuses to hand back a result it cannot vouch for (Issue #1151).

Two roles independently shipped the same silent no-op against board 9 within 24
hours, because **every** failure mode of the obvious lookup returns empty with
exit 0: the key does not exist, the item is not carded, or the read did not cover
everything. All three are indistinguishable at the call site, so a resolver that
returns ``None`` hands a write a value nobody has vouched for.

This module resolves ``issue -> board-item id`` and **raises** instead, with three
distinguishable errors so a caller can tell "not carded" from "I could not tell".

Why ``Issue.projectItems`` and not a board scan
-----------------------------------------------
``ProjectV2.items`` defaults to ``archivedStates: [NOT_ARCHIVED]`` **and** its
``totalCount`` is scoped to that same subset. A board scan therefore omits every
archived card *and* passes a ``seen == totalCount`` completeness assertion while
doing it. Measured on board 9 (2026-07-23): **289 unarchived, 1,163 archived**.
The first draft of this module scanned the board and reported

    BoardItemNotFound: issue #569 is not on project Jin-HoMLee/9 (complete read of 289 items)

for a card that plainly exists - the exact defect this module was written to
prevent, certified sound by its own guard.

``Issue.projectItems`` defaults ``includeArchived: true`` and is scoped to one
issue, so there is no board-wide read to truncate and no archived blind spot.
Both failure modes become *unrepresentable* rather than asserted-against, which
is the stronger fix. It is also one API call instead of three.

The repo is part of the query, not decoration: board 9 spans two repos with
independent numbering, and a complete read (1,453 items, archived included)
finds **56** numbers that exist as an Issue in *both* repos - #37 among them. A
number-only match returns whichever card page-order yields first.

Do not evidence that collision by scanning board content numbers: issues and
PRs share one number space, so such a scan overcounts. #183 looks like a
collision and is not one - it is an Issue in the pipeline repo and a PULL
REQUEST in personas.
"""

import subprocess

OWNER = "Jin-HoMLee"
PROJECT_NUMBER = 9
DEFAULT_REPO = "Jin-HoMLee/splice-neoepitope-pipeline"

# includeArchived is explicit rather than defaulted: it is load-bearing, and a
# reader must not have to know the schema default to see that archived cards are
# in scope. first:20 comfortably covers an issue's project memberships; the
# hasNextPage guard below turns "more than that" into a loud error, not a wrong
# answer.
# `kind` is interpolated into the query as a field name, so it is allow-listed
# rather than trusted. Both content types expose projectItems, so one resolver
# serves an Issue and a PR's own card.
KINDS = ("issue", "pullRequest")

QUERY_TEMPLATE = """
query($owner: String!, $name: String!, $number: Int!) {
  repository(owner: $owner, name: $name) {
    %s(number: $number) {
      projectItems(first: 20, includeArchived: true) {
        pageInfo { hasNextPage endCursor }
        nodes {
          id
          isArchived
          project { number }
          fieldValues(first: 20) {
            nodes {
              ... on ProjectV2ItemFieldSingleSelectValue {
                name
                field { ... on ProjectV2FieldCommon { name } }
              }
            }
          }
        }
      }
    }
  }
}
"""


class BoardLookupError(Exception):
    """Base: the lookup cannot vouch for its answer, so it refuses to return one."""


class BoardItemNotFound(BoardLookupError):
    """A COMPLETE read showed the issue has no card on this board.

    The only trustworthy "no". Still an exception, not a ``None``, so it cannot
    be accidentally threaded into a mutation.
    """


class BoardReadIncomplete(BoardLookupError):
    """The read did not cover every project item, so absence cannot be concluded.

    Deliberately distinct from :class:`BoardItemNotFound`: conflating "I did not
    see it" with "it is not there" is the defect this module exists to stop.
    """


_GH_CACHE = None


def _default_gh():
    """The hardened ``gh`` wrapper, loaded lazily so tests never need a live ``gh``.

    Cached: the hook resolves one item per linked issue in a loop, and
    re-``exec_module``-ing gh_client on every call is pure waste.

    Loaded by explicit path rather than ``sys.path.insert``: this module's own
    directory contains names that collide with other suites' modules (e.g.
    ``recheck_parent_status``), and prepending it would shadow them for the rest
    of the process. An importing test suite must not be able to break an
    unrelated one.
    """
    global _GH_CACHE
    if _GH_CACHE is not None:
        return _GH_CACHE

    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parent / "gh_client.py"
    spec = importlib.util.spec_from_file_location("gh_client", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _GH_CACHE = module.gh
    return _GH_CACHE


def resolve_board_item_id(
    issue_number,
    *,
    repo=DEFAULT_REPO,
    kind="issue",
    project_number=PROJECT_NUMBER,
    _gh=None,
):
    """Return the board-item id for ``issue_number`` in ``repo``, or raise.

    Never returns ``None`` or an empty string. See :func:`resolve_board_item` for
    the full contract; this is the id-only convenience wrapper.
    """
    return resolve_board_item(
        issue_number, repo=repo, kind=kind, project_number=project_number, _gh=_gh
    )[0]


def resolve_board_item(
    issue_number,
    *,
    repo=DEFAULT_REPO,
    kind="issue",
    project_number=PROJECT_NUMBER,
    _gh=None,
):
    """Return ``(board_item_id, status_name)`` for ``issue_number`` in ``repo``, or raise.

    ``status_name`` is ``None`` when the card has no Status set; the *id* is
    never ``None``. Raises :class:`BoardItemNotFound` only when a complete read
    proves there is no card on ``project_number``; :class:`BoardReadIncomplete`
    when the read was partial; :class:`BoardLookupError` for a missing
    issue/PR, an absent ``id``, a failed ``gh`` call, or an unexpected shape.
    """
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
    if _gh is None:
        _gh = _default_gh()

    owner, _, name = repo.partition("/")
    if not owner or not name:
        raise BoardLookupError(f"repo must be 'owner/name', got {repo!r}")

    try:
        resp = _gh(
            "api", "graphql",
            "-f", f"query={QUERY_TEMPLATE % kind}",
            "-F", f"owner={owner}",
            "-F", f"name={name}",
            "-F", f"number={issue_number}",
        )
    except subprocess.CalledProcessError as exc:
        # gh_client.GhError subclasses CalledProcessError. A GraphQL NOT_FOUND
        # (wrong repo, or a PR number - repository.issue is Issue-only) exits
        # non-zero, so without this the caller gets a CalledProcessError that is
        # NOT catchable as BoardLookupError, breaking the module's one promise:
        # every outcome is either an id or a BoardLookupError.
        #
        # KNOWN COST, accepted (PR #1300 review, finding 1): a GraphQL NOT_FOUND
        # is an `errors` array on HTTP *200*, so gh_client's transient-detector
        # (which matches an `HTTP 4xx` token in stderr) does not recognise it and
        # RETRIES ~4x / ~14s of backoff before we get here. Deliberately not
        # "fixed" by widening that detector: it is shared by every board script,
        # and teaching it a GraphQL-specific string is the kind of denylist widen
        # that does not converge. The error path is off the hot path - the hook's
        # linked issues come from same-repo closingIssuesReferences and resolve -
        # so the cost is bounded and paid only when something is already wrong.
        stderr = (getattr(exc, "stderr", "") or "")
        hint = ""
        if "Could not resolve to an Issue" in stderr:
            hint = (
                " - #{n} does not resolve to an Issue in {repo} (a PR shares the "
                "number space but is not an Issue; this resolver handles Issues)"
            ).format(n=issue_number, repo=repo)
        raise BoardLookupError(
            f"board lookup failed for #{issue_number} in {repo}{hint}"
        ) from exc

    try:
        node = (resp["data"]["repository"] or {})[kind]
    except (KeyError, TypeError) as exc:
        raise BoardLookupError(
            f"unexpected response shape resolving #{issue_number} in {repo}: {exc}"
        ) from exc

    if node is None:
        # Not the same as "no card": the issue itself could not be resolved (wrong
        # repo, or a PR number - repository.issue is Issue-only).
        raise BoardLookupError(
            f"#{issue_number} does not resolve to an Issue in {repo} "
            "(wrong repo, or a PR number?)"
        )

    try:
        items = node["projectItems"]
        nodes = items["nodes"]
        # pageInfo is tolerated-absent on purpose. A HIT needs no completeness
        # information, and treating a missing optional as fatal made the migrated
        # hook inert against a stub that modelled the pre-migration query. Absence
        # is handled fail-SAFE below: it makes a negative conclusion unprovable,
        # never a confident "no card".
        page_info = items.get("pageInfo") or {}
    except (KeyError, TypeError) as exc:
        raise BoardLookupError(
            f"unexpected projectItems shape for #{issue_number} in {repo}: {exc}"
        ) from exc

    # A hit is trustworthy even from a partial list - the id matched is real.
    # Only the negative conclusion below needs the read to be complete.
    for item in nodes:
        if ((item.get("project") or {}).get("number")) != project_number:
            continue
        item_id = item.get("id")
        if not item_id:
            raise BoardLookupError(
                f"board item for #{issue_number} has no 'id' - the schema "
                "changed; refusing to return an unusable value"
            )
        status = None
        for fv in (item.get("fieldValues") or {}).get("nodes") or []:
            if fv and (fv.get("field") or {}).get("name") == "Status":
                status = fv.get("name")
                break
        return item_id, status

    if page_info.get("hasNextPage"):
        raise BoardReadIncomplete(
            f"#{issue_number} in {repo} has more project items than were read; "
            f"cannot conclude it has no card on project {project_number}"
        )
    if "hasNextPage" not in page_info:
        raise BoardReadIncomplete(
            f"#{issue_number} in {repo}: response carried no pageInfo, so a "
            f"complete read of its project items cannot be proven; refusing to "
            f"report 'no card on project {project_number}'"
        )

    raise BoardItemNotFound(
        f"#{issue_number} in {repo} has no card on project {project_number} "
        f"(complete read of {len(nodes)} project item(s), archived included)"
    )
