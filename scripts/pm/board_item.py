#!/usr/bin/env python3
"""Board-item lookup that refuses to hand back a result it cannot vouch for (Issue #1151)."""

OWNER = "Jin-HoMLee"
PROJECT_NUMBER = 9
DEFAULT_REPO = "Jin-HoMLee/splice-neoepitope-pipeline"


class BoardLookupError(Exception):
    """Base: the lookup cannot vouch for its answer, so it refuses to return one."""


class BoardItemNotFound(BoardLookupError):
    """A COMPLETE board read did not contain the issue. The only trustworthy 'no'."""


class BoardReadIncomplete(BoardLookupError):
    """The board read did not deliver every item, so absence cannot be concluded.

    Deliberately distinct from :class:`BoardItemNotFound`: conflating "I did not
    see it" with "it is not there" is the exact defect this module exists to stop.
    """

QUERY = """
query($owner: String!, $number: Int!, $after: String) {
  user(login: $owner) {
    projectV2(number: $number) {
      items(first: 100, after: $after) {
        totalCount
        pageInfo { hasNextPage endCursor }
        nodes {
          id
          content {
            ... on Issue { number repository { nameWithOwner } }
            ... on PullRequest { number repository { nameWithOwner } }
          }
        }
      }
    }
  }
}
"""


def _default_gh():
    """The hardened ``gh`` wrapper, loaded lazily so tests never need a live ``gh``.

    Loaded by explicit path rather than ``sys.path.insert``: this module's own
    directory contains names that collide with other suites' modules (e.g.
    ``recheck_parent_status``), and prepending it would shadow them for the rest
    of the process. An importing test suite must not be able to break an
    unrelated one.
    """
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parent / "gh_client.py"
    spec = importlib.util.spec_from_file_location("gh_client", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.gh


def resolve_board_item_id(
    issue_number,
    *,
    repo=DEFAULT_REPO,
    owner=OWNER,
    project_number=PROJECT_NUMBER,
    _gh=None,
):
    """Resolve ``issue_number`` in ``repo`` to its board-item id, or raise.

    ``repo`` is load-bearing, not decoration: board 9 spans two repos with
    independent numbering and the collision is live (#183 and #193 exist in
    both, both carded). Matching on number alone returns a page-order-dependent
    card - a silent wrong answer of exactly the kind this module exists to stop.
    """
    if _gh is None:
        _gh = _default_gh()
    cursor = None
    seen = 0
    total = None
    while True:
        args = [
            "api", "graphql",
            "-f", f"query={QUERY}",
            "-F", f"owner={owner}",
            "-F", f"number={project_number}",
        ]
        if cursor is not None:
            args += ["-F", f"after={cursor}"]
        resp = _gh(*args)

        try:
            items = resp["data"]["user"]["projectV2"]["items"]
            total = items["totalCount"]
            nodes = items["nodes"]
            page_info = items["pageInfo"]
        except (KeyError, TypeError) as exc:
            raise BoardLookupError(
                f"unexpected board response shape for {owner}/{project_number} "
                f"(wrong owner or project number?): {exc}"
            ) from exc
        seen += len(nodes)

        # A hit is trustworthy from a partial read: the id we matched is real.
        # Only the NEGATIVE conclusion below needs the read to be complete.
        for node in nodes:
            content = node.get("content") or {}
            node_repo = (content.get("repository") or {}).get("nameWithOwner")
            if content.get("number") == issue_number and node_repo == repo:
                item_id = node.get("id")
                if not item_id:
                    raise BoardLookupError(
                        f"board item for #{issue_number} has no 'id' - the schema "
                        "changed; refusing to return an unusable value"
                    )
                return item_id

        if not page_info["hasNextPage"]:
            break
        cursor = page_info["endCursor"]

    if seen != total:
        raise BoardReadIncomplete(
            f"read {seen} of {total} items on project {owner}/{project_number}; "
            f"cannot conclude whether #{issue_number} is absent"
        )

    raise BoardItemNotFound(
        f"issue #{issue_number} is not on project {owner}/{project_number} "
        f"(complete read of {total} items)"
    )
