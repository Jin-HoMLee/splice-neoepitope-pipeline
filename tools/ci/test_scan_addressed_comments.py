"""Unit tests for scripts/pm/scan_addressed_comments.py fetch_comments (Issue #1055).

These two scan scripts have no other unit tests; the migration off `--jq` to
`--paginate --slurp` + a Python filter rests on the live CLI-contract diff, which is
only as strong as the window's data. This locks in the two things that diff might not
have exercised: the array-of-arrays flatten and the #1011 null-`.user` guard.
"""

import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import scan_addressed_comments as sac


def _stub_gh(pages):
    """Replace sac.gh with a stub returning `pages` (or raising if it's an Exception)."""
    def fake_gh(*args, **kwargs):
        if isinstance(pages, Exception):
            raise pages
        return pages
    return fake_gh


class TestFetchComments:
    def test_flattens_pages_and_guards_null_user(self, monkeypatch):
        # --slurp yields a list of per-page lists; the second comment has a null
        # `.user` (deleted account, the #1011 case).
        pages = [
            [
                {"body": "hi", "created_at": "2026-07-01T00:00:00Z",
                 "user": {"login": "alice"}},
                {"body": "ghost", "created_at": "2026-07-02T00:00:00Z",
                 "user": None},
            ],
            [
                {"body": "page2", "created_at": "2026-07-03T00:00:00Z",
                 "user": {"login": "bob"}},
            ],
        ]
        monkeypatch.setattr(sac, "gh", _stub_gh(pages))
        out = sac.fetch_comments(42)
        assert [c["user"]["login"] for c in out] == ["alice", "?", "bob"]  # flatten + guard
        assert [c["body"] for c in out] == ["hi", "ghost", "page2"]
        assert out[1]["created_at"] == "2026-07-02T00:00:00Z"

    def test_missing_user_key_also_guards(self, monkeypatch):
        pages = [[{"body": "x", "created_at": "t"}]]  # no `user` key at all
        monkeypatch.setattr(sac, "gh", _stub_gh(pages))
        assert sac.fetch_comments(1) == [
            {"body": "x", "created_at": "t", "user": {"login": "?"}}
        ]

    def test_flat_page_list_defensive_fallback(self, monkeypatch):
        # If a page ever comes back as a bare dict rather than a list, the
        # isinstance guard still yields it rather than iterating its keys.
        pages = [{"body": "x", "created_at": "t", "user": {"login": "z"}}]
        monkeypatch.setattr(sac, "gh", _stub_gh(pages))
        assert sac.fetch_comments(1) == [
            {"body": "x", "created_at": "t", "user": {"login": "z"}}
        ]

    def test_gh_error_returns_empty(self, monkeypatch):
        monkeypatch.setattr(sac, "gh", _stub_gh(sac.GhError(1, ["gh"], stderr="HTTP 503")))
        assert sac.fetch_comments(1) == []

    def test_json_decode_error_isolated_returns_empty(self, monkeypatch):
        # A 0-exit non-JSON body raises JSONDecodeError from inside gh(); fetch_comments
        # must isolate it (per-item), not abort the whole scan (Issue #1055 review).
        monkeypatch.setattr(sac, "gh", _stub_gh(json.JSONDecodeError("bad", "doc", 0)))
        assert sac.fetch_comments(1) == []
