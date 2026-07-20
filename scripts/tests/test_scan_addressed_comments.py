"""Tests for `scripts/pm/scan_addressed_comments.py` (Issue #901).

Pure-helper unit tests + a network-free orchestration test that monkeypatches the
two `gh` I/O functions with fixtures. No live `gh` is reached.
"""

import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

SCRIPTS_PM = Path(__file__).parent.parent / "pm"
sys.path.insert(0, str(SCRIPTS_PM))

import scan_addressed_comments as s  # noqa: E402

UTC = timezone.utc


class TestDisplayNames:
    def test_known_roles(self):
        assert s.display_names("pm") == ["PM"]
        assert "Developer" in s.display_names("developer")
        assert "MM" in s.display_names("memory_manager")

    def test_unknown_role_raises(self):
        try:
            s.display_names("nobody")
            assert False, "expected ValueError"
        except ValueError:
            pass


class TestBodyAddressesRole:
    def test_simple_bold_to(self):
        body = "**From:** Developer -> **To:** PM\n\nplease review."
        assert s.body_addresses_role(body, s.display_names("pm")) is True

    def test_multi_recipient_list_matches_each(self):
        body = "**From:** PM -> **To:** Scientist, Developer\n\nheads up."
        assert s.body_addresses_role(body, s.display_names("scientist")) is True
        assert s.body_addresses_role(body, s.display_names("developer")) is True
        assert s.body_addresses_role(body, s.display_names("pm")) is False

    def test_short_alias_dev(self):
        body = "**From:** PM -> **To:** Dev\n\nlook at this."
        assert s.body_addresses_role(body, s.display_names("developer")) is True

    def test_memory_manager_two_word_phrase(self):
        body = "**From:** PM -> **To:** Memory Manager\n\ncommit please."
        assert s.body_addresses_role(body, s.display_names("memory_manager")) is True

    def test_role_mentioned_off_the_to_line_not_matched(self):
        # "PM" appears in the body but NOT in a To: field -> not a ping.
        body = "**From:** Developer -> **To:** Scientist\n\nask the PM later."
        assert s.body_addresses_role(body, s.display_names("pm")) is False

    def test_substring_inside_word_not_matched(self):
        body = "**From:** X -> **To:** PMgroup\n\nnope."
        assert s.body_addresses_role(body, s.display_names("pm")) is False

    def test_broadcast_all_not_matched(self):
        # `To: all` is a broadcast; a named addressee owes an ack, a broadcast
        # does not, so it is deliberately not surfaced per role.
        body = "**From:** PM -> **To:** all\n\nconvention is live."
        assert s.body_addresses_role(body, s.display_names("pm")) is False

    def test_no_to_line(self):
        assert s.body_addresses_role("just a normal comment about PM stuff", s.display_names("pm")) is False

    def test_none_body_safe(self):
        assert s.body_addresses_role(None, s.display_names("pm")) is False


class TestComputeSince:
    def test_marker_present_applies_overlap(self):
        marker = {"last_session_end_utc": "2026-07-04T12:00:00Z"}
        now = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        got = s.compute_since(marker, now)
        assert got == datetime(2026, 7, 3, 12, 0, tzinfo=UTC)  # minus 1-day overlap

    def test_marker_absent_uses_floor(self):
        now = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        assert s.compute_since(None, now) == now - timedelta(days=s.FLOOR_DAYS)

    def test_marker_malformed_ts_uses_floor(self):
        now = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        got = s.compute_since({"last_session_end_utc": "not-a-date"}, now)
        assert got == now - timedelta(days=s.FLOOR_DAYS)

    def test_marker_missing_key_uses_floor(self):
        now = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        assert s.compute_since({"schema": 1}, now) == now - timedelta(days=s.FLOOR_DAYS)


class TestSelectPings:
    def _c(self, body, created, login="dev"):
        return {"body": body, "created_at": created, "user": {"login": login}}

    def test_filters_by_since_and_role(self):
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        comments = [
            self._c("**To:** PM\n\nnew ping", "2026-07-04T10:00:00Z", "alice"),
            self._c("**To:** PM\n\nold ping", "2026-07-01T10:00:00Z", "bob"),      # too old
            self._c("**To:** Scientist\n\nnot for pm", "2026-07-04T11:00:00Z"),    # wrong role
        ]
        pings = s.select_pings(comments, s.display_names("pm"), since)
        assert len(pings) == 1
        assert pings[0][0] == "alice"

    def test_unparseable_timestamp_kept(self):
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        comments = [self._c("**To:** PM\n\nkeep me", None, "carol")]
        pings = s.select_pings(comments, s.display_names("pm"), since)
        assert len(pings) == 1 and pings[0][0] == "carol"

    def test_empty_comments_safe(self):
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        assert s.select_pings([], s.display_names("pm"), since) == []
        assert s.select_pings(None, s.display_names("pm"), since) == []


class TestSelectPingsDeclaredRaiser:
    """Issue #1240: the ping surfaces the declared raiser role, not author.login."""

    def _c(self, body, login="Jin-HoMLee"):
        return {"body": body, "created_at": "2026-07-04T10:00:00Z", "user": {"login": login}}

    def test_raiser_role_from_declared_from_field(self):
        # AC2: the raiser (index 1) is parsed from the body's From: line.
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        pings = s.select_pings(
            [self._c("**From:** Developer -> **To:** PM\n\nreview")],
            s.display_names("pm"), since,
        )
        assert pings[0][0] == "Jin-HoMLee"  # GitHub login unchanged (index 0)
        assert pings[0][1] == "Developer"   # declared raiser (index 1)

    def test_matched_pair_flip(self):
        # Same GitHub author, only the From: role flipped -> reported raiser flips.
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        p_sci = s.select_pings(
            [self._c("**From:** Scientist -> **To:** PM")], s.display_names("pm"), since)
        p_dev = s.select_pings(
            [self._c("**From:** Developer -> **To:** PM")], s.display_names("pm"), since)
        assert p_sci[0][1] == "Scientist"
        assert p_dev[0][1] == "Developer"

    def test_undeclared_raiser_is_none(self):
        # AC4: no From:/Created by: -> raiser is None (caller falls back to login).
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        pings = s.select_pings([self._c("**To:** PM\n\nno attribution")], s.display_names("pm"), since)
        assert pings[0][1] is None


class TestSnippet:
    def test_short_untouched(self):
        assert s.snippet("hello world") == "hello world"

    def test_collapses_whitespace(self):
        assert s.snippet("a\n\n  b   c") == "a b c"

    def test_truncates_long(self):
        out = s.snippet("x" * 200, maxlen=50)
        assert len(out) == 50 and out.endswith("...")


class TestFetchComments:
    """fetch_comments: the `--jq` -> `--paginate --slurp` + Python-filter rewrite (#1055).

    Locks in what the live CLI-contract diff could not guarantee it exercised - the
    array-of-arrays flatten and the #1011 null-`.user` guard - by feeding a fake `gh`.
    """

    @staticmethod
    def _stub(pages):
        def fake_gh(*args, **kwargs):
            if isinstance(pages, Exception):
                raise pages
            return pages
        return fake_gh

    def test_flattens_pages_and_guards_null_user(self, monkeypatch):
        # --slurp yields a list of per-page lists; the 2nd comment has a null
        # `.user` (deleted account, the #1011 case).
        pages = [
            [
                {"body": "hi", "created_at": "2026-07-01T00:00:00Z", "user": {"login": "alice"}},
                {"body": "ghost", "created_at": "2026-07-02T00:00:00Z", "user": None},
            ],
            [
                {"body": "page2", "created_at": "2026-07-03T00:00:00Z", "user": {"login": "bob"}},
            ],
        ]
        monkeypatch.setattr(s, "gh", self._stub(pages))
        out = s.fetch_comments(42)
        assert [c["user"]["login"] for c in out] == ["alice", "?", "bob"]  # flatten + guard
        assert [c["body"] for c in out] == ["hi", "ghost", "page2"]

    def test_missing_user_key_also_guards(self, monkeypatch):
        monkeypatch.setattr(s, "gh", self._stub([[{"body": "x", "created_at": "t"}]]))
        assert s.fetch_comments(1) == [{"body": "x", "created_at": "t", "user": {"login": "?"}}]

    def test_flat_page_defensive_fallback(self, monkeypatch):
        # A page that is a bare dict (not a list) is yielded, not iterated by key.
        monkeypatch.setattr(s, "gh", self._stub([{"body": "x", "created_at": "t", "user": {"login": "z"}}]))
        assert s.fetch_comments(1) == [{"body": "x", "created_at": "t", "user": {"login": "z"}}]

    def test_gh_error_returns_empty(self, monkeypatch):
        monkeypatch.setattr(s, "gh", self._stub(s.GhError(1, ["gh"], stderr="HTTP 503")))
        assert s.fetch_comments(1) == []

    def test_json_decode_error_isolated_returns_empty(self, monkeypatch):
        # A 0-exit non-JSON body raises JSONDecodeError from inside gh(); isolate it
        # per-item rather than aborting the whole scan (#1055 review).
        import json
        monkeypatch.setattr(s, "gh", self._stub(json.JSONDecodeError("bad", "doc", 0)))
        assert s.fetch_comments(1) == []


class TestOrchestrationNoNetwork:
    def test_scan_groups_pings(self, monkeypatch):
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        monkeypatch.setattr(s, "recently_updated", lambda d: [
            {"number": 10, "title": "Alpha", "kind": "Issue"},
            {"number": 11, "title": "Beta", "kind": "PR"},
            {"number": 12, "title": "Gamma", "kind": "Issue"},
        ])

        def fake_comments(number):
            # New contract (#1055): fetch_comments returns a list of comment dicts
            # (was NDJSON text parsed by the now-removed _parse_jsonl).
            data = {
                10: [{"body": "**To:** PM\n\nreview", "created_at": "2026-07-04T10:00:00Z", "user": {"login": "dev"}}],
                11: [{"body": "**To:** Scientist\n\nnot pm", "created_at": "2026-07-04T10:00:00Z", "user": {"login": "sci"}}],
                12: [],  # comment-less issue
            }
            return data[number]

        monkeypatch.setattr(s, "fetch_comments", fake_comments)
        groups = s.scan("pm", since)
        assert [g[0]["number"] for g in groups] == [10]
        assert groups[0][1][0][0] == "dev"

    def test_render_none(self):
        since = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        out = s.render("pm", [], since)
        assert "(none)" in out

    def test_render_groups_shows_declared_raiser(self):
        # AC2: render surfaces the declared raiser role, keeping the login visible.
        since = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        groups = [({"number": 10, "title": "Alpha", "kind": "Issue"},
                   [("Jin-HoMLee", "Developer", "2026-07-04T10:00:00Z", "review please")])]
        out = s.render("pm", groups, since)
        assert "#10 [Issue] Alpha" in out
        assert "Developer" in out and "@Jin-HoMLee" in out

    def test_render_undeclared_raiser_fallback(self):
        # AC4: no declared role -> show the login plus an "undeclared role" marker.
        since = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        groups = [({"number": 10, "title": "Alpha", "kind": "Issue"},
                   [("Jin-HoMLee", None, "2026-07-04T10:00:00Z", "review please")])]
        out = s.render("pm", groups, since)
        assert "@Jin-HoMLee" in out and "undeclared role" in out


class TestCliGuards:
    def _run(self, *args):
        import subprocess
        return subprocess.run(
            [sys.executable, str(SCRIPTS_PM / "scan_addressed_comments.py"), *args],
            capture_output=True, text=True, timeout=10,
        )

    def test_negative_days_rejected(self):
        r = self._run("--role", "pm", "--days", "-3")
        assert r.returncode == 2 and "non-negative" in r.stderr

    def test_unknown_role_rejected(self):
        r = self._run("--role", "nobody")
        assert r.returncode == 2

    def test_missing_role_rejected(self):
        r = self._run("--days", "1")
        assert r.returncode == 2
