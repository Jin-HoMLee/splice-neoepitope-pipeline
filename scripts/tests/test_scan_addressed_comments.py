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


class TestSnippet:
    def test_short_untouched(self):
        assert s.snippet("hello world") == "hello world"

    def test_collapses_whitespace(self):
        assert s.snippet("a\n\n  b   c") == "a b c"

    def test_truncates_long(self):
        out = s.snippet("x" * 200, maxlen=50)
        assert len(out) == 50 and out.endswith("...")


class TestParseJsonl:
    def test_parses_and_skips_bad(self):
        text = '{"a": 1}\n\nnot json\n{"b": 2}\n'
        assert s._parse_jsonl(text) == [{"a": 1}, {"b": 2}]

    def test_empty(self):
        assert s._parse_jsonl("") == []
        assert s._parse_jsonl(None) == []


class TestOrchestrationNoNetwork:
    def test_scan_groups_pings(self, monkeypatch):
        since = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)
        monkeypatch.setattr(s, "recently_updated", lambda d: [
            {"number": 10, "title": "Alpha", "kind": "Issue"},
            {"number": 11, "title": "Beta", "kind": "PR"},
            {"number": 12, "title": "Gamma", "kind": "Issue"},
        ])

        def fake_comments(number):
            data = {
                10: '{"body": "**To:** PM\\n\\nreview", "created_at": "2026-07-04T10:00:00Z", "user": {"login": "dev"}}',
                11: '{"body": "**To:** Scientist\\n\\nnot pm", "created_at": "2026-07-04T10:00:00Z", "user": {"login": "sci"}}',
                12: "",  # comment-less issue
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

    def test_render_groups(self):
        since = datetime(2026, 7, 4, 15, 0, tzinfo=UTC)
        groups = [({"number": 10, "title": "Alpha", "kind": "Issue"},
                   [("dev", "2026-07-04T10:00:00Z", "review please")])]
        out = s.render("pm", groups, since)
        assert "#10 [Issue] Alpha" in out and "@dev" in out
