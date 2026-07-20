"""Tests for `scripts/pm/declared_role.py` (Issue #1240).

Every persona posts on GitHub under Jin-Ho's one account, so `author.login` is
always `Jin-HoMLee` and carries zero role signal. The source of truth for who
raised a coordination item is the declared `**From:** <Role>` / `**Created by:**
<Role>` line in the body. `parse_declared_role` is the pure parser both the inbox
scan and the dispatch digest use so neither ever eyeballs `author.login`.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pm"))
import declared_role as dr  # noqa: E402


class TestMatchedPairFalsifier:
    """AC3: same GitHub author, one body variable flipped, reported role flips."""

    def test_from_scientist_reports_scientist(self):
        body = "**From:** Scientist -> **To:** PM\n\nplease review."
        assert dr.parse_declared_role(body) == "Scientist"

    def test_flip_body_role_flips_result(self):
        # Identical shape, only the From: role changed -> the reported raiser
        # must change with it. If it didn't, the parser is reading something
        # other than the declared role (the falsifier).
        body = "**From:** Developer -> **To:** PM\n\nplease review."
        assert dr.parse_declared_role(body) == "Developer"


class TestFromField:
    def test_to_role_is_not_leaked_as_raiser(self):
        # The To: recipient (PM) must never be reported as the raiser.
        body = "**From:** Scientist -> **To:** PM"
        assert dr.parse_declared_role(body) == "Scientist"

    def test_multi_recipient_to_does_not_confuse_from(self):
        body = "**From:** PM -> **To:** Scientist, Developer\n\nheads up."
        assert dr.parse_declared_role(body) == "PM"

    def test_short_alias_dev(self):
        assert dr.parse_declared_role("**From:** Dev -> **To:** PM") == "Developer"

    def test_two_word_memory_manager(self):
        assert dr.parse_declared_role("**From:** Memory Manager -> **To:** PM") == "Memory Manager"

    def test_mm_short_alias(self):
        assert dr.parse_declared_role("**From:** MM -> **To:** Developer") == "Memory Manager"


class TestCreatedByFallback:
    def test_created_by_when_no_from(self):
        body = "Some issue body.\n\n**Created by:** Scientist"
        assert dr.parse_declared_role(body) == "Scientist"

    def test_from_takes_precedence_over_created_by(self):
        # A coordination comment can carry both; From: is the raiser.
        body = "**From:** Developer -> **To:** PM\n\n**Created by:** PM"
        assert dr.parse_declared_role(body) == "Developer"


class TestUndeclared:
    def test_no_declared_role_returns_none(self):
        assert dr.parse_declared_role("just a plain comment, no attribution") is None

    def test_empty_body_returns_none(self):
        assert dr.parse_declared_role("") is None
        assert dr.parse_declared_role(None) is None

    def test_unrecognized_role_returns_none(self):
        # A From: line naming something outside the role vocabulary is undeclared,
        # not a false match.
        assert dr.parse_declared_role("**From:** Marketing -> **To:** PM") is None
