"""Tests for `.agents/hooks/_shell_parse.py` (Issue #1130).

The two `gh`-matching PostToolUse hooks never fired on a heredoc-created PR -
which is how essentially every PR with a real body is opened - because heredoc
bodies tokenized into the command stream and a newline was not a separator. Two
shipped automations (board-add/Status, and the auto-requested bot review) were
silently dead on the dominant path.

The tension these tests hold: the fix must make the heredoc form match WITHOUT
reintroducing the PR #558 false positive, where a `gh pr create` quoted *inside*
a comment body must not read as an invocation.
"""

import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
sys.path.insert(0, str(HOOKS_DIR))

import _shell_parse as sp  # noqa: E402


class TestStripHeredocBodies:
    def test_body_and_delimiter_removed_command_line_kept(self):
        cmd = "cat > b.md <<'EOF'\nBody text here.\nEOF\ngh pr create --fill"
        out = sp.strip_heredoc_bodies(cmd)
        assert "Body text here." not in out
        assert "cat > b.md" in out          # the opening command survives
        assert "gh pr create --fill" in out  # the command AFTER the heredoc survives

    def test_body_prose_cannot_be_mistaken_for_a_command(self):
        # The whole point of removing (not tokenizing) the body.
        cmd = "cat > b.md <<'EOF'\ngh pr create --fill\nEOF\necho done"
        out = sp.strip_heredoc_bodies(cmd)
        assert "gh pr create" not in out

    def test_unquoted_and_double_quoted_and_dash_forms(self):
        for opener in ("<<EOF", "<<'EOF'", '<<"EOF"', "<<-EOF"):
            cmd = f"cat > b.md {opener}\nsecret body\nEOF\ngh pr create --fill"
            out = sp.strip_heredoc_bodies(cmd)
            assert "secret body" not in out, opener
            assert "gh pr create" in out, opener

    def test_custom_delimiter(self):
        cmd = "cat > b.md <<'XEOF'\nbody\nXEOF\ngh pr create --fill"
        assert "body" not in sp.strip_heredoc_bodies(cmd)

    def test_no_heredoc_is_a_passthrough(self):
        assert sp.strip_heredoc_bodies("gh pr create --fill") == "gh pr create --fill"

    def test_unterminated_heredoc_consumes_the_rest(self):
        # No closing delimiter: the rest is body. Reading it as body is correct -
        # and fail-safe, since it can only suppress a match, never invent one.
        cmd = "cat > b.md <<'EOF'\nbody line\nmore body"
        out = sp.strip_heredoc_bodies(cmd)
        assert "body line" not in out and "more body" not in out


class TestNewlinesToSeparators:
    def test_unquoted_newline_becomes_a_separator(self):
        assert sp.newlines_to_separators("echo hi\ngh pr create") == "echo hi;gh pr create"

    def test_newline_inside_single_quotes_is_preserved(self):
        assert sp.newlines_to_separators("echo 'a\nb'") == "echo 'a\nb'"

    def test_newline_inside_double_quotes_is_preserved(self):
        assert sp.newlines_to_separators('echo "a\nb"') == 'echo "a\nb"'

    def test_no_newline_is_a_passthrough(self):
        assert sp.newlines_to_separators("gh pr create --fill") == "gh pr create --fill"


class TestNormalizeCommand:
    def test_the_real_pr_creation_shape(self):
        """The exact shape that has been silently missing the hook."""
        cmd = (
            'S=/tmp/x\n'
            'cat > "$S/pr.md" <<\'EOF\'\n'
            'Some body (with parens) and prose.\n'
            'EOF\n'
            'gh pr create --title "t" --body-file "$S/pr.md" --base main'
        )
        out = sp.normalize_command(cmd)
        assert "Some body" not in out       # heredoc body gone
        assert ";gh pr create" in out       # and gh is now at a command start

    def test_none_is_safe(self):
        assert sp.normalize_command(None) == ""
