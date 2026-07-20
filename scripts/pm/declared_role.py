#!/usr/bin/env python3
"""Parse the declared raiser role from a GitHub body (Issue #1240).

Every persona posts under Jin-Ho's single GitHub account, so `author.login` is
always `Jin-HoMLee` and carries no role signal. The authoritative source of who
raised a coordination item is the declared line in the body:

    **From:** <Role> -> **To:** <Role(s)>      (coordination comments/Discussions)
    **Created by:** <Role>                      (issue/PR/comment attribution)

`parse_declared_role` is the single pure parser the inbox scan
(`scan_addressed_comments.py`) and the dispatch digest (`dispatch_digest.py`) both
use, so neither ever reads `author.login` as the raiser. `ROLE_ALIASES` is defined
here as the one canonical role vocabulary; `scan_addressed_comments` imports it.
"""

# Canonical role slug -> the display aliases that appear in a `**From:**` /
# `**Created by:**` / `**To:**` field. Full names plus the short board forms
# (Sci / Dev / MM). This is the single source of the role vocabulary.
ROLE_ALIASES = {
    "pm": ["PM"],
    "scientist": ["Scientist", "Sci"],
    "developer": ["Developer", "Dev"],
    "memory_manager": ["Memory Manager", "MM"],
}

# Canonical slug -> display name reported to the reader. Derived from the single
# source (the first, canonical alias) so it cannot drift from ROLE_ALIASES.
ROLE_DISPLAY = {slug: aliases[0] for slug, aliases in ROLE_ALIASES.items()}

# (alias_lower, slug) pairs, longest alias first so a multi-word alias
# ("memory manager") is matched before any shorter one could partial-hit.
_ALIAS_TO_SLUG = sorted(
    ((alias.lower(), slug) for slug, aliases in ROLE_ALIASES.items() for alias in aliases),
    key=lambda pair: len(pair[0]),
    reverse=True,
)

_FROM_MARKER = "From:**"
_CREATED_BY_MARKER = "Created by:**"
# Separators after which the `**From:**` field's value ends (the To: recipient
# must never leak in as the raiser).
_FROM_TERMINATORS = ("->", "→", "To:**", "**To")
# A `**Created by:**` footer names one role; cut at a trailing parenthetical so
# `**Created by:** PM (reviewed by Developer)` reports PM, not the longest match.
_CREATED_BY_TERMINATORS = ("(",)


def _flatten(text):
    """Lowercase, punctuation-to-space, space-padded - for whole-phrase matching."""
    words = "".join(c if c.isalnum() else " " for c in (text or "").lower()).split()
    return f" {' '.join(words)} "


def _match_role(text):
    """Canonical slug of the first role alias appearing as a whole phrase, or None."""
    flat = _flatten(text)
    for alias_lower, slug in _ALIAS_TO_SLUG:
        if f" {alias_lower} " in flat:
            return slug
    return None


def _field_value(body, marker, terminators=()):
    """Text after `marker` on the first line carrying it, truncated at any terminator.

    Returns None when no line carries the marker.
    """
    for line in (body or "").splitlines():
        idx = line.find(marker)
        if idx == -1:
            continue
        tail = line[idx + len(marker):]
        cut = len(tail)
        for sep in terminators:
            j = tail.find(sep)
            if j != -1:
                cut = min(cut, j)
        return tail[:cut]
    return None


def parse_declared_role(body):
    """Display name of the declared raiser, or None when undeclared/unrecognized.

    Precedence: a present `**From:**` is authoritative - it is the raiser slot, so
    it short-circuits (an unrecognized From: role returns None rather than falling
    through to the author-attribution field, which would misreport). Only when no
    From: field is present does the `**Created by:**` attribution stand in. Returns
    None when neither field is present or the resolved value is outside the role
    vocabulary, so callers can fall back to `author.login` with an "undeclared
    role" marker.
    """
    for marker, terminators in (
        (_FROM_MARKER, _FROM_TERMINATORS),
        (_CREATED_BY_MARKER, _CREATED_BY_TERMINATORS),
    ):
        value = _field_value(body, marker, terminators)
        if value is not None:
            slug = _match_role(value)
            return ROLE_DISPLAY[slug] if slug is not None else None
    return None
