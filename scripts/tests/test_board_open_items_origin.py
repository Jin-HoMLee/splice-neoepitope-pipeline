# scripts/tests/test_board_open_items_origin.py
#
# Issue #999 - board #9 aggregates the project repo AND the personas repo, whose
# numbers collide. The text table must disambiguate a same-numbered project vs
# personas item so a reader can't misattribute one (the misread that produced
# this Issue's own original, incorrect filing).
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import board_open_items as b  # noqa: E402

PROJECT_URL = "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/71"
PERSONAS_URL = "https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/issues/71"


class TestOriginFromUrl:
    def test_project(self):
        assert b.origin_from_url(PROJECT_URL) == "project"

    def test_personas_not_shadowed_by_substring(self):
        # PERSONAS_REPO contains PROJECT_REPO as a substring; personas must win.
        assert b.origin_from_url(PERSONAS_URL) == "personas"

    def test_pr_url(self):
        assert b.origin_from_url(
            "https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/pull/12"
        ) == "personas"

    def test_unknown_repo_is_other(self):
        assert b.origin_from_url("https://github.com/someone/other-repo/issues/3") == "other"

    def test_empty_or_none_is_other(self):
        assert b.origin_from_url("") == "other"
        assert b.origin_from_url(None) == "other"


class TestRefCell:
    def test_project_is_bare(self):
        assert b.ref_cell({"number": 71, "origin": "project"}) == "71"

    def test_personas_is_tagged(self):
        assert b.ref_cell({"number": 71, "origin": "personas"}) == "pers#71"

    def test_other_is_tagged(self):
        assert b.ref_cell({"number": 71, "origin": "other"}) == "ext#71"

    def test_missing_origin_falls_back_to_bare(self):
        # Defensive: an item dict without an origin key renders the bare number.
        assert b.ref_cell({"number": 71}) == "71"


def test_same_number_project_and_personas_are_distinguishable():
    """The core #999 guarantee: two items sharing #71 render distinct refs."""
    proj = {"number": 71, "origin": b.origin_from_url(PROJECT_URL)}
    pers = {"number": 71, "origin": b.origin_from_url(PERSONAS_URL)}
    assert b.ref_cell(proj) != b.ref_cell(pers)
    assert b.ref_cell(proj) == "71" and b.ref_cell(pers) == "pers#71"


def test_format_table_tags_personas_row_only():
    """End-to-end: the rendered table shows a bare project ref and a tagged
    personas ref for two same-numbered items."""
    items = [
        {"number": 71, "url": PROJECT_URL, "origin": "project", "title": "project item",
         "kind": "Issue", "is_draft": False, "is_parent": False, "status": "Ready",
         "priority": "P2", "size": "S", "role": "role:developer", "arc": None,
         "arc_phase": None, "updated_at": None},
        {"number": 71, "url": PERSONAS_URL, "origin": "personas", "title": "personas item",
         "kind": "Issue", "is_draft": False, "is_parent": False, "status": "Ready",
         "priority": "P2", "size": "S", "role": "role:memory_manager", "arc": None,
         "arc_phase": None, "updated_at": None},
    ]
    out = b.format_table(items)
    assert "pers#71" in out          # personas row tagged
    assert " 71 " in out or out.count("71") >= 2  # project row keeps a bare 71
    # header renamed # -> Ref
    assert "Ref" in out.splitlines()[0]
