"""Tests for crosscheck_annotate_flag.py — CI canary cross-checking the
home-rolled ``annotated`` flag against ``regtools junctions annotate``.

[Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377):
A second independent source of truth (regtools' own annotator) makes the
class of semantic-misread bug behind [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370)
(anchor outers vs. intron donor/acceptor) loud instead of silent.

These cover the *pure* comparison logic — no regtools binary required, so they
run in the lightweight ``pipeline-pytest`` job. The end-to-end run against the
real ``regtools`` binary lives in the ``annotate-flag-canary`` CI job.

Coordinate facts (verified empirically, regtools 1.0.0):
- ``regtools junctions annotate`` emits ``start`` = intron donor (0-based) and
  ``end`` = acceptor-exclusive **+ 1**, so the home-rolled key (``chrom``,
  ``start``, ``end_exclusive``, ``strand``) is recovered as
  ``(chrom, start, end - 1, strand)``.
- ``known_junction`` is column 14 (1-indexed): ``1`` = annotated, ``0`` = novel.
"""

import pytest

from crosscheck_annotate_flag import (
    crosscheck,
    homerolled_flags,
    parse_regtools_annotate,
)

# A faithful regtools annotate output header + rows (subset of columns is fine
# for parsing; the parser keys off position, not trailing columns).
_HEADER = (
    "chrom\tstart\tend\tname\tscore\tstrand\tsplice_site\tacceptors_skipped"
    "\texons_skipped\tdonors_skipped\tanchor\tknown_donor\tknown_acceptor"
    "\tknown_junction\tgene_names\tgene_ids\ttranscripts"
)


def _row(chrom, start, end, strand, known):
    return (
        f"{chrom}\t{start}\t{end}\tJX\t10\t{strand}\tGT-AG\t0\t0\t0\tDA"
        f"\t1\t1\t{known}\tGENE\tENSG\tTX"
    )


class TestParseRegtoolsAnnotate:
    def test_maps_end_minus_one_and_known_flag(self):
        # regtools end=10527853 -> home-rolled exclusive end 10527852
        lines = [_HEADER, _row("chr22", 10524446, 10527853, "-", 1)]
        flags = parse_regtools_annotate(lines)
        assert flags == {("chr22", 10524446, 10527852, "-"): 1}

    def test_novel_junction_flag_zero(self):
        lines = [_HEADER, _row("chr22", 11000000, 11001001, "+", 0)]
        flags = parse_regtools_annotate(lines)
        assert flags == {("chr22", 11000000, 11001000, "+"): 0}

    def test_skips_header_blank_and_comment_lines(self):
        lines = [
            _HEADER,
            "",
            "# a comment",
            _row("chr22", 100, 201, "+", 1),
        ]
        flags = parse_regtools_annotate(lines)
        assert flags == {("chr22", 100, 200, "+"): 1}

    def test_empty_input_returns_empty(self):
        assert parse_regtools_annotate([_HEADER]) == {}


class TestHomerolledFlags:
    def test_membership_drives_flag(self):
        ref = frozenset({("chr22", 10524446, 10527852, "-")})
        ids = [
            "chr22:10524447:10527852:-",  # in ref  -> 1
            "chr22:11000001:11001000:+",  # not in ref -> 0
        ]
        flags = homerolled_flags(ids, ref)
        assert flags == {
            ("chr22", 10524446, 10527852, "-"): 1,
            ("chr22", 11000000, 11001000, "+"): 0,
        }

    def test_unparseable_id_skipped(self):
        flags = homerolled_flags(["garbage", "chr22:100:200:+"], frozenset())
        assert flags == {("chr22", 99, 200, "+"): 0}


class TestCrosscheck:
    def test_perfect_agreement(self):
        rt = {("chr22", 100, 200, "+"): 1, ("chr22", 300, 400, "-"): 0}
        hr = dict(rt)
        res = crosscheck(rt, hr)
        assert res.n_compared == 2
        assert res.agreement == 1.0
        assert res.coverage == 1.0
        assert res.disagreements == []

    def test_single_disagreement_listed(self):
        rt = {("chr22", 100, 200, "+"): 1, ("chr22", 300, 400, "-"): 1}
        hr = {("chr22", 100, 200, "+"): 1, ("chr22", 300, 400, "-"): 0}
        res = crosscheck(rt, hr)
        assert res.n_compared == 2
        assert res.agreement == 0.5
        assert (("chr22", 300, 400, "-"), 1, 0) in res.disagreements

    def test_partial_coverage_key_only_one_side(self):
        rt = {("chr22", 100, 200, "+"): 1, ("chr22", 300, 400, "-"): 1}
        hr = {("chr22", 100, 200, "+"): 1}  # missing the 2nd key
        res = crosscheck(rt, hr)
        # intersection is 1 of 2 keys -> coverage 0.5, agreement 1.0 over it
        assert res.n_compared == 1
        assert res.agreement == 1.0
        assert res.coverage == 0.5

    def test_disjoint_keys_no_zero_division(self):
        rt = {("chr22", 100, 200, "+"): 1}
        hr = {("chr22", 999, 1000, "+"): 1}
        res = crosscheck(rt, hr)
        assert res.n_compared == 0
        # vacuous intersection must NOT read as healthy agreement
        assert res.agreement == 0.0
        assert res.coverage == 0.0

    def test_both_empty_no_zero_division(self):
        res = crosscheck({}, {})
        assert res.n_compared == 0
        assert res.agreement == 0.0
        assert res.coverage == 0.0
