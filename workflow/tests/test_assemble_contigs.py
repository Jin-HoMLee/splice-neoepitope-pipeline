"""Tests for assemble_contigs.py — contig assembly helpers."""

import pandas as pd

from assemble_contigs import (
    _build_downstream_bed,
    _build_upstream_bed,
    _has_soft_clip,
    _parse_fasta,
)


class TestHasSoftClip:
    def test_all_uppercase_returns_false(self):
        assert _has_soft_clip("ATCGATCGATCG") is False

    def test_lowercase_bases_return_true(self):
        assert _has_soft_clip("ATCGatcg") is True

    def test_all_lowercase_returns_true(self):
        assert _has_soft_clip("atcgatcg") is True

    def test_empty_string_returns_false(self):
        assert _has_soft_clip("") is False


class TestBuildUpstreamBed:
    def _make_junctions(self):
        return pd.DataFrame({
            "chrom": ["chr22"],
            "start": [200],   # 0-based junction start
            "end":   [300],
            "strand": ["+"],
            "junction_id": ["chr22:201:300:+"],
        })

    def test_upstream_end_equals_junction_start(self):
        df = self._make_junctions()
        bed = _build_upstream_bed(df, upstream_nt=26)
        assert bed.iloc[0]["end"] == 200

    def test_upstream_start_is_junction_start_minus_nt(self):
        df = self._make_junctions()
        bed = _build_upstream_bed(df, upstream_nt=26)
        assert bed.iloc[0]["start"] == 174

    def test_upstream_start_clipped_at_zero(self):
        df = pd.DataFrame({
            "chrom": ["chr22"], "start": [10], "end": [200],
            "strand": ["+"], "junction_id": ["chr22:11:200:+"],
        })
        bed = _build_upstream_bed(df, upstream_nt=26)
        assert bed.iloc[0]["start"] == 0


class TestBuildDownstreamBed:
    def _make_junctions(self):
        return pd.DataFrame({
            "chrom": ["chr22"],
            "start": [200],
            "end":   [300],   # 0-based junction end
            "strand": ["+"],
            "junction_id": ["chr22:201:300:+"],
        })

    def test_downstream_start_equals_junction_end(self):
        df = self._make_junctions()
        bed = _build_downstream_bed(df, downstream_nt=24)
        assert bed.iloc[0]["start"] == 300

    def test_downstream_end_is_junction_end_plus_nt(self):
        df = self._make_junctions()
        bed = _build_downstream_bed(df, downstream_nt=24)
        assert bed.iloc[0]["end"] == 324


class TestParseFasta:
    def test_parses_single_record(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">header1\nATCGATCG\n")
        result = _parse_fasta(fa)
        assert result == {"header1": "ATCGATCG"}

    def test_parses_multiple_records(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">seq1\nAAAA\n>seq2\nCCCC\n")
        result = _parse_fasta(fa)
        assert result["seq1"] == "AAAA"
        assert result["seq2"] == "CCCC"

    def test_strips_bedtools_coordinate_suffix(self, tmp_path):
        # bedtools getfasta -name appends ::chrom:start-end to the header
        fa = tmp_path / "test.fa"
        fa.write_text(">chr22:201:300:+::chr22:174-200\nATCG\n")
        result = _parse_fasta(fa)
        assert "chr22:201:300:+" in result

    def test_empty_file_returns_empty_dict(self, tmp_path):
        fa = tmp_path / "empty.fa"
        fa.write_text("")
        assert _parse_fasta(fa) == {}
