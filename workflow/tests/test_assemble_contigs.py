"""Tests for assemble_contigs.py — contig assembly helpers."""

import pandas as pd

from assemble_contigs import (
    _build_downstream_bed,
    _build_upstream_bed,
    _has_soft_clip,
    _parse_fasta,
    _write_zero_stats,
    assemble_contigs,
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


class TestEmptyPipelineEarlyReturns:
    """Issue #215 follow-up: when the pipeline produces zero candidates, the
    early-return paths must still write a zero-count stats TSV so the
    cross-step aggregator does not fail with missing input.
    """

    _NOVEL_COLUMNS = [
        "junction_id", "chrom", "start", "end", "strand",
        "junction_origin", "sample_id", "sample_type",
    ]

    def _expected_zero_stats(self):
        return pd.DataFrame(
            [
                {"category": "contigs_written", "count": 0},
                {"category": "skipped_softclip", "count": 0},
                {"category": "skipped_length", "count": 0},
            ],
            columns=["category", "count"],
        )

    def test_write_zero_stats_helper_emits_three_rows(self, tmp_path):
        out = tmp_path / "stats.tsv"
        _write_zero_stats(out)
        got = pd.read_csv(out, sep="\t")
        pd.testing.assert_frame_equal(got, self._expected_zero_stats())

    def test_write_zero_stats_helper_no_path_is_noop(self, tmp_path):
        # No exception, no file created
        _write_zero_stats(None)
        assert list(tmp_path.iterdir()) == []

    def test_empty_novel_junctions_writes_zero_stats(self, tmp_path):
        novel_tsv = tmp_path / "novel.tsv"
        pd.DataFrame(columns=self._NOVEL_COLUMNS).to_csv(
            novel_tsv, sep="\t", index=False
        )
        fasta_out = tmp_path / "contigs.fa"
        stats_out = tmp_path / "contig_stats.tsv"

        assemble_contigs(
            novel_junctions_tsv=novel_tsv,
            genome_fasta=tmp_path / "fake.fa",  # never read in early-return path
            output_fasta=fasta_out,
            stats_output_path=stats_out,
        )

        assert fasta_out.exists()
        got = pd.read_csv(stats_out, sep="\t")
        pd.testing.assert_frame_equal(got, self._expected_zero_stats())

    def test_only_normal_shared_writes_zero_stats(self, tmp_path):
        novel_tsv = tmp_path / "novel.tsv"
        pd.DataFrame({
            "junction_id": ["j1", "j2"],
            "chrom": ["chr22", "chr22"],
            "start": [100, 200],
            "end": [200, 300],
            "strand": ["+", "+"],
            "junction_origin": ["normal_shared", "normal_shared"],
            "sample_id": ["s1", "s1"],
            "sample_type": ["tumor", "tumor"],
        }).to_csv(novel_tsv, sep="\t", index=False)
        fasta_out = tmp_path / "contigs.fa"
        stats_out = tmp_path / "contig_stats.tsv"

        assemble_contigs(
            novel_junctions_tsv=novel_tsv,
            genome_fasta=tmp_path / "fake.fa",
            output_fasta=fasta_out,
            stats_output_path=stats_out,
        )

        assert fasta_out.exists()
        got = pd.read_csv(stats_out, sep="\t")
        pd.testing.assert_frame_equal(got, self._expected_zero_stats())
