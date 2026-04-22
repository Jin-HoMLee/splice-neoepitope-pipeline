"""Tests for proteome_filter.py — filtering peptides against the human proteome via k-mer index."""

import csv
from pathlib import Path

from proteome_filter import _build_kmer_index, _parse_accession, proteome_filter


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_fasta(tmp_path: Path, records: list[tuple[str, str]]) -> Path:
    """Write a minimal Swiss-Prot-style FASTA to tmp_path/proteome.fasta."""
    fasta = tmp_path / "proteome.fasta"
    with fasta.open("w") as fh:
        for accession, seq in records:
            fh.write(f">sp|{accession}|PROT_HUMAN Description\n{seq}\n")
    return fasta


def _make_peptides_tsv(tmp_path: Path, rows: list[dict]) -> Path:
    tsv = tmp_path / "peptides.tsv"
    with tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["contig_key", "start_nt", "peptide"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    return tsv


def _read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# _parse_accession
# ---------------------------------------------------------------------------

class TestParseAccession:
    def test_swissprot_header(self):
        assert _parse_accession("sp|P12345|PROT_HUMAN Description") == "P12345"

    def test_trembl_header(self):
        assert _parse_accession("tr|A0A000|UNKNOWN_HUMAN Protein") == "A0A000"

    def test_nonstandard_header_falls_back_to_first_token(self):
        assert _parse_accession("CUSTOM_ACC some description") == "CUSTOM_ACC"


# ---------------------------------------------------------------------------
# _build_kmer_index
# ---------------------------------------------------------------------------

class TestBuildKmerIndex:
    def test_contains_kmers_from_protein(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "ACDEFGHIKLM")])
        index = _build_kmer_index(fasta, [9])
        assert "ACDEFGHIK" in index
        assert "CDEFGHIKL" in index
        assert "DEFGHIKLM" in index

    def test_multi_length(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "ACDEFGHIKLM")])
        index = _build_kmer_index(fasta, [8, 9, 10])
        assert "ACDEFGHI" in index       # 8-mer
        assert "ACDEFGHIK" in index      # 9-mer
        assert "ACDEFGHIKL" in index     # 10-mer

    def test_absent_sequence_not_in_index(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "ACDEFGHIKLM")])
        index = _build_kmer_index(fasta, [9])
        assert "NOVELPEPBD" not in index

    def test_first_accession_recorded(self, tmp_path):
        fasta = _make_fasta(tmp_path, [
            ("P11111", "ACDEFGHIKLM"),
            ("P22222", "ACDEFGHIKLM"),
        ])
        index = _build_kmer_index(fasta, [9])
        assert index["ACDEFGHIK"][0] == "P11111"

    def test_all_accessions_recorded(self, tmp_path):
        fasta = _make_fasta(tmp_path, [
            ("P11111", "ACDEFGHIKLM"),
            ("P22222", "ACDEFGHIKLM"),
        ])
        index = _build_kmer_index(fasta, [9])
        assert set(index["ACDEFGHIK"]) == {"P11111", "P22222"}

    def test_no_duplicate_accessions_per_kmer(self, tmp_path):
        # Same accession appears twice via different proteins with shared k-mer
        fasta = _make_fasta(tmp_path, [
            ("P11111", "ACDEFGHIKLM"),
            ("P11111", "XACDEFGHIKLY"),  # same accession, different context
        ])
        index = _build_kmer_index(fasta, [9])
        assert index["ACDEFGHIK"].count("P11111") == 1

    def test_empty_fasta_returns_empty_index(self, tmp_path):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        assert _build_kmer_index(fasta, [9]) == {}

    def test_protein_shorter_than_kmer_length_skipped(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "ACDE")])
        index = _build_kmer_index(fasta, [9])
        assert len(index) == 0


# ---------------------------------------------------------------------------
# proteome_filter
# ---------------------------------------------------------------------------

class TestProteomeFilter:
    def test_splits_novel_and_excluded(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "SELFPEPAACD")])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "NOVELPEPB"},
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[9],
        )
        novel = _read_tsv(tmp_path / "novel.tsv")
        excluded = _read_tsv(tmp_path / "excluded.tsv")
        assert len(novel) == 1
        assert novel[0]["peptide"] == "NOVELPEPB"
        assert len(excluded) == 1
        assert excluded[0]["peptide"] == "SELFPEPAA"

    def test_matched_accessions_column(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "SELFPEPAACD")])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[9],
        )
        excluded = _read_tsv(tmp_path / "excluded.tsv")
        assert excluded[0]["matched_accessions"] == "P11111"

    def test_multiple_accessions_semicolon_separated(self, tmp_path):
        fasta = _make_fasta(tmp_path, [
            ("P11111", "SELFPEPAACD"),
            ("P22222", "SELFPEPAAEF"),
        ])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[9],
        )
        excluded = _read_tsv(tmp_path / "excluded.tsv")
        accessions = excluded[0]["matched_accessions"].split(";")
        assert set(accessions) == {"P11111", "P22222"}

    def test_all_novel_when_no_matches(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "AAAAAAAAA")])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "NOVELPEPB"},
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[9],
        )
        assert len(_read_tsv(tmp_path / "novel.tsv")) == 1
        assert len(_read_tsv(tmp_path / "excluded.tsv")) == 0

    def test_all_excluded(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "SELFPEPAACD")])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[9],
        )
        assert len(_read_tsv(tmp_path / "novel.tsv")) == 0
        assert len(_read_tsv(tmp_path / "excluded.tsv")) == 1

    def test_same_peptide_from_two_contigs_both_excluded(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "SELFPEPAACD")])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "SELFPEPAA"},
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[9],
        )
        excluded = _read_tsv(tmp_path / "excluded.tsv")
        assert len(excluded) == 2
        assert excluded[0]["contig_key"] == "c1"
        assert excluded[1]["contig_key"] == "c2"

    def test_output_headers(self, tmp_path):
        fasta = _make_fasta(tmp_path, [("P11111", "SELFPEPAACD")])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "NOVELPEPB"},
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[9],
        )
        novel_header = (tmp_path / "novel.tsv").read_text().splitlines()[0]
        excl_header = (tmp_path / "excluded.tsv").read_text().splitlines()[0]
        assert novel_header == "contig_key\tstart_nt\tpeptide"
        assert excl_header == "contig_key\tstart_nt\tpeptide\tmatched_accessions"

    def test_multi_length_filter(self, tmp_path):
        # proteome contains both an 8-mer and a 10-mer; a 9-mer stays novel
        fasta = _make_fasta(tmp_path, [("P11111", "ACDEFGHICD" "KLMNPQRSTU")])
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "ACDEFGHI"},    # 8-mer, in proteome
            {"contig_key": "c2", "start_nt": "0", "peptide": "NOVELPEPB"},   # 9-mer, not in proteome
            {"contig_key": "c3", "start_nt": "0", "peptide": "KLMNPQRSTU"},  # 10-mer, in proteome
        ])
        proteome_filter(
            peptides_tsv=peptides_tsv,
            proteome_fasta=fasta,
            novel_tsv=tmp_path / "novel.tsv",
            excluded_tsv=tmp_path / "excluded.tsv",
            peptide_lengths=[8, 9, 10],
        )
        novel = _read_tsv(tmp_path / "novel.tsv")
        excluded = _read_tsv(tmp_path / "excluded.tsv")
        assert len(novel) == 1
        assert novel[0]["peptide"] == "NOVELPEPB"
        assert len(excluded) == 2
