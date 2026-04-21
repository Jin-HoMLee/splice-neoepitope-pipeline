#!/usr/bin/env python3
"""blastp_filter.py — Filter translated 9-mers against the human reference proteome.

Runs blastp against the UniProt Swiss-Prot human proteome and excludes any
peptide with a perfect full-length match (100% identity, full query coverage).
Peptides that already exist in the canonical human proteome are self-peptides
the immune system is tolerized to and cannot be neoantigens.

Output:
  peptides_novel.tsv    — same schema as peptides.tsv, exact matches removed
  peptides_excluded.tsv — excluded peptides with matching UniProt accession
  blastp_hits.tsv       — raw blastp output (outfmt 6), kept for audit

Usage (Snakemake):
  Called automatically by the blastp_filter_peptides rule.

Usage (standalone):
  python blastp_filter.py \\
      --peptides-tsv results/.../peptides/peptides.tsv \\
      --db-prefix resources/blastp_db/human_proteome \\
      --novel-tsv results/.../peptides/peptides_novel.tsv \\
      --excluded-tsv results/.../peptides/peptides_excluded.tsv \\
      --blastp-hits results/.../peptides/blastp_hits.tsv
"""

import argparse
import csv
import logging
import subprocess
import tempfile
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def _write_peptide_fasta(peptides: list[dict], fasta_path: Path) -> None:
    """Write unique peptide sequences to FASTA for blastp input."""
    seen = set()
    with fasta_path.open("w") as fh:
        for row in peptides:
            pep = row["peptide"]
            if pep not in seen:
                fh.write(f">{pep}\n{pep}\n")
                seen.add(pep)
    log.info("Wrote %d unique peptides to %s", len(seen), fasta_path)


def _run_blastp(
    fasta_path: Path,
    db_prefix: str,
    hits_path: Path,
    threads: int,
) -> None:
    """Run blastp with settings optimised for short peptides."""
    cmd = [
        "blastp",
        "-query", str(fasta_path),
        "-db", db_prefix,
        "-task", "blastp-short",
        "-word_size", "2",
        "-evalue", "200000",
        "-qcov_hsp_perc", "100",
        "-outfmt", "6 qseqid sseqid pident",
        "-num_threads", str(threads),
        "-out", str(hits_path),
    ]
    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"blastp failed:\n{result.stderr}")
    log.info("blastp completed")


def _parse_exact_matches(hits_path: Path) -> dict[str, str]:
    """Return {peptide: first_matching_accession} for full-length exact hits.

    Full query coverage is guaranteed by -qcov_hsp_perc 100 at search time.
    Only 100% identity hits are retained (int conversion avoids float equality
    issues with BLAST's decimal output, e.g. "100.000").
    """
    matches: dict[str, str] = {}
    if not hits_path.exists() or hits_path.stat().st_size == 0:
        return matches
    with hits_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            qseqid, sseqid, pident = parts[:3]
            if int(float(pident)) == 100:
                if qseqid not in matches:
                    matches[qseqid] = sseqid
    return matches


def blastp_filter(
    peptides_tsv: Path,
    db_prefix: str,
    novel_tsv: Path,
    excluded_tsv: Path,
    hits_path: Path,
    threads: int = 1,
) -> None:
    peptides_tsv = Path(peptides_tsv)
    novel_tsv = Path(novel_tsv)
    excluded_tsv = Path(excluded_tsv)
    hits_path = Path(hits_path)

    for p in (novel_tsv, excluded_tsv, hits_path):
        p.parent.mkdir(parents=True, exist_ok=True)

    with peptides_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        peptides = list(reader)

    log.info("Loaded %d peptides from %s", len(peptides), peptides_tsv)

    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False) as tmp:
        tmp_fasta = Path(tmp.name)

    try:
        _write_peptide_fasta(peptides, tmp_fasta)
        _run_blastp(tmp_fasta, db_prefix, hits_path, threads)
    finally:
        tmp_fasta.unlink(missing_ok=True)

    exact_matches = _parse_exact_matches(hits_path)
    log.info("Found %d peptide sequences with exact proteome matches", len(exact_matches))

    fieldnames = reader.fieldnames  # preserve whatever columns the input has
    excluded_fieldnames = list(reader.fieldnames) + ["matched_accession"]

    with novel_tsv.open("w", newline="") as novel_fh, \
         excluded_tsv.open("w", newline="") as excl_fh:
        novel_writer = csv.DictWriter(novel_fh, fieldnames=fieldnames, delimiter="\t")
        excl_writer = csv.DictWriter(excl_fh, fieldnames=excluded_fieldnames, delimiter="\t")
        novel_writer.writeheader()
        excl_writer.writeheader()

        n_novel = n_excluded = 0
        for row in peptides:
            pep = row["peptide"]
            if pep in exact_matches:
                excl_writer.writerow({**row, "matched_accession": exact_matches[pep]})
                n_excluded += 1
            else:
                novel_writer.writerow(row)
                n_novel += 1

    log.info(
        "Result: %d novel peptides passed, %d excluded (%.1f%% reduction)",
        n_novel, n_excluded,
        100 * n_excluded / len(peptides) if peptides else 0,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry points
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    sm = snakemake  # type: ignore[name-defined]  # noqa: F821
    log_file = sm.log[0]
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    blastp_filter(
        peptides_tsv=sm.input.peptides_tsv,
        db_prefix=sm.params.db_prefix,
        novel_tsv=sm.output.novel_tsv,
        excluded_tsv=sm.output.excluded_tsv,
        hits_path=sm.output.blastp_hits,
        threads=sm.threads,
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter translated peptides against the human proteome via blastp."
    )
    parser.add_argument("--peptides-tsv", required=True)
    parser.add_argument("--db-prefix", required=True)
    parser.add_argument("--novel-tsv", required=True)
    parser.add_argument("--excluded-tsv", required=True)
    parser.add_argument("--blastp-hits", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    blastp_filter(
        peptides_tsv=args.peptides_tsv,
        db_prefix=args.db_prefix,
        novel_tsv=args.novel_tsv,
        excluded_tsv=args.excluded_tsv,
        hits_path=args.blastp_hits,
        threads=args.threads,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
