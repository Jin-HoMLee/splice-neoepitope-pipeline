#!/usr/bin/env python3
"""run_star_alignment.py — Generate splice junction quantification from raw
FASTQ files using STAR aligner.

This script aligns RNA-Seq FASTQ files using STAR and produces splice junction
quantification files in the format expected by the downstream filtering steps.

STAR outputs a file named `SJ.out.tab` containing splice junctions:
  - Column 1: chromosome
  - Column 2: intron start (1-based)
  - Column 3: intron end (1-based, inclusive)
  - Column 4: strand (0=undefined, 1=+, 2=-)
  - Column 5: intron motif (0=non-canonical, 1=GT/AG, 2=CT/AC, etc.)
  - Column 6: annotated (0=unannotated, 1=annotated)
  - Column 7: number of uniquely mapping reads crossing junction
  - Column 8: number of multi-mapping reads crossing junction
  - Column 9: maximum spliced alignment overhang

We convert this to the format expected by filter_junctions.py:
  junction_id     mapped_reads
  chr:start:end:strand    <count>

Usage (standalone):
  python run_star_alignment.py \\
      --fastq1 sample_R1.fastq.gz \\
      --fastq2 sample_R2.fastq.gz \\
      --genome-dir resources/star_index \\
      --output-dir results/raw_data/local/files \\
      --sample-id sample1 \\
      --threads 8

Usage (Snakemake):
  Called automatically by the `star_align` rule.
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def build_star_index(
    genome_fasta: str | Path,
    gtf_file: str | Path,
    index_dir: str | Path,
    threads: int = 8,
    genomeSAindexNbases: int = 14,
) -> None:
    """Build STAR genome index if it doesn't exist.

    Args:
        genome_fasta: Path to genome FASTA file.
        gtf_file: Path to GTF annotation file.
        index_dir: Output directory for STAR index.
        threads: Number of threads.
        genomeSAindexNbases: Genome SA index size (smaller for smaller genomes).
    """
    index_dir = Path(index_dir)
    if (index_dir / "SA").exists():
        log.info("STAR index already exists at %s", index_dir)
        return

    log.info("Building STAR index in %s...", index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)

    # Handle gzipped GTF
    gtf_path = Path(gtf_file)
    if str(gtf_path).endswith('.gz'):
        log.info("GTF file is gzipped, STAR will read it directly")
        # STAR can read gzipped files directly with --readFilesCommand zcat
        gtf_cmd = str(gtf_path)
    else:
        gtf_cmd = str(gtf_path)

    cmd = [
        "STAR",
        "--runMode", "genomeGenerate",
        "--runThreadN", str(threads),
        "--genomeDir", str(index_dir),
        "--genomeFastaFiles", str(genome_fasta),
        "--sjdbGTFfile", gtf_cmd,
        "--sjdbOverhang", "100",  # read length - 1, 100 is a good default
        "--genomeSAindexNbases", str(genomeSAindexNbases),
    ]

    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("STAR index generation failed:\n%s", result.stderr)
        raise RuntimeError(f"STAR indexing failed with code {result.returncode}")
    log.info("STAR index built successfully")


def run_star_alignment(
    fastq1: str | Path,
    fastq2: str | Path | None,
    genome_dir: str | Path,
    output_dir: str | Path,
    sample_id: str,
    threads: int = 8,
) -> Path:
    """Run STAR alignment and return path to junction quantification file.

    Args:
        fastq1: Path to read 1 FASTQ (can be gzipped).
        fastq2: Path to read 2 FASTQ for paired-end (None for single-end).
        genome_dir: Path to STAR genome index directory.
        output_dir: Output directory for alignment results.
        sample_id: Sample identifier for output naming.
        threads: Number of threads.

    Returns:
        Path to the converted junction quantification TSV file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = output_dir / f"{sample_id}_"

    cmd = [
        "STAR",
        "--runMode", "alignReads",
        "--runThreadN", str(threads),
        "--genomeDir", str(genome_dir),
        "--readFilesIn", str(fastq1),
        "--outFileNamePrefix", str(prefix),
        "--outSAMtype", "None",  # Don't output BAM, we only need junctions
        "--outSJfilterReads", "Unique",  # Only count uniquely mapped reads
        "--outSJfilterCountUniqueMin", "1", "1", "1", "1",
        "--outSJfilterCountTotalMin", "1", "1", "1", "1",
    ]

    # Add read 2 if paired-end
    if fastq2:
        cmd[cmd.index(str(fastq1)) + 1:cmd.index(str(fastq1)) + 1] = [str(fastq2)]
        # Fix the command - readFilesIn takes both files
        cmd_idx = cmd.index("--readFilesIn")
        cmd[cmd_idx + 1] = str(fastq1)
        cmd.insert(cmd_idx + 2, str(fastq2))

    # Handle gzipped input
    if str(fastq1).endswith('.gz'):
        cmd.extend(["--readFilesCommand", "zcat"])

    log.info("Running STAR alignment for sample %s...", sample_id)
    log.debug("Command: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("STAR alignment failed:\n%s", result.stderr)
        raise RuntimeError(f"STAR alignment failed with code {result.returncode}")

    # Convert STAR SJ.out.tab to the format expected by filter_junctions.py
    star_sj_file = output_dir / f"{sample_id}_SJ.out.tab"
    output_tsv = output_dir / f"{sample_id}.tsv"

    _convert_star_junctions(star_sj_file, output_tsv)
    log.info("Junction quantification saved to %s", output_tsv)

    return output_tsv


def _convert_star_junctions(star_sj_file: Path, output_tsv: Path) -> None:
    """Convert STAR SJ.out.tab format to pipeline junction format.

    STAR SJ.out.tab columns:
        1. chromosome
        2. intron start (1-based)
        3. intron end (1-based, inclusive)
        4. strand (0=undefined, 1=+, 2=-)
        5. intron motif
        6. annotated (0/1)
        7. unique reads
        8. multi-map reads
        9. max overhang

    Output format (tab-separated):
        junction_id     mapped_reads
        chr:start:end:strand    <unique_reads>
    """
    with star_sj_file.open() as fin, output_tsv.open("w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue

            chrom = parts[0]
            start = parts[1]  # Keep 1-based from STAR
            end = parts[2]
            strand_code = int(parts[3])
            unique_reads = int(parts[6])

            # Convert strand code
            if strand_code == 1:
                strand = "+"
            elif strand_code == 2:
                strand = "-"
            else:
                strand = "."

            # Only include junctions with reads
            if unique_reads > 0:
                junction_id = f"{chrom}:{start}:{end}:{strand}"
                fout.write(f"{junction_id}\t{unique_reads}\n")


def create_sample_manifest(
    samples: list[dict],
    manifest_path: str | Path,
) -> None:
    """Create a manifest TSV file for local samples.

    Args:
        samples: List of dicts with keys: sample_id, sample_type, fastq1, fastq2
        manifest_path: Output path for manifest TSV.
    """
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with manifest_path.open("w") as fh:
        fh.write("file_id\tfile_name\tsample_type\tproject_id\n")
        for sample in samples:
            fh.write(
                f"{sample['sample_id']}\t{sample['sample_id']}.tsv\t"
                f"{sample.get('sample_type', 'Unknown')}\tlocal\n"
            )
    log.info("Manifest written to %s (%d samples)", manifest_path, len(samples))


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    """Entry point when called from a Snakemake rule."""
    rule_name = snakemake.rule  # type: ignore[name-defined]  # noqa: F821

    if rule_name == "star_index":
        build_star_index(
            genome_fasta=snakemake.input.genome,  # type: ignore[name-defined]  # noqa: F821
            gtf_file=snakemake.input.gtf,  # type: ignore[name-defined]  # noqa: F821
            index_dir=snakemake.output.index_dir,  # type: ignore[name-defined]  # noqa: F821
            threads=snakemake.threads,  # type: ignore[name-defined]  # noqa: F821
        )
    elif rule_name == "star_align":
        run_star_alignment(
            fastq1=snakemake.input.fastq1,  # type: ignore[name-defined]  # noqa: F821
            fastq2=snakemake.input.get("fastq2"),  # type: ignore[name-defined]  # noqa: F821
            genome_dir=snakemake.input.index_dir,  # type: ignore[name-defined]  # noqa: F821
            output_dir=snakemake.params.output_dir,  # type: ignore[name-defined]  # noqa: F821
            sample_id=snakemake.wildcards.sample,  # type: ignore[name-defined]  # noqa: F821
            threads=snakemake.threads,  # type: ignore[name-defined]  # noqa: F821
        )
        Path(snakemake.output.done).touch()  # type: ignore[name-defined]  # noqa: F821


def _cli_main() -> None:
    """Command-line entry point for standalone use."""
    parser = argparse.ArgumentParser(
        description="Generate splice junction quantification using STAR aligner."
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    # index sub-command
    idx = sub.add_parser("index", help="Build STAR genome index")
    idx.add_argument("--genome-fasta", required=True, help="Genome FASTA file")
    idx.add_argument("--gtf", required=True, help="GTF annotation file")
    idx.add_argument("--output-dir", required=True, help="Output index directory")
    idx.add_argument("--threads", type=int, default=8, help="Number of threads")

    # align sub-command
    aln = sub.add_parser("align", help="Run STAR alignment")
    aln.add_argument("--fastq1", required=True, help="Read 1 FASTQ file")
    aln.add_argument("--fastq2", default=None, help="Read 2 FASTQ file (paired-end)")
    aln.add_argument("--genome-dir", required=True, help="STAR index directory")
    aln.add_argument("--output-dir", required=True, help="Output directory")
    aln.add_argument("--sample-id", required=True, help="Sample identifier")
    aln.add_argument("--threads", type=int, default=8, help="Number of threads")

    # manifest sub-command
    man = sub.add_parser("manifest", help="Create sample manifest from config")
    man.add_argument("--samples-tsv", required=True, help="Samples TSV file")
    man.add_argument("--output", required=True, help="Output manifest TSV")

    args = parser.parse_args()

    if args.mode == "index":
        build_star_index(
            genome_fasta=args.genome_fasta,
            gtf_file=args.gtf,
            index_dir=args.output_dir,
            threads=args.threads,
        )
    elif args.mode == "align":
        run_star_alignment(
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            genome_dir=args.genome_dir,
            output_dir=args.output_dir,
            sample_id=args.sample_id,
            threads=args.threads,
        )
    elif args.mode == "manifest":
        import csv
        with open(args.samples_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            samples = list(reader)
        create_sample_manifest(samples, args.output)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
