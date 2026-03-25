#!/usr/bin/env python3
"""run_netmhcpan.py — Wrapper to run NetMHCPan 4.1 and parse the output.

NetMHCPan 4.1 must be installed separately under an academic licence from DTU
Bioinformatics.  See README.md for installation instructions.

The script:
  1. Runs ``netMHCpan`` with the 16-mer peptide FASTA as input, requesting all
     possible 9-mer sub-peptides (``-l 9``) and IC50 output (``-BA``).
  2. Parses the raw NetMHCPan output into a structured TSV.
  3. Classifies each 9-mer as a strong binder (IC50 < 50 nM), weak binder
     (IC50 < 500 nM), or non-binder.

Output TSV columns:
  peptide_16mer  position  peptide_9mer  allele  ic50_nM  rank  binder_class
  source_header

Usage (standalone):
  python run_netmhcpan.py \\
      --peptides-fasta results/peptides/TCGA-BRCA/peptides.fa \\
      --output results/predictions/TCGA-BRCA/predictions.tsv \\
      --executable netMHCpan \\
      --allele HLA-A02:01 \\
      --peptide-length 9 \\
      --ic50-strong 50 \\
      --ic50-weak 500

Usage (Snakemake):
  Called automatically by the ``run_netmhcpan`` rule.
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# NetMHCPan runner
# ---------------------------------------------------------------------------

def _run_netmhcpan(
    executable: str,
    peptides_fasta: str | Path,
    allele: str,
    peptide_length: int,
    output_raw: str | Path,
) -> None:
    """Invoke NetMHCPan 4.1 and save the raw output.

    Args:
        executable:     Path or name of the netMHCpan binary.
        peptides_fasta: FASTA of 16-mer input peptides.
        allele:         HLA allele in NetMHCPan format (e.g. ``HLA-A02:01``).
        peptide_length: Sliding-window peptide length (9 for 9-mers).
        output_raw:     File to capture the raw NetMHCPan stdout.

    Raises:
        FileNotFoundError:          if the executable cannot be found.
        subprocess.CalledProcessError: if NetMHCPan returns non-zero.
    """
    cmd = [
        executable,
        "-f", str(peptides_fasta),
        "-a", allele,
        "-l", str(peptide_length),
        "-BA",         # include binding affinity (IC50 nM)
        "-inptype", "0",  # input is FASTA
    ]
    log.info("Running: %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        log.error(
            "NetMHCPan executable not found: %r. "
            "Please install NetMHCPan 4.1 and ensure it is on your PATH "
            "or set netmhcpan.executable in config/config.yaml.",
            executable,
        )
        raise
    except subprocess.CalledProcessError as exc:
        log.error("NetMHCPan failed (exit %d):\n%s", exc.returncode, exc.stderr)
        raise

    Path(output_raw).write_text(result.stdout)
    if result.stderr:
        log.warning("NetMHCPan stderr:\n%s", result.stderr)


# ---------------------------------------------------------------------------
# Output parser
# ---------------------------------------------------------------------------

# Pattern matching a data line in NetMHCPan 4.1 output:
# Pos  HLA  Peptide  ... nM  Rank  ... BindLevel
_DATA_LINE_RE = re.compile(
    r"^\s*(\d+)\s+(HLA-\S+)\s+(\w+)"  # pos  allele  peptide
    r".*?"
    r"(\d+\.\d+)\s+"                  # IC50 nM
    r"(\d+\.\d+)"                     # %Rank
    r".*?(?:(SB|WB))?\s*$",           # optional binder label
    re.MULTILINE,
)

# Simpler fallback: try to extract columns by position (NetMHCPan output is
# fixed-width).  We try both parsing strategies.
_HEADER_SENTINEL = "Pos"


def _parse_netmhcpan_output(
    raw_output: str,
    ic50_strong: float,
    ic50_weak: float,
) -> list[dict]:
    """Parse NetMHCPan 4.1 tabular output into a list of dicts.

    NetMHCPan prints a fixed-width table.  We parse it line by line, looking
    for lines that start with a position integer.

    Args:
        raw_output:  String containing the full NetMHCPan stdout.
        ic50_strong: IC50 threshold (nM) for strong binders.
        ic50_weak:   IC50 threshold (nM) for weak binders.

    Returns:
        List of dicts with keys:
        peptide_16mer, position, peptide_9mer, allele, ic50_nM, rank,
        binder_class, source_header.
    """
    records: list[dict] = []
    current_source = ""

    for line in raw_output.splitlines():
        # Track which peptide is currently being predicted
        if line.startswith("# Sequence number"):
            # Format: # Sequence number X : <header>
            parts = line.split(":", 1)
            if len(parts) == 2:
                current_source = parts[1].strip()
            continue

        # Try regex match on data lines
        line_stripped = line.strip()
        if not line_stripped or line_stripped.startswith("#"):
            continue

        # Split by whitespace — NetMHCPan 4.1 output columns:
        # Pos HLA Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL
        #     Score_BA Affinity(nM) %Rank_BA BindLevel
        parts = line_stripped.split()
        if len(parts) < 14:
            continue
        try:
            pos = int(parts[0])
        except ValueError:
            continue

        try:
            allele   = parts[1]
            peptide  = parts[2]
            # Score_BA is at index 13, Affinity(nM) at 14, %Rank_BA at 15
            # (column positions may vary slightly between versions; we look
            # for the nM value as a large float)
            ic50_nM  = float(parts[14]) if len(parts) > 14 else float("inf")
            rank     = float(parts[15]) if len(parts) > 15 else float("inf")
            label    = parts[16] if len(parts) > 16 else ""
        except (ValueError, IndexError):
            continue

        if ic50_nM <= ic50_strong:
            binder_class = "strong"
        elif ic50_nM <= ic50_weak:
            binder_class = "weak"
        else:
            binder_class = "non"

        records.append(
            {
                "source_header": current_source,
                "peptide_16mer": current_source.split("|")[0] if current_source else "",
                "position": pos,
                "peptide_9mer": peptide,
                "allele": allele,
                "ic50_nM": ic50_nM,
                "rank": rank,
                "binder_class": binder_class,
            }
        )

    return records


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_prediction(
    peptides_fasta: str | Path,
    output_tsv: str | Path,
    executable: str = "netMHCpan",
    allele: str = "HLA-A02:01",
    peptide_length: int = 9,
    ic50_strong: float = 50.0,
    ic50_weak: float = 500.0,
) -> None:
    """Run NetMHCPan and write parsed results to a TSV file.

    Args:
        peptides_fasta: FASTA of 16-mer peptides.
        output_tsv:     Destination TSV file.
        executable:     NetMHCPan 4.1 binary name or path.
        allele:         HLA allele (NetMHCPan format).
        peptide_length: Sliding-window peptide length.
        ic50_strong:    Strong-binder IC50 threshold (nM).
        ic50_weak:      Weak-binder IC50 threshold (nM).
    """
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as tmp:
        raw_output_path = tmp.name

    try:
        _run_netmhcpan(executable, peptides_fasta, allele, peptide_length, raw_output_path)
        raw_output = Path(raw_output_path).read_text()
    finally:
        try:
            os.unlink(raw_output_path)
        except OSError:
            pass

    records = _parse_netmhcpan_output(raw_output, ic50_strong, ic50_weak)
    df = pd.DataFrame(
        records,
        columns=[
            "source_header", "peptide_16mer", "position", "peptide_9mer",
            "allele", "ic50_nM", "rank", "binder_class",
        ],
    )
    df.to_csv(output_tsv, sep="\t", index=False)
    log.info(
        "Predictions: %d total, %d strong, %d weak binders → %s",
        len(df),
        (df["binder_class"] == "strong").sum(),
        (df["binder_class"] == "weak").sum(),
        output_tsv,
    )


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    run_prediction(
        peptides_fasta=snakemake.input.peptides_fasta,  # type: ignore[name-defined]  # noqa: F821
        output_tsv=snakemake.output.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        executable=snakemake.params.executable,  # type: ignore[name-defined]  # noqa: F821
        allele=snakemake.params.hla_allele,  # type: ignore[name-defined]  # noqa: F821
        peptide_length=snakemake.params.peptide_length,  # type: ignore[name-defined]  # noqa: F821
        ic50_strong=float(snakemake.params.ic50_strong),  # type: ignore[name-defined]  # noqa: F821
        ic50_weak=float(snakemake.params.ic50_weak),  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Run NetMHCPan 4.1 epitope prediction on 16-mer peptides."
    )
    parser.add_argument("--peptides-fasta", required=True, help="Input peptides FASTA")
    parser.add_argument("--output", required=True, help="Output predictions TSV")
    parser.add_argument("--executable", default="netMHCpan", help="NetMHCPan binary")
    parser.add_argument("--allele", default="HLA-A02:01", help="HLA allele")
    parser.add_argument("--peptide-length", type=int, default=9)
    parser.add_argument("--ic50-strong", type=float, default=50.0)
    parser.add_argument("--ic50-weak", type=float, default=500.0)
    args = parser.parse_args()

    run_prediction(
        peptides_fasta=args.peptides_fasta,
        output_tsv=args.output,
        executable=args.executable,
        allele=args.allele,
        peptide_length=args.peptide_length,
        ic50_strong=args.ic50_strong,
        ic50_weak=args.ic50_weak,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
