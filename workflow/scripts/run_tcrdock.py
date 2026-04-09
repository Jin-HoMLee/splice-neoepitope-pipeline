#!/usr/bin/env python3
"""run_tcrdock.py — Run TCRdock structural prediction on top neoepitope candidates.

TCRdock (Bradley et al.) predicts the 3D structure of a TCR-peptide-MHC
ternary complex using a modified AlphaFold v2 multimer backend adapted
specifically for TCR:pMHC complexes.

Requirements
------------
  - Linux x64 with NVIDIA GPU
  - TCRdock installed (MIT license, no registration required)
    https://github.com/phbradley/TCRdock
  - AlphaFold v2 parameters (~3.5 GB, CC BY 4.0)

This script is NOT compatible with macOS arm64 (M1/M2). It is intended to
run on a GCP Spot GPU VM (n1-standard-4 + NVIDIA T4). The Snakemake rule
that calls this script is only included when config[tcrdock][enabled] is true,
so local / CPU-only runs are unaffected.

Candidate selection
-------------------
Currently selects the top N strong binders by IC50 (MHCflurry). Once
TRUST4 + ProTCR (#24) is implemented, selection should switch to the top
ProTCR-ranked candidates instead.

Inputs
------
  predictions_tsv: MHCflurry predictions TSV
                   (columns: contig_key, start_nt, peptide, allele,
                    ic50_nM, percentile_rank, binder_class)

Outputs
-------
  top_candidate.pdb    — predicted TCR-peptide-MHC ternary complex (PDB format)
  docking_scores.tsv   — docking geometry metrics per candidate
                         (columns: peptide, allele, va_gene, vb_gene,
                          cdr3a, cdr3b, rmsd, tcr_pdb_rmsd, docking_geometry)

Usage (standalone, Linux + GPU only):
  python run_tcrdock.py \\
      --predictions-tsv results/predictions/local/predictions.tsv \\
      --output-pdb results/predictions/local/tcrdock/top_candidate.pdb \\
      --output-scores results/predictions/local/tcrdock/docking_scores.tsv \\
      --tcrdock-dir /opt/tcrdock \\
      --alphafold-params /opt/af_params

Usage (Snakemake):
  Called automatically by the run_tcrdock rule when tcrdock.enabled is true.
"""

import argparse
import csv
import json
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Candidate selection
# ---------------------------------------------------------------------------

def select_top_candidates(
    predictions_tsv: str | Path,
    n_candidates: int = 1,
) -> pd.DataFrame:
    """Select the top N strong binders by IC50 from MHCflurry predictions.

    TODO (#24): Once ProTCR scores are available, replace IC50 ranking with
    ProTCR score ranking to select the most immunogenic candidates rather than
    just the best MHC binders.

    Args:
        predictions_tsv: MHCflurry predictions TSV.
        n_candidates:    Number of candidates to return.

    Returns:
        DataFrame of top candidates, sorted by ic50_nM ascending.
    """
    df = pd.read_csv(predictions_tsv, sep="\t")
    strong = df[df["binder_class"] == "strong"].sort_values("ic50_nM")
    if strong.empty:
        log.warning("No strong binders found — trying weak binders as fallback.")
        strong = df[df["binder_class"] == "weak"].sort_values("ic50_nM")
    if strong.empty:
        log.error("No strong or weak binders found. Cannot run TCRdock.")
        return pd.DataFrame()
    return strong.head(n_candidates).reset_index(drop=True)


# ---------------------------------------------------------------------------
# TCRdock input preparation
# ---------------------------------------------------------------------------

def build_tcrdock_input(
    candidates: pd.DataFrame,
    fallback_hla: dict,
    fallback_tcr: dict,
    output_dir: Path,
) -> Path:
    """Write TCRdock input TSV for the selected candidates.

    TCRdock expects a TSV with columns:
        organism  mhc_class  mhc  peptide  va_gene  ja_gene  cdr3a
        vb_gene  jb_gene  cdr3b

    Fallback HLA uses the A-family allele (HLA-A*02:01 by default) since
    MHCflurry predictions are currently run for a single allele. When HLA
    typing (#23) is implemented, the per-sample allele should be used instead.

    Args:
        candidates:   DataFrame of top candidates from select_top_candidates().
        fallback_hla: Dict with keys A, B, C (allele strings).
        fallback_tcr: Dict with va_gene, ja_gene, cdr3a, vb_gene, jb_gene, cdr3b.
        output_dir:   Directory to write the input TSV.

    Returns:
        Path to the written input TSV.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    input_tsv = output_dir / "tcrdock_input.tsv"

    rows = []
    for _, row in candidates.iterrows():
        # Use the allele from the prediction; fall back to config A-family allele
        # Convert MHCflurry format (HLA-A*02:01) to TCRdock format (HLA-A*02:01)
        allele = row.get("allele", fallback_hla["A"])
        rows.append({
            "organism":  "human",
            "mhc_class": "1",
            "mhc":       allele,
            "peptide":   row["peptide"],
            "va_gene":   fallback_tcr["va_gene"],
            "ja_gene":   fallback_tcr["ja_gene"],
            "cdr3a":     fallback_tcr["cdr3a"],
            "vb_gene":   fallback_tcr["vb_gene"],
            "jb_gene":   fallback_tcr["jb_gene"],
            "cdr3b":     fallback_tcr["cdr3b"],
        })

    with input_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    log.info("TCRdock input written to %s (%d candidates)", input_tsv, len(rows))
    return input_tsv


# ---------------------------------------------------------------------------
# TCRdock execution
# ---------------------------------------------------------------------------

def run_tcrdock(
    input_tsv: Path,
    output_dir: Path,
    tcrdock_dir: str | Path,
    alphafold_params_dir: str | Path,
) -> Path:
    """Run TCRdock and return the directory containing output PDB files.

    Args:
        input_tsv:            TCRdock input TSV (from build_tcrdock_input).
        output_dir:           Directory to write TCRdock outputs.
        tcrdock_dir:          Path to TCRdock installation.
        alphafold_params_dir: Path to AlphaFold v2 params/ directory.

    Returns:
        Path to TCRdock output directory.
    """
    tcrdock_dir = Path(tcrdock_dir).expanduser()
    alphafold_params_dir = Path(alphafold_params_dir).expanduser()
    predict_script = tcrdock_dir / "predict.py"

    if not predict_script.exists():
        raise FileNotFoundError(
            f"TCRdock predict.py not found at {predict_script}. "
            f"Check config[tcrdock][tcrdock_dir]."
        )

    cmd = [
        sys.executable, str(predict_script),
        "--targets_tsvfile", str(input_tsv),
        "--outfile_prefix", str(output_dir / "tcrdock_out"),
        "--params_dir", str(alphafold_params_dir),
    ]

    log.info("Running TCRdock: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("TCRdock failed:\n%s", result.stderr)
        raise RuntimeError(f"TCRdock exited with code {result.returncode}")

    log.info("TCRdock completed successfully.")
    return output_dir


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

def collect_outputs(
    tcrdock_output_dir: Path,
    candidates: pd.DataFrame,
    output_pdb: Path,
    output_scores: Path,
    fallback_tcr: dict,
) -> None:
    """Collect TCRdock PDB and docking scores into pipeline output files.

    TCRdock writes one PDB per candidate. This function:
      1. Copies the top candidate PDB to output_pdb.
      2. Parses docking geometry from the TCRdock score file and writes
         a clean TSV to output_scores.

    Args:
        tcrdock_output_dir: Directory containing TCRdock output files.
        candidates:         DataFrame of candidates (same order as input TSV).
        output_pdb:         Destination PDB path (top candidate only).
        output_scores:      Destination docking scores TSV.
        fallback_tcr:       TCR config dict (for annotating the scores TSV).
    """
    # TCRdock names PDB files by row index: tcrdock_out_0.pdb, tcrdock_out_1.pdb ...
    top_pdb = tcrdock_output_dir / "tcrdock_out_0.pdb"
    if not top_pdb.exists():
        # Try alternative naming conventions
        pdbs = sorted(tcrdock_output_dir.glob("*.pdb"))
        if not pdbs:
            raise FileNotFoundError(
                f"No PDB files found in {tcrdock_output_dir}. TCRdock may have failed."
            )
        top_pdb = pdbs[0]

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_bytes(top_pdb.read_bytes())
    log.info("Top candidate PDB written to %s", output_pdb)

    # Parse docking geometry from TCRdock JSON/TSV score output
    score_files = list(tcrdock_output_dir.glob("*score*.json")) + \
                  list(tcrdock_output_dir.glob("*geometry*.tsv"))

    score_rows = []
    for i, (_, cand) in enumerate(candidates.iterrows()):
        row = {
            "peptide":          cand["peptide"],
            "allele":           cand.get("allele", ""),
            "ic50_nM":          cand["ic50_nM"],
            "va_gene":          fallback_tcr["va_gene"],
            "ja_gene":          fallback_tcr["ja_gene"],
            "cdr3a":            fallback_tcr["cdr3a"],
            "vb_gene":          fallback_tcr["vb_gene"],
            "jb_gene":          fallback_tcr["jb_gene"],
            "cdr3b":            fallback_tcr["cdr3b"],
            "pdb_file":         f"tcrdock_out_{i}.pdb",
            "docking_geometry": "",
        }
        # If a JSON score file exists for this candidate, extract geometry
        json_file = tcrdock_output_dir / f"tcrdock_out_{i}_scores.json"
        if json_file.exists():
            try:
                scores = json.loads(json_file.read_text())
                row["docking_geometry"] = json.dumps(scores)
            except Exception:
                pass
        score_rows.append(row)

    pd.DataFrame(score_rows).to_csv(output_scores, sep="\t", index=False)
    log.info("Docking scores written to %s", output_scores)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_structural_validation(
    predictions_tsv: str | Path,
    output_pdb: str | Path,
    output_scores: str | Path,
    tcrdock_dir: str | Path,
    alphafold_params_dir: str | Path,
    n_candidates: int = 1,
    fallback_hla: dict | None = None,
    fallback_tcr: dict | None = None,
) -> None:
    """Run TCRdock structural validation on top neoepitope candidates.

    Args:
        predictions_tsv:      MHCflurry predictions TSV.
        output_pdb:           Output PDB file (top candidate).
        output_scores:        Output docking scores TSV.
        tcrdock_dir:          Path to TCRdock installation.
        alphafold_params_dir: Path to AlphaFold v2 params/.
        n_candidates:         Number of top candidates to model.
        fallback_hla:         Fallback HLA alleles dict (keys: A, B, C).
        fallback_tcr:         Fallback TCR sequences dict.
    """
    if fallback_hla is None:
        fallback_hla = {"A": "HLA-A*02:01", "B": "HLA-B*07:02", "C": "HLA-C*07:02"}
    if fallback_tcr is None:
        fallback_tcr = {
            "va_gene": "TRAV12-2", "ja_gene": "TRAJ21", "cdr3a": "CAVNFGGGKLI",
            "vb_gene": "TRBV6-5",  "jb_gene": "TRBJ2-7", "cdr3b": "CASSLAGGRPEQYF",
        }

    output_pdb = Path(output_pdb)
    output_scores = Path(output_scores)
    work_dir = output_pdb.parent / "tcrdock_workdir"

    # Step 1: Select top candidates
    candidates = select_top_candidates(predictions_tsv, n_candidates)
    if candidates.empty:
        log.error("No candidates available. Skipping TCRdock.")
        return

    log.info(
        "Top %d candidate(s) selected for TCRdock:\n%s",
        len(candidates),
        candidates[["peptide", "allele", "ic50_nM"]].to_string(index=False),
    )

    # Step 2: Build TCRdock input TSV
    input_tsv = build_tcrdock_input(candidates, fallback_hla, fallback_tcr, work_dir)

    # Step 3: Run TCRdock
    tcrdock_output_dir = run_tcrdock(input_tsv, work_dir, tcrdock_dir, alphafold_params_dir)

    # Step 4: Collect outputs
    collect_outputs(tcrdock_output_dir, candidates, output_pdb, output_scores, fallback_tcr)


# ---------------------------------------------------------------------------
# Snakemake / CLI entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))

    run_structural_validation(
        predictions_tsv=snakemake.input.predictions_tsv,  # type: ignore[name-defined]  # noqa: F821
        output_pdb=snakemake.output.pdb,  # type: ignore[name-defined]  # noqa: F821
        output_scores=snakemake.output.scores_tsv,  # type: ignore[name-defined]  # noqa: F821
        tcrdock_dir=snakemake.params.tcrdock_dir,  # type: ignore[name-defined]  # noqa: F821
        alphafold_params_dir=snakemake.params.alphafold_params_dir,  # type: ignore[name-defined]  # noqa: F821
        n_candidates=snakemake.params.n_candidates,  # type: ignore[name-defined]  # noqa: F821
        fallback_hla=snakemake.params.fallback_hla,  # type: ignore[name-defined]  # noqa: F821
        fallback_tcr=snakemake.params.fallback_tcr,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main() -> None:
    parser = argparse.ArgumentParser(
        description="Run TCRdock structural validation on top neoepitope candidates."
    )
    parser.add_argument("--predictions-tsv", required=True)
    parser.add_argument("--output-pdb", required=True)
    parser.add_argument("--output-scores", required=True)
    parser.add_argument("--tcrdock-dir", required=True)
    parser.add_argument("--alphafold-params", required=True)
    parser.add_argument("--n-candidates", type=int, default=1)
    parser.add_argument("--fallback-hla", type=json.loads,
                        default={"A": "HLA-A*02:01", "B": "HLA-B*07:02", "C": "HLA-C*07:02"})
    parser.add_argument("--fallback-tcr", type=json.loads,
                        default={"va_gene": "TRAV12-2", "ja_gene": "TRAJ21",
                                 "cdr3a": "CAVNFGGGKLI", "vb_gene": "TRBV6-5",
                                 "jb_gene": "TRBJ2-7", "cdr3b": "CASSLAGGRPEQYF"})
    args = parser.parse_args()

    run_structural_validation(
        predictions_tsv=args.predictions_tsv,
        output_pdb=args.output_pdb,
        output_scores=args.output_scores,
        tcrdock_dir=args.tcrdock_dir,
        alphafold_params_dir=args.alphafold_params,
        n_candidates=args.n_candidates,
        fallback_hla=args.fallback_hla,
        fallback_tcr=args.fallback_tcr,
    )


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
