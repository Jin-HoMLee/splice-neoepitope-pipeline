# =============================================================================
# Splice Neoepitope Pipeline — Main Snakefile
# =============================================================================
#
# Modernised reimplementation of the 2015 neoepitope prediction pipeline
# (Jin-Ho Lee, Seoul National University).
#
# Data source modes
# ─────────────────
#   "gdc"   — Download pre-computed junction files from GDC (requires dbGaP access)
#   "local" — Align your own FASTQ files using STAR or HISAT2 (open access)
#
# Aligner options (for local mode)
# ────────────────────────────────
#   "star"   — Full accuracy, requires ~32 GB RAM (default)
#   "hisat2" — Lower memory (~8 GB), good for laptops/small servers
#
# Workflow steps
# ──────────────
#   1. download/align — fetch TCGA data via GDC API OR align local FASTQ
#   2. filter         — classify junctions by origin (tumor_specific / patient_specific)
#   3. assemble       — build 50 nt contigs around tumor_specific junctions
#   4. translate      — in-silico translation into 16-mer peptides (3 reading frames)
#   5. predict        — MHCflurry 2.x epitope prediction
#   6. report         — junction origin summary + top binders HTML report
#
# Usage
# ──────
#   snakemake --cores <N> --use-conda
#
# Configuration
# ─────────────
#   Edit config/config.yaml to set data_source mode and other parameters.
#
# =============================================================================

configfile: "config/config.yaml"

# ── convenience aliases ──────────────────────────────────────────────────────
DATA_SOURCE = config.get("data_source", "local")
PATIENT_ID  = config["patient_id"]
OUT         = config["output"]


# ── include rule modules ─────────────────────────────────────────────────────
include: "workflow/rules/download.smk"
include: "workflow/rules/local_alignment.smk"
include: "workflow/rules/hisat2_alignment.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/assemble.smk"
include: "workflow/rules/translate.smk"
include: "workflow/rules/predict.smk"
include: "workflow/rules/analysis.smk"
include: "workflow/rules/tcrdock.smk"

# When TCRdock is enabled, prefer the structure report over the plain report.
if config.get("tcrdock", {}).get("enabled", False):
    ruleorder: generate_report_with_structure > generate_report


# ── final target ─────────────────────────────────────────────────────────────
rule all:
    input:
        f"{OUT['reports']}/{PATIENT_ID}/report.html",

