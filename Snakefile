# =============================================================================
# Splice Neoepitope Pipeline — Main Snakefile
# =============================================================================
#
# Modernised reimplementation of the 2015 neoepitope prediction pipeline
# (Jin-Ho Lee, Seoul National University).
#
# Workflow steps
# ──────────────
#   1. download   – fetch TCGA splice-junction quantification data via GDC API
#   2. filter     – remove low-read and known (reference) junctions
#   3. assemble   – build 50 nt contigs around each novel junction
#   4. translate  – in-silico translation into 16-mer peptides (3 reading frames)
#   5. predict    – NetMHCPan 4.1 epitope prediction
#   6. analysis   – Fisher's exact test, summary statistics, report generation
#
# Usage
# ──────
#   snakemake --cores <N> --use-conda
#
# Configuration
# ─────────────
#   Edit config/config.yaml to change cancer types, file paths, thresholds, etc.
#
# =============================================================================

configfile: "config/config.yaml"

# ── convenience aliases ──────────────────────────────────────────────────────
CANCER_TYPES = config["cancer_types"]
OUT           = config["output"]

# ── include rule modules ─────────────────────────────────────────────────────
include: "workflow/rules/download.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/assemble.smk"
include: "workflow/rules/translate.smk"
include: "workflow/rules/predict.smk"
include: "workflow/rules/analysis.smk"


# ── final target ─────────────────────────────────────────────────────────────
rule all:
    input:
        # Per-cancer type final report
        expand(
            "{reports}/{cancer_type}/report.html",
            reports=OUT["reports"],
            cancer_type=CANCER_TYPES,
        ),
        # Cross-cancer summary table
        os.path.join(OUT["reports"], "summary_table.tsv"),
