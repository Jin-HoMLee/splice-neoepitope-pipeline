# =============================================================================
# Splice Neoepitope Pipeline — Main Snakefile
# =============================================================================
#
# Modernised reimplementation of the 2015 neoepitope prediction pipeline
# (Jin-Ho Lee, Seoul National University).
#
# Aligner options
# ───────────────
#   "star"   — Full accuracy, requires ~32 GB RAM
#   "hisat2" — Lower memory (~8 GB), good for laptops/small servers
#
# Workflow steps
# ──────────────
#   1. rna-align            — align FASTQ files using STAR or HISAT2
#   2. junction-filter      — classify junctions by origin (tumor_exclusive / normal_shared)
#   3. hla-typing           — patient HLA allele prediction with OptiType (optional)
#   4. contig-assemble      — build junction-spanning contigs from tumor-exclusive junctions
#   5. peptide-translate    — in-silico translation into junction-spanning peptides
#   5b. proteome-filter     — remove self-peptides present in UniProt Swiss-Prot (optional)
#   6. mhc-affinity         — MHC-I binding affinity prediction with MHCflurry
#   7. TCR-pMHC-structure   — TCR–pMHC complex structure prediction with TCRdock (optional, GPU)
#   8. report               — HTML report: junction summary, HLA QC, top binder candidates
#
# Usage
# ──────
#   snakemake --cores <N> --use-conda
#
# Configuration
# ─────────────
#   Edit config/config.yaml to set the aligner and other parameters.
#
# =============================================================================

configfile: "config/config.yaml"

# ── include rule modules ─────────────────────────────────────────────────────
include: "workflow/rules/common.smk"
include: "workflow/rules/download.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/hla_typing.smk"
include: "workflow/rules/filter_junctions.smk"
include: "workflow/rules/assemble_contigs.smk"
include: "workflow/rules/translate_peptides.smk"
include: "workflow/rules/proteome_filter.smk"
include: "workflow/rules/mhc_affinity.smk"
include: "workflow/rules/analysis.smk"
include: "workflow/rules/structure.smk"

# When TCRdock is enabled, prefer the structure report over the plain report.
if config.get("tcrdock", {}).get("enabled", False):
    ruleorder: generate_report_with_structure > generate_report

# ── final target ─────────────────────────────────────────────────────────────
# _read_samples_tsv (common.smk) validates exactly one patient per run.
PATIENT_ID = _read_samples_tsv(config["samples_tsv"])[0]["patient_id"]

rule all:
    input:
        f"{_RES}/{PATIENT_ID}/reports/report.html",
