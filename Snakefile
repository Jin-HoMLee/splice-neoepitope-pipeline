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
#   1. align    — align FASTQ files using STAR or HISAT2
#   2. filter   — classify junctions by origin (tumor_exclusive / normal_shared)
#   3. hla      — OptiType HLA typing (optional)
#   4. assemble — build 50 nt contigs around tumor_exclusive junctions
#   5. translate — in-silico translation into junction-spanning 9-mers
#   6. predict  — MHCflurry 2.x epitope prediction
#   7. report   — junction origin summary + HLA QC + top binders HTML report
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

import csv
import os

# ── convenience aliases ──────────────────────────────────────────────────────
# All patient outputs are rooted at _RES/{patient_id}/<step>/
_RES  = config["output"]["results_dir"]
_LOGS = config["output"]["logs"]

# Derive patient IDs from the samples TSV — the TSV is the single source of
# truth for what to run. All rows sharing the same patient_id are treated as
# a matched set (e.g. tumor + normal for one patient).
def _read_patient_ids(samples_tsv):
    patient_ids = set()
    with open(samples_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = (row.get("patient_id") or "").strip()
            if not pid or pid.startswith("#"):
                continue
            patient_ids.add(pid)
    return sorted(patient_ids)

PATIENT_IDS = _read_patient_ids(config["samples_tsv"])

wildcard_constraints:
    patient_id="[^/]+",
    sample="[^/]+",


# ── include rule modules ─────────────────────────────────────────────────────
include: "workflow/rules/common.smk"
include: "workflow/rules/download.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/hla_typing.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/assemble.smk"
include: "workflow/rules/translate.smk"
include: "workflow/rules/mhcflurry.smk"
include: "workflow/rules/analysis.smk"
include: "workflow/rules/tcrdock.smk"

# When TCRdock is enabled, prefer the structure report over the plain report.
if config.get("tcrdock", {}).get("enabled", False):
    ruleorder: generate_report_with_structure > generate_report


# ── final target ─────────────────────────────────────────────────────────────
rule all:
    input:
        expand(os.path.join(_RES, "{patient_id}", "reports", "report.html"), patient_id=PATIENT_IDS),

