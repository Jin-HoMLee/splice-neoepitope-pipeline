#!/usr/bin/env bash
# scripts/visualize_dag.sh
#
# Render a Snakemake rule-graph SVG via snakevision and open it.
# Workflow: snakemake --rulegraph → .dot file → snakevision → SVG → open.
#
# snakevision uses dagviz's "metro map" style: top-to-bottom data flow with
# horizontal lane branching (subway-diagram metaphor). The defaults below
# (scale, node_radius, edge_stroke_width) are tuned for readability on
# multi-rule pipelines — override via STYLE_ARGS env if needed.
#
# Usage:
#   bash scripts/visualize_dag.sh [--clean] [configfile] [extra snakemake args...]
#
# Examples:
#   bash scripts/visualize_dag.sh                                                          # default: config/test_config.yaml
#   bash scripts/visualize_dag.sh --clean                                                  # full DAG (per-sample rules included)
#   bash scripts/visualize_dag.sh config/test_config.yaml
#   bash scripts/visualize_dag.sh config/config.yaml --config samples_tsv=config/samples/patient_002.tsv
#
# --clean flag: snakemake's --rulegraph prunes rules whose outputs already exist,
# even with --forceall. Per-sample rules (run_optitype, hisat2_align, etc.) drop
# off the DAG once the test pipeline has been run. --clean renders against a
# symlink-only temp workspace (all dirs except results/ are symlinked from $PWD)
# so snakemake sees a "no outputs yet" state and emits the full per-sample DAG.
#
# Override styling (full list: dagviz.style.metro.StyleConfig fields):
#   STYLE_ARGS="--style scale=20 --style node_radius=12" bash scripts/visualize_dag.sh
#
# Prerequisite: snakevision installed in the active conda env. setup_local.sh installs
# it into the 'snakemake' env; manually: `conda activate snakemake && pip install snakevision`.
#
set -euo pipefail

CLEAN_MODE=false
if [[ "${1:-}" == "--clean" ]]; then
    CLEAN_MODE=true
    shift
fi

CONFIGFILE="${1:-config/test_config.yaml}"
shift || true
OUTPUT="${OUTPUT:-dag.svg}"
STYLE_ARGS="${STYLE_ARGS:---style scale=15 --style node_radius=10 --style edge_stroke_width=3}"

DOT_FILE=$(mktemp)
TMPWORK=""

cleanup() {
    rm -f "$DOT_FILE"
    [[ -n "$TMPWORK" && -d "$TMPWORK" ]] && rm -rf "$TMPWORK"
}
trap cleanup EXIT

if ! command -v snakevision &>/dev/null; then
    echo "Error: snakevision not found on PATH." >&2
    echo "Activate the snakemake env and install: conda activate snakemake && pip install snakevision" >&2
    exit 1
fi

if [[ ! -f "$CONFIGFILE" ]]; then
    echo "Error: configfile '$CONFIGFILE' not found." >&2
    exit 1
fi

if [[ "$CLEAN_MODE" == true ]]; then
    echo "Generating rule graph in clean workspace (full DAG, all per-sample rules)..."
    TMPWORK=$(mktemp -d)
    # Symlink read-only project files; omit results/ + data/ so snakemake sees a
    # fresh state. data/ is populated below with empty placeholders for source
    # FASTQs so input-existence checks pass without requiring real test data.
    for f in Snakefile workflow config scripts resources .snakemake; do
        [[ -e "$f" ]] && ln -s "$PWD/$f" "$TMPWORK/$f"
    done

    # Resolve samples_tsv from the configfile (override via --config samples_tsv=...
    # passed in "$@" is intentionally ignored here for simplicity — pass the right
    # configfile explicitly if you need a different samples list).
    SAMPLES_TSV=$(awk -F'[: ]+' '/^samples_tsv:/ {gsub(/"/, "", $2); print $2; exit}' "$CONFIGFILE")
    if [[ -n "$SAMPLES_TSV" && -f "$SAMPLES_TSV" ]]; then
        # For each non-URL FASTQ path in the samples TSV, create an empty placeholder
        # at the same relative path inside TMPWORK. URLs (gs://, https://) are produced
        # by the download_fastq rule and don't need placeholders.
        awk -F'\t' 'NR>1 && $1 !~ /^#/ {
            for (i=4; i<=5; i++) {
                if ($i != "" && $i !~ /^(https?|gs):\/\//) print $i
            }
        }' "$SAMPLES_TSV" | while IFS= read -r fq; do
            [[ -z "$fq" ]] && continue
            mkdir -p "$TMPWORK/$(dirname "$fq")"
            : > "$TMPWORK/$fq"
        done
    fi

    (cd "$TMPWORK" && snakemake --configfile "$CONFIGFILE" "$@" --rulegraph --forceall) > "$DOT_FILE"
else
    echo "Generating rule graph (configfile: $CONFIGFILE; pass --clean for full per-sample DAG)..."
    snakemake --configfile "$CONFIGFILE" "$@" --rulegraph --forceall > "$DOT_FILE"
fi

echo "Rendering with snakevision (metro style, top-to-bottom flow)..."
# shellcheck disable=SC2086 # STYLE_ARGS is intentionally word-split
snakevision $STYLE_ARGS -o "$OUTPUT" "$DOT_FILE"

echo "Wrote $OUTPUT"
if [[ "$OSTYPE" == "darwin"* ]]; then
    open "$OUTPUT"
fi
