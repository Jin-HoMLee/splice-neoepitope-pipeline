# =============================================================================
# Rule module: Steps 1b + 2 — Reference junction construction and novel-junction
# filtering
# =============================================================================

import csv
import os
from pathlib import Path


rule build_reference_junctions:
    """Parse the GENCODE GRCh38 GTF and extract all annotated donor/acceptor
    splice-junction pairs.  Output is a sorted BED-like file:
      chrom  start  end  strand
    where start/end are the 0-based half-open coordinates of the junction gap
    (i.e., last nt of exon N to first nt of exon N+1)."""
    input:
        gtf=config["reference"]["gencode_gtf"],
    output:
        bed=config["reference"]["junction_bed"],
    log:
        os.path.join(_LOGS, "filter_junctions", "build_reference_junctions.log"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/build_reference_junctions.py"


def _get_junction_files_input(wildcards):
    """Return junction TSV paths for all samples of this patient."""
    samples_file = config.get("samples_tsv")
    if not samples_file:
        return []
    samples_path = Path(samples_file)
    if not samples_path.exists():
        return []
    sample_ids = []
    with samples_path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = (row.get("patient_id") or "").strip()
            sid = (row.get("sample_id") or "").strip()
            if pid == wildcards.patient_id and sid and not sid.startswith("#"):
                sample_ids.append(sid)
    return expand(
        os.path.join(_RES, wildcards.patient_id, "alignment", "{sample_id}", "raw_junctions.tsv"),
        sample_id=sample_ids,
    )


def _get_manifest_input(wildcards):
    """Return the manifest TSV path for this patient."""
    return os.path.join(_RES, wildcards.patient_id, "alignment", "manifest.tsv")


def _get_gtex_bed_input(wildcards):
    """GTEx pan-tissue blacklist BED for the filter step, or [] when disabled.

    Returns a list so an empty result cleanly omits the optional named input
    (the established repo idiom — see _get_junction_files_input). The path is
    single-sourced via _gtex_blacklist_bed() in common.smk so it matches the
    download_gtex_pan_tissue_bed rule's output exactly. When the BED is a gs://
    object the returned local path is produced by that download rule; when it is
    a local fixture (chr22 tests) the path is consumed in place."""
    bed = _gtex_blacklist_bed()
    return [bed] if bed else []


rule filter_junctions:
    """For each raw junction quantification file:
      1. Keep junctions with read count > mean read count in that file
         (removes low-read background noise).
      2. Discard junctions present in the reference (GENCODE) junction list.
      3. Classify unannotated junctions by origin: normal_shared (matched
         normal) → gtex_pantissue_shared (GTEx pan-tissue population blacklist,
         when gtex_filter.enabled) → tumor_exclusive. Only tumor_exclusive
         junctions proceed to neoepitope prediction; the others are retained in
         the TSV for transparency.
      4. Annotate each junction with its canonical CDS reading frame derived from
         the GENCODE GTF (informational; all three frames are still translated
         downstream).
    The output is a TSV of novel junctions per sample.

    Works with both HISAT2 and STAR alignments."""
    input:
        junction_files=_get_junction_files_input,
        manifest=_get_manifest_input,
        reference_junctions=config["reference"]["junction_bed"],
        gencode_gtf=config["reference"]["gencode_gtf"],
        # Optional: GTEx pan-tissue population-normal blacklist (Issue #211/#212).
        # Empty (filter disabled) → omitted; the script treats absence as a no-op.
        gtex_bed=_get_gtex_bed_input,
    output:
        novel_junctions=os.path.join(
            _RES, "{patient_id}", "junctions", "novel_junctions.tsv"
        ),
        stats=os.path.join(
            _RES, "{patient_id}", "junctions", "junction_filter_stats.tsv"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "filter_junctions", "filter.log"),
    params:
        min_normal_reads=config["filtering"]["min_normal_reads"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_junctions.py"
