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
        os.path.join(OUT["logs"], "reference", "build_reference_junctions.log"),
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
        os.path.join(OUT["alignment"], wildcards.patient_id, "{sample_id}", "junctions.tsv"),
        sample_id=sample_ids,
    )


def _get_manifest_input(wildcards):
    """Return the manifest TSV path for this patient."""
    return os.path.join(OUT["alignment"], wildcards.patient_id, "manifest.tsv")


rule filter_junctions:
    """For each raw junction quantification file:
      1. Keep junctions with read count > mean read count in that file
         (removes low-read background noise).
      2. Remove junctions present in the reference (GENCODE) junction list.
    The output is a TSV of novel junctions per sample.

    Works with both HISAT2 and STAR alignments."""
    input:
        junction_files=_get_junction_files_input,
        manifest=_get_manifest_input,
        reference_junctions=config["reference"]["junction_bed"],
    output:
        novel_junctions=os.path.join(
            OUT["junctions"], "{patient_id}", "novel_junctions.tsv"
        ),
    log:
        os.path.join(OUT["logs"], "filter", "{patient_id}_filter.log"),
    params:
        min_normal_reads=config["filtering"]["min_normal_reads"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_junctions.py"
