# =============================================================================
# Rule module: Steps 1b + 2 — Reference junction construction and novel-junction
# filtering
# =============================================================================

import csv
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


def _junction_files_gdc(wildcards):
    """Return the list of raw junction quantification files for this cancer
    type, discovered after the GDC download checkpoint resolves."""
    checkpoint_output = checkpoints.download_gdc_files.get(**wildcards).output.data_dir
    return glob_wildcards(
        os.path.join(checkpoint_output, "{file_id}.tsv")
    ).file_id


def _junction_files_local(wildcards):
    """Return the list of junction files for local alignment mode."""
    samples_file = config.get("samples_tsv")
    if not samples_file:
        return []
    samples_path = Path(samples_file)
    if not samples_path.exists():
        return []
    sample_ids = []
    with samples_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if not row["sample_id"].startswith("#"):
                sample_ids.append(row["sample_id"])
    return sample_ids


def _get_junction_files_input(wildcards):
    """Get junction files based on data source mode."""
    if config.get("data_source") == "local":
        sample_ids = _junction_files_local(wildcards)
        return expand(
            os.path.join(OUT["raw_data"], wildcards.patient_id, "files", "{sample_id}.tsv"),
            sample_id=sample_ids,
        )
    else:
        file_ids = _junction_files_gdc(wildcards)
        return expand(
            os.path.join(OUT["raw_data"], wildcards.patient_id, "files", "{file_id}.tsv"),
            file_id=file_ids,
        )


def _get_manifest_input(wildcards):
    """Get manifest file based on data source mode."""
    return os.path.join(OUT["raw_data"], wildcards.patient_id, "manifest.tsv")


rule filter_junctions:
    """For each raw junction quantification file:
      1. Keep junctions with read count > mean read count in that file
         (removes low-read background noise).
      2. Remove junctions present in the reference (GENCODE) junction list.
    The output is a TSV of novel junctions per sample.

    Works with both GDC downloads and local STAR alignments."""
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
