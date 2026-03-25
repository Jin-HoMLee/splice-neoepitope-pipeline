# =============================================================================
# Rule module: Steps 1b + 2 — Reference junction construction and novel-junction
# filtering
# =============================================================================

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


def _junction_files(wildcards):
    """Return the list of raw junction quantification files for this cancer
    type, discovered after the download checkpoint resolves."""
    checkpoint_output = checkpoints.download_gdc_files.get(**wildcards).output.data_dir
    return glob_wildcards(
        os.path.join(checkpoint_output, "{file_id}.tsv")
    ).file_id


rule filter_junctions:
    """For each raw junction quantification file:
      1. Keep junctions with read count > mean read count in that file
         (removes low-read background noise).
      2. Remove junctions present in the reference (GENCODE) junction list.
    The output is a TSV of novel junctions per sample."""
    input:
        junction_files=lambda wc: expand(
            os.path.join(OUT["raw_data"], wc.cancer_type, "files", "{file_id}.tsv"),
            file_id=_junction_files(wc),
        ),
        manifest=os.path.join(OUT["raw_data"], "{cancer_type}", "manifest.tsv"),
        reference_junctions=config["reference"]["junction_bed"],
    output:
        novel_junctions=os.path.join(
            OUT["junctions"], "{cancer_type}", "novel_junctions.tsv"
        ),
    log:
        os.path.join(OUT["logs"], "filter", "{cancer_type}_filter.log"),
    params:
        strategy=config["filtering"]["strategy"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_junctions.py"
