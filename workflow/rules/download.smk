# =============================================================================
# Rule module: Step 1 — Download splice-junction quantification data from GDC
# =============================================================================

rule download_gdc_manifest:
    """Query the GDC API and build a manifest of splice-junction files for one
    cancer type.  The manifest is a TSV with columns: file_id, file_name,
    sample_type, project_id."""
    output:
        manifest=os.path.join(OUT["raw_data"], "{cancer_type}", "manifest.tsv"),
    log:
        os.path.join(OUT["logs"], "download", "{cancer_type}_manifest.log"),
    params:
        gdc_files_endpoint=config["gdc"]["files_endpoint"],
        max_files=config["gdc"]["max_files_per_project"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/download_gdc_data.py"


checkpoint download_gdc_files:
    """Download every file listed in the manifest.  Uses a checkpoint so
    Snakemake re-evaluates the DAG once the exact set of downloaded files is
    known.

    **Authentication**: TCGA data is controlled-access and requires a GDC
    token. See README.md section 'GDC Authentication' for instructions.
    Set `gdc.token_file` in config/config.yaml to the path of your token file.
    """
    input:
        manifest=rules.download_gdc_manifest.output.manifest,
    output:
        done=os.path.join(OUT["raw_data"], "{cancer_type}", "download.done"),
        data_dir=directory(os.path.join(OUT["raw_data"], "{cancer_type}", "files")),
    log:
        os.path.join(OUT["logs"], "download", "{cancer_type}_files.log"),
    params:
        gdc_data_endpoint=config["gdc"]["data_endpoint"],
        gdc_token_file=config["gdc"].get("token_file", None),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/download_gdc_data.py"
