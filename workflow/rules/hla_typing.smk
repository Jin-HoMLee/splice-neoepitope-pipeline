# =============================================================================
# Rule module: HLA typing with OptiType
# =============================================================================
#
# Runs OptiType on per-sample FASTQs to predict patient HLA-A/B/C alleles,
# then aggregates results per patient using a normal-first policy.
#
# The resulting alleles.tsv is consumed by run_mhcflurry to make predictions
# patient-specific instead of relying on the hardcoded fallback allele.
#
# Enabled via config.hla.enabled: true.
# If disabled, run_mhcflurry uses config.mhcflurry.fallback_alleles.
#
# =============================================================================

import os


if config.get("hla", {}).get("enabled", False):

    _HLA_TYPING_DIR = os.path.join(os.path.dirname(OUT["alignment"]), "hla_typing")
    _OPTITYPE_SOLVER = str(config.get("hla", {}).get("solver", "cbc")).lower()
    _OPTITYPE_ILP_THREADS = int(config.get("hla", {}).get("ilp_threads", 1))

    if _OPTITYPE_SOLVER not in {"glpk", "cbc", "cplex"}:
        raise ValueError(
            f"Unsupported HLA solver '{_OPTITYPE_SOLVER}'. "
            "Set config.hla.solver to one of: glpk, cbc, cplex."
        )


    def _hla_sample_fastqs(wildcards):
        """Return list of FASTQ path(s) for a given sample."""
        for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id):
            if s["sample_id"] == wildcards.sample:
                fastqs = [s["fastq1"]]
                if (s.get("fastq2") or "").strip():
                    fastqs.append(s["fastq2"])
                return fastqs
        return []


    def _aggregate_hla_inputs(wildcards):
        sample_ids = [
            s["sample_id"]
            for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id)
        ]
        return {
            "result_tsvs": expand(
                os.path.join(
                    _HLA_TYPING_DIR, wildcards.patient_id,
                    "{sample}", "{sample}_result.tsv",
                ),
                sample=sample_ids,
            ),
        }


    rule run_optitype:
        """Run OptiType HLA typing on a single sample FASTQ.

        OptiType aligns reads to an HLA-specific reference using razers3,
        then solves an ILP to call the most likely HLA-A/B/C genotype.

        Uses --rna mode for RNA-seq data (lenient edit-distance tolerance).
        A config.ini is generated at runtime so razers3 is resolved from the
        active conda environment rather than relying on a hardcoded path.

        If OptiType produces no result (too few HLA reads), an empty TSV
        with the correct header is written so aggregate_hla_alleles can
        apply the configured fallback alleles downstream.
        """
        input:
            fastqs=_hla_sample_fastqs,
        output:
            result_tsv=os.path.join(
                _HLA_TYPING_DIR, "{patient_id}", "{sample}", "{sample}_result.tsv"
            ),
            coverage_plot=os.path.join(
                _HLA_TYPING_DIR, "{patient_id}", "{sample}", "{sample}_coverage_plot.pdf"
            ),
        log:
            os.path.join(OUT["logs"], "hla_typing", "{patient_id}_{sample}.log"),
        params:
            outdir=lambda w: os.path.join(_HLA_TYPING_DIR, w.patient_id, w.sample),
            solver=_OPTITYPE_SOLVER,
            ilp_threads=_OPTITYPE_ILP_THREADS,
        threads: config.get("hla", {}).get("threads", 4)
        conda:
            "../envs/optitype.yaml"
        shell:
            r"""
            set -uo pipefail
            mkdir -p {params.outdir} $(dirname {log})

            # Generate a config.ini that resolves razers3 from the conda env PATH.
            # OptiType needs this to find the aligner at runtime.
            RAZERS3=$(which razers3 2>/dev/null || echo "razers3")
            cat > {params.outdir}/config.ini << EOF
[mapping]
razers3=$RAZERS3
threads={threads}

[ilp]
solver={params.solver}
threads={params.ilp_threads}

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
EOF

            echo "[run_optitype] Sample: {wildcards.sample}" | tee -a {log}
            echo "[run_optitype] Input FASTQs: {input.fastqs}" | tee -a {log}
            echo "[run_optitype] Solver: {params.solver} (ILP threads={params.ilp_threads}, mapping threads={threads})" | tee -a {log}

            PYTHONUNBUFFERED=1 OptiTypePipeline.py \
                -i {input.fastqs} \
                --rna \
                --prefix {wildcards.sample} \
                --outdir {params.outdir} \
                -c {params.outdir}/config.ini \
                >> {log} 2>&1 \
                || echo "[run_optitype] WARN: OptiType failed — fallback will be used" | tee -a {log}

            # Guarantee declared outputs exist so Snakemake is satisfied.
            # aggregate_hla_alleles handles empty TSVs via the fallback path.
            if [ ! -f "{output.result_tsv}" ]; then
                printf '\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n' > {output.result_tsv}
                echo "[run_optitype] WARN: no result TSV produced — wrote empty header" | tee -a {log}
            fi
            if [ ! -f "{output.coverage_plot}" ]; then
                touch {output.coverage_plot}
            fi
            """


    rule aggregate_hla_alleles:
        """Aggregate per-sample OptiType results into a per-patient alleles TSV.

        Applies a normal-first policy: HLA alleles are germline, so calls from
        Solid Tissue Normal / Blood Derived Normal are preferred over tumor
        (which may carry LOH at the HLA locus). Falls back to
        config.mhcflurry.fallback_alleles when no sample passes the min_reads threshold.

        Writes a side hla_qc.tsv flagging the source of each locus call
        (normal / tumor / fallback) and any normal/tumor discrepancies.
        """
        input:
            unpack(_aggregate_hla_inputs),
        output:
            alleles=os.path.join(_HLA_TYPING_DIR, "{patient_id}", "alleles.tsv"),
            qc=os.path.join(_HLA_TYPING_DIR, "{patient_id}", "hla_qc.tsv"),
        log:
            os.path.join(OUT["logs"], "hla_typing", "aggregate_{patient_id}.log"),
        params:
            samples_tsv=config["samples_tsv"],
            min_reads=config.get("hla", {}).get("min_reads_per_locus", 30),
            fallback=config.get("mhcflurry", {}).get("fallback_alleles", {}),
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/aggregate_hla_alleles.py"
