# Configuration Reference

All pipeline parameters live in `config/config.yaml`. For TCRdock (GPU step),
an overlay file `config/tcrdock_gpu.yaml` is merged at runtime.

---

## Reference Genome

```yaml
reference:
  genome_fasta: "resources/GRCh38.primary_assembly.genome.fa"
  gencode_gtf:  "resources/gencode.v47.annotation.gtf.gz"
  junction_bed: "resources/reference_junctions.bed"   # auto-generated
```

See [`resources/README.md`](../resources/README.md) for download commands.

---

## Alignment

```yaml
alignment:
  aligner: "hisat2"                           # "hisat2" (8 GB) or "star" (32 GB)
  hisat2_index_dir: "resources/hisat2_index"
  star_index_dir:   "resources/star_index"
  threads: 8
```

See [`docs/data_preparation.md`](data_preparation.md) for aligner comparison.

---

## Junction Filtering

```yaml
filtering:
  min_normal_reads: 2   # min reads in normal to classify a junction as normal_shared
```

---

## Contig Assembly

```yaml
assembly:
  upstream_nt:   26   # nt upstream of the splice junction
  downstream_nt: 24   # nt downstream of the splice junction
  contig_length: 50   # total contig length
```

---

## HLA Typing (OptiType)

```yaml
hla:
  enabled: true
  min_reads_per_locus: 30   # min OptiType read count to trust a call
  threads: 5                # razers3 mapping threads per sample; ≥5 forces sequential
                            # execution on 8-core VMs to prevent concurrent OOM
  solver: "cbc"             # ILP solver: cbc (recommended), glpk, or cplex
  ilp_threads: 4
```

When `hla.enabled: false`, the fallback alleles below are used for all predictions.

---

## Epitope Prediction (MHCflurry)

```yaml
mhcflurry:
  fallback_alleles:         # used when HLA typing is disabled or a locus has too few reads
    A: "HLA-A*02:01"
    B: "HLA-B*07:02"
    C: "HLA-C*07:02"
  ic50_strong:  50          # nM — strong binder threshold
  ic50_weak:   500          # nM — weak binder threshold
  threads: 8
  prediction_mode: "affinity"   # "affinity" | "presentation" | "processing"
```

---

## TCRdock Structural Validation (GPU)

Disabled by default. Enable via the `config/tcrdock_gpu.yaml` overlay:

```bash
snakemake --configfile config/config.yaml config/tcrdock_gpu.yaml ...
```

Key parameters (in `config/tcrdock_gpu.yaml`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tcrdock.enabled` | `true` | Enables the TCRdock step |
| `tcrdock.docker_image` | `tcrdock:latest` | Docker image built by `setup_tcrdock_vm.sh` |
| `tcrdock.n_candidates` | `1` | Number of top candidates to model |
| `tcrdock.fallback_tcr` | DMF5 TCR | TCR sequences used until TRUST4 integration (#24) |

See [`docs/google_cloud_guide.md`](google_cloud_guide.md) for GPU VM setup.
