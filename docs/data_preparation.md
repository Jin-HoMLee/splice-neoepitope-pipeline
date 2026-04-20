# Data Preparation

## Aligner Selection

The pipeline supports two RNA-seq aligners. Set `alignment.aligner` in
`config/config.yaml`:

| Aligner | RAM Required | Index Size | Best For |
|---------|--------------|------------|----------|
| **HISAT2** (default) | ~8 GB | ~8 GB | Laptops, small servers, cloud CPU VMs |
| **STAR** | ~32 GB | ~30 GB | Full accuracy, high-memory systems |

```yaml
alignment:
  aligner: "hisat2"    # or "star"
```

HISAT2 is the default and is sufficient for most use cases. STAR is available
for production runs where maximum sensitivity is required.

---

## Obtaining RNA-Seq FASTQs

### ENA (recommended — no account required)

ENA provides direct HTTPS FASTQ download for SRA accessions. No `sra-tools`
needed:

```bash
# Example: gastric cancer tumor (SRR9143066)
curl -L "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9143066/SRR9143066.fastq.gz" \
    -o data/SRR9143066.fastq.gz
```

The path pattern is:
`https://ftp.sra.ebi.ac.uk/vol1/fastq/<SRR_prefix6>/<SRR_prefix3>/<accession>/<accession>.fastq.gz`

For paired-end: append `_1.fastq.gz` / `_2.fastq.gz`.

### Other public sources

| Source | Description |
|--------|-------------|
| **GEO** | Gene Expression Omnibus |
| **ENCODE** | High-quality RNA-seq from cell lines |
| **GTEx** | Normal tissue RNA-seq (open access) |

### Remote paths in the sample manifest

`gs://` and `https://` paths in `fastq1`/`fastq2` are downloaded automatically
by the pipeline before alignment. No manual staging needed.

---

## Sample Manifest

Each patient has its own TSV under `config/samples/`. Pass it at runtime:

```bash
snakemake ... --config samples_tsv=config/samples/patient_001.tsv
```

### Format

```tsv
patient_id	sample_id	sample_type	fastq1	fastq2
patient_001	SRR9143066	Primary Tumor	https://ftp.sra.ebi.ac.uk/.../SRR9143066.fastq.gz
patient_001	SRR9143065	Solid Tissue Normal	https://ftp.sra.ebi.ac.uk/.../SRR9143065.fastq.gz
```

| Column | Required | Description |
|--------|----------|-------------|
| `patient_id` | Yes | Patient identifier — all rows with the same ID are treated as a matched set |
| `sample_id` | Yes | Unique sample identifier |
| `sample_type` | Yes | `"Primary Tumor"`, `"Solid Tissue Normal"`, or `"Blood Derived Normal"` |
| `fastq1` | Yes | Read 1 FASTQ — local path, `gs://` URI, or `https://` URL |
| `fastq2` | No | Read 2 FASTQ for paired-end data; leave empty for single-end |

### Sample type roles

| Sample type | Junction filtering | HLA typing |
|-------------|-------------------|------------|
| `Primary Tumor` | Source of candidate junctions | Yes |
| `Solid Tissue Normal` | Filters out junctions present in normal | Yes |
| `Blood Derived Normal` | Not used for junction filtering | Yes (preferred for germline HLA calls) |

When no normal sample is present, all unannotated junctions are labelled
`tumor_exclusive` with a warning — the pipeline still runs.
