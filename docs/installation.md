# Installation

## System Requirements

| Requirement | Minimum | Notes |
|-------------|---------|-------|
| OS | Linux (x86-64) or macOS | Windows is not supported |
| CPU | 4 cores | More cores speed up alignment and prediction |
| RAM | **8 GB** | Using HISAT2 (default); 32 GB for STAR |
| Disk | 50 GB free | Reference genome + data files |
| Python | 3.11+ | Managed automatically via conda |
| Git | any | For cloning the repository |

> For contributors: best-practice local GitHub CLI authentication (Keychain-primary, fine-grained tokens) is documented in [gh auth setup](gh_auth_setup.md).

---

## 1. Install Conda (Miniforge — recommended)

Miniforge is a minimal Conda installer that defaults to `conda-forge` and ships
with the faster `mamba` solver.

```bash
# Linux
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
  -o Miniforge3.sh

# macOS (Apple Silicon)
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh \
  -o Miniforge3.sh

# macOS (Intel)
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh \
  -o Miniforge3.sh

bash Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init bash   # or: conda init zsh
```

> **Already have Anaconda or Miniconda?** That works too. Make sure the channels
> are set:
> ```bash
> conda config --add channels bioconda
> conda config --add channels conda-forge
> conda config --set channel_priority strict
> ```

---

## 2. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda \
  "snakemake>=8.0,<9" \
  python=3.11 \
  -y

conda activate snakemake
snakemake --version   # expected: 8.x.x
```

Rule-specific environments (`workflow/envs/*.yaml`) are created automatically on
first use — you do **not** need to install biopython, bedtools, HISAT2, etc. manually.

---

## 3. Clone the Repository

```bash
git clone https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git
cd splice-neoepitope-pipeline
```

---

## 4. Download Reference Data

Follow the instructions in [`docs/resources.md`](resources.md) to download
the GRCh38 FASTA and GENCODE v47 GTF into `references/`. These are required
for the full pipeline; the chr22 test subset is handled automatically by
`scripts/prepare_test_data.sh`.

---

## 5. Sanity Check

```bash
conda activate snakemake
snakemake --cores 1 --use-conda -n \
    --configfile config/config.yaml \
    --config samples_tsv=config/samples/patient_001.tsv
```

You should see a list of planned jobs (or "Nothing to be done" if outputs already
exist). A Python traceback indicates a missing config file or inactive conda env.

---

## 6. DAG Visualization (dev tool)

After editing Snakemake rules, render the rule graph to confirm the DAG matches
intent — `scripts/setup_local.sh` installs [snakevision](https://pypi.org/project/snakevision/)
into the `snakemake` env automatically. Wrapper:

```bash
conda activate snakemake
bash scripts/visualize_dag.sh                                                       # default: config/test_config.yaml
bash scripts/visualize_dag.sh config/config.yaml --config samples_tsv=config/samples/patient_002.tsv
```

The wrapper runs `snakemake --rulegraph --forceall`, pipes the `.dot` to
`snakevision`, writes `dag.svg`, and (on macOS) opens it. Override the output
path with `OUTPUT=foo.svg bash scripts/visualize_dag.sh`.

**`--clean` flag — full per-sample DAG.** Snakemake's `--rulegraph` prunes
rules whose outputs already exist, even with `--forceall`. Once the test
pipeline has run once, per-sample rules (`run_optitype`, `hisat2_align`,
`extract_junctions`, `download_fastq`, …) drop off the rendered DAG. To see
the complete pipeline architecture, use `--clean`:

```bash
bash scripts/visualize_dag.sh --clean
```

This renders against a symlink-only temp workspace (no `results/`), so snakemake
sees a "no outputs yet" state and emits every rule. Real `results/` are untouched.

> **Caveat:** `--clean` reads `samples_tsv` from the configfile only. Overrides
> via `--config samples_tsv=…` extra args are passed through to snakemake but
> NOT applied to placeholder FASTQ generation, so the DAG can come out pruned or
> fail. Use a configfile that already points to the desired samples TSV (the
> wrapper warns at runtime if it detects `--config samples_tsv=…`).

Recommended workflow when changing rules: edit → dry-run → `bash scripts/visualize_dag.sh`
→ open the SVG → confirm intent before opening the PR. Use `--clean` when adding
or restructuring per-sample rules.
