# SNAF - open-only install recipe (Linux / free-GPU; NOT local arm64)

SNAF (`frankligy/SNAF`, MIT) is an open-only **GO**, but it does not run on our macOS arm64 / no-Docker box: Python 3.7.12 + TensorFlow 2.3.0 are Linux-x86 only, and junction extraction uses AltAnalyze's **amd64-only** Docker image. Run it on Linux (RunPod/free-GPU) or under Docker/amd64 emulation.

## Open-only key: MHCflurry is a native flag, not a swap

```python
import snaf
snaf.initialize(
    df=junction_count_matrix,      # from AltAnalyze
    db_dir='/path/to/download',    # the 2.72 GB reference bundle
    binding_method='MHCflurry',    # <-- open-only: no software_path, no netMHCpan
)
```
Both engines register results under the same internal attribute, so nothing downstream changes. Never pass a `software_path` (that is the netMHCpan branch).

## Steps (Linux x86_64)

```bash
# 1. Env (Linux-x86 lockfile ships in the SNAF upstream repo, frankligy/SNAF)
conda env create -f snaf_linux_py37_env.yml && conda activate snaf37
mhcflurry-downloads fetch                      # open MHCflurry models

# 2. Reference bundle (2.72 GB)
curl -o download.tar.gz http://altanalyze.org/SNAF/download.tar.gz
tar -xzf download.tar.gz

# 3. Junction extraction (AltAnalyze, amd64 Docker) on >=2 BAMs
docker run --platform linux/amd64 -v $PWD:/data frankligy123/altanalyze:0.7.0.1 \
    <altanalyze bam-to-junction args>          # -> ExpressionInput/counts.original.pruned.txt

# 4. SNAF T-antigen run (see snaf.initialize above), binding_method='MHCflurry'
```

## Smoke note
No tiny fixture ships; AltAnalyze is cohort-level (needs >=2 BAMs). Minimal proof = run steps 3-4 on 2 small BAMs (e.g. the chr22 test pair SRR9143066/65). Genome-wide junction UIDs mean a chr22 subset proves the plumbing runs, not biological recall. This smoke is deferred to the leaf-B Linux/free-GPU slice.

## Openness
- SNAF: MIT. MHCflurry: Apache-2.0. AltAnalyze: open (amd64 container). Reference bundle: open.
- No non-redistributable binary is installed anywhere in this recipe (netMHCpan/TMHMM deliberately omitted).
