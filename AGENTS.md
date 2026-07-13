# Splice Neoepitope Pipeline — Project Notes

## Instructions for Claude
Update this file only when something is **not derivable from the code or git history** — non-obvious decisions, known gotchas, workarounds, and infrastructure facts. Do not use it as a change log; that's what `git log` is for.

## Project Overview
Modernised reimplementation of a 2015 cancer neoepitope prediction pipeline (Jin-Ho Lee, Seoul National University). Identifies tumor-specific splice junctions from RNA-Seq data and predicts MHC-binding neoepitopes.

## Infrastructure
- **⚠️ GCP DECOMMISSIONED 2026-06-26 — $0-budget keep-alive posture.** The GCP stack (VMs `neoepitope-pipeline` + `neoepitope-orchestrator`, their disks, and the `gs://splice-neoepitope-project` bucket) was deleted to stop charges before the free-trial lapse; project data was preserved on Cloudflare R2 first ([#854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854), verified). Current compute = **local CPU-only core** (chr22, MHCflurry CPU fallback); the optional TCRdock structural validation runs opportunistically on free GPU (Kaggle/Colab). Rationale: project memory `gcp-trial-expiry-self-funded-infra`.
- **GCP/P100 operational recipe is ARCHIVED, not deleted** → [`docs/legacy/gcp_p100_setup.md`](docs/legacy/gcp_p100_setup.md): the VM specs, the NVIDIA closed-driver pin (the `-open`-fails-on-Pascal-GSP saga), the cu126/Pascal torch pin, and the zone-exhaustion history all live there, so a **funded revival is a checklist, not a reconstruction**. `scripts/run_cloud_gpu.sh` + `docs/google_cloud_guide.md` remain in-repo for that path. Forward paid-provider plan (RunPod L4/Ada + R2): [`docs/migration_runbook.md`](docs/migration_runbook.md).
- Pipeline is run with `snakemake --cores $(nproc) --use-conda --rerun-triggers mtime` inside a `tmux` session
- GitHub project board: user project #9 ("JH M Lee Lab") under user `Jin-HoMLee` — query via `gh api graphql` with `user(login: "Jin-HoMLee") { projectV2(number: 9) { ... } }` (it's a user project, not org)
- `main` branch protection: required CI checks (`pipeline-pytest`, `pipeline-snakemake-dry-run`) + squash-merge default. The "Require branches to be up to date before merging" rule was **removed 2026-05-09** — it fired on every PR cut from a worktree branch lagging `main`, even without real conflicts. Don't suggest `gh pr update-branch` or `--admin` workarounds for "branch is behind main" anymore (the rule is gone). Real file-level conflicts still need manual resolution.
- Remote scheduled/one-shot agents (board-hygiene sweeps; overnight Issue→draft-PR dispatch) run in an isolated cloud sandbox via the `/schedule` skill — see [`docs/remote_routines.md`](docs/remote_routines.md) for the sandbox env facts (allowlisted egress, no conda/snakemake, the `pip --upgrade pyyaml` trap, `claude/*`-only push) and the hardened dispatch prompt checklist. The CCR sandbox token is a *different* identity (`Claude <noreply@anthropic.com>`) with repo+project scope.

## Pipeline Design Decisions

### Junction origin classification
Tumor junctions are filtered against matched normal at the **junction** level, not the prediction level. Labels used in code/TSVs: `annotated` (in GENCODE) and `normal_shared` (also in normal) are discarded; `tumor_exclusive` (unannotated, absent in normal) goes to neoepitope prediction. With no normal sample, all unannotated junctions fall back to `tumor_exclusive` with a warning. How/why (the clinical rationale + the removed Fisher's exact test): [`docs/design/junction-origin-classification.md`](docs/design/junction-origin-classification.md).

### TCRdock via Docker
TCRdock runs in a Docker container (`docker/Dockerfile.pipeline`), not a conda env; the image bundles CUDA 11.8, cuDNN 8, Python 3.10, JAX 0.3.25, AlphaFold params, and BLAST, so the host needs only the NVIDIA Container Toolkit. Rationale (the conda-conflict history + forward-compatibility guarantee): [`docs/adr/0001-tcrdock-via-docker.md`](docs/adr/0001-tcrdock-via-docker.md).

### PDB chain relabelling
AlphaFold outputs all residues as a single chain (A). `relabel_pdb_chains()` in `run_tcrdock.py` reassigns chain IDs (A=MHC, B=peptide, C=TCR-α, D=TCR-β) using per-chain sequence lengths from TCRdock's `alphafold_setup/targets.tsv`. The report injects PDB COMPND records so Mol* displays meaningful chain names in the sequence panel instead of "Polymer 1/2/3/4".

**Authoring an experiment notebook or Quarto slide deck?** Layout, cross-experiment sharing, size bands, and the 3 deck tiers (experiment / eval / research-decision) live in [`docs/research_artifact_conventions.md`](docs/research_artifact_conventions.md).

## MHC Presentation Vocabulary

This pipeline uses **`Class1PresentationPredictor`** (MHCflurry 2.x), which scores *presentation likelihood* — a combined estimate of binding affinity + antigen processing. It is distinct from the older `Class1AffinityPredictor` (affinity-only).

Relevant output columns (use these names in code, reports, and prose):

- `presentation_class` — `strong | weak | non`
- `presentation_score`, `presentation_percentile`, `best_presentation_percentile`
- `genotype_presentation_score`
- `n_strong_alleles`

Use **"presenter" / "top presenters" / "presentation percentile"** throughout. Avoid **"binder" / "top binders" / "binding affinity threshold"** — those refer to the affinity-only predictor we do not use as the primary ranker. IC50 (`ic50_nM`) is still emitted for reference but is a secondary metric.

## Snakemake Conda Activation

Always activate the environment explicitly before invoking Snakemake:

```bash
conda activate snakemake
snakemake --cores $(nproc) --use-conda ...
```

**Do not** use `conda run -n snakemake snakemake ...` — `conda run` buffers all stdout, hiding real-time log output during the run.

**Conda env cleanup after `workflow/envs/*.yaml` changes:** Automatic cleanup was removed from `run_cloud_gpu.sh`. When any `workflow/envs/*.yaml` file changes, manually delete the affected old environment on the VM before running — old envs will not be rebuilt automatically and stale cached packages will be used instead.

## Python environments

The project has **4 functional Python environments**, each with a different use case. **No root-level `.venv` exists or should exist** — any reference to `.venv/bin/python` at the project root is wrong.

| Use case | Canonical path | Setup recipe | Documented at |
|---|---|---|---|
| Pytest test suite | `workflow/tests/.venv` | `pyenv local 3.13.5 && uv venv --seed --python 3.13.5 workflow/tests/.venv && uv pip install --python workflow/tests/.venv/bin/python -r workflow/tests/requirements-test.txt` | [workflow/tests/README.md](workflow/tests/README.md) |
| Research / Jupyter notebooks | `research/.venv` | `cd research/ && pyenv local 3.14.4 && uv venv --seed --python 3.14.4 .venv && uv pip install --python .venv/bin/python -r requirements.txt` | [research/README.md](research/README.md) |
| Snakemake orchestration (and ad-hoc `python -c` quick checks) | `conda activate snakemake` | `bash scripts/setup_local.sh` | "Snakemake Conda Activation" above |
| Per-rule Python (`workflow/scripts/*.py` invoked by Snakemake rules) | rule's own `--use-conda` env | auto-created from `workflow/envs/*.yaml` on first `snakemake --use-conda` run | implicit; each rule declares `conda: "envs/<env>.yaml"` |

**Per-clone setup.** The two pyenv venvs (`workflow/tests/.venv` and `research/.venv`) are per-clone and gitignored. After the 2026-05-14 separate-clone migration, each clone needs its own one-time setup — neither is auto-created by `scripts/setup_local.sh`. The `.python-version` files are also gitignored so clones stay independent.

**Picking an interpreter:**
- Pytest invocations: `workflow/tests/.venv/bin/python -m pytest workflow/tests/...` (never `.venv/bin/python -m pytest`).
- Ad-hoc `python -c "..."` quick checks (yaml parse, file compile-check, one-off imports): `conda activate snakemake` first, then call bare `python -c "..."` — no separate path.
- Running `workflow/scripts/*.py` directly outside Snakemake is rarely needed; if you must, activate the rule's conda env explicitly (e.g. `conda activate .snakemake/conda/<hash>` after Snakemake has built it).

The `snakemake` conda env stays pristine on purpose — adding test or notebook deps risks dep drift in the workhorse env ([workflow/tests/README.md "Why this split"](workflow/tests/README.md)).

## Snakemake 8 Gotchas

### `--configfile` flag collapsing
In Snakemake 8, passing `--configfile` as **separate flags** (`--configfile A --configfile B`) causes the second invocation to replace the first due to argparse `nargs="+"` semantics. Only the last file is loaded.
**Fix:** pass multiple config files in a **single** `--configfile` invocation:
```bash
snakemake --configfile config/test_config.yaml config/gpu_config.yaml   # correct
# NOT: --configfile config/test_config.yaml --configfile config/gpu_config.yaml
```

**Companion form: positional target after `--configfile`.** The same `nargs="+"` argparse semantics also swallow a *target* if it follows `--configfile`. The target gets loaded as a (non-existent) config file and the run fails with `FileNotFoundError` on the target path:
```bash
# FAILS — junctions.tsv is interpreted as a second --configfile argument
snakemake --cores 4 --use-conda --configfile config/test_config.yaml results/.../junctions.tsv

# OK (canonical) — `--` terminates the configfile list; order-independent
snakemake --cores 4 --use-conda --configfile config/test_config.yaml -- results/.../junctions.tsv

# OK — put the target before --configfile (order-dependent; fragile if more flags get added later)
snakemake --cores 4 --use-conda results/.../junctions.tsv --configfile config/test_config.yaml
```
Dry-runs (`-n` flag between `--configfile` and the target) accidentally hide this bug because the flag breaks the nargs sequence. Caught 2026-05-17 while running the chr22 test pipeline for [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224).

### `srcdir()` not available at .smk module level
`srcdir(path)` was a Snakemake-7 helper for resolving paths relative to the calling Snakefile. In Snakemake 8 it is no longer exposed at the `.smk` module scope and a top-level call (`sys.path.insert(0, srcdir("../scripts"))`) fails with `NameError: name 'srcdir' is not defined` at parse time — CI's `snakemake -n` dry-run will catch it.
**Fix:** use `workflow.basedir` instead, which is always exposed and points at the Snakefile's directory (the repo root in this project):
```python
sys.path.insert(0, os.path.join(workflow.basedir, "workflow", "scripts"))
```
Caught in [PR #358](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/358) when a code-review suggestion swapped the working `workflow.basedir` form for the deprecated `srcdir()`.

### `from __future__ import annotations` in `script:`-invoked Python files
Snakemake's `PythonScript.write_script` (`snakemake/script/__init__.py:807`) unconditionally prepends its `snakemake = pickle.loads(...)` preamble before the source. There's a `PY_PREAMBLE_RE` regex on line 44 with a TODO to "use this to find the right place for inserting the preamble", but the regex was never wired up. So **any** workflow script using `from __future__ import annotations` (or any other `__future__` import) at the source level will fail at execute time with `SyntaxError: from __future__ imports must occur at the beginning of the file`. CI's `snakemake -n` dry-run does NOT catch this — the wrapper is only generated at execute time.

**Fix:** drop `from __future__ import annotations` from any file invoked via `script:`. If the file uses PEP 604 syntax (`X | None`) or generic builtins (`list[int]`), rewrite those to `Optional[X]` from `typing` (Python 3.9-compatible). Other workflow scripts may share this latent issue — tracked in [Issue #461](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/461). Caught on [PR #457](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/457) during the chr22 end-to-end run for [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204).

### What `snakemake -n` (dry-run) does NOT catch
`pipeline-snakemake-dry-run` in CI walks the DAG without (a) building conda envs or (b) executing scripts. Bug classes structurally invisible to dry-run:
- **Conda solver failure** — e.g. `IMGTgeneDL>=0.7.0` pinned when PyPI max is `0.6.1`, or upstream channel drift (the bioconda libdeflate/htslib transition that re-broke `star.yaml`, [Issue #629](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/629)). Envs are built lazily on first execute. **Backstopped by the `pipeline-conda-env-solve` CI job ([Issue #646](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/646)), which dry-solves every `workflow/envs/*.yaml` on linux-64; it runs on *every* PR (not path-filtered) so it also catches upstream drift with no env-file change.** Watch the conda prerelease-ordering gotcha when pinning: `2.7.11b` sorts *below* `2.7.11` (a `b` suffix is a beta prerelease).
- **Subprocess CLI typos** — e.g. `stitchr -species HUMAN` (correct flag is `-s` / `--species`). Only the real subprocess argparses the args.
- **NaN flowing through pandas into subprocesses** — curated unit-test fixtures hide single-chain rows that the production VDJdb TSV contains.
- **Snakemake `script:` wrapper incompat** (preceding section).
- **Path discovery silent no-ops** — `if [ -d "$DISCOVERED_PATH" ]; then ...; fi` skips silently when the path is wrong; sentinel still touches.

For new rules, run the chr22 integration locally before merge (see `feedback_integration_run_for_new_rules.md` in role memory). All 5 classes above hit [PR #457](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/457) and were caught only by the first integration run, after the PR had been bot-reviewed + CI-green + reviewer-comments-addressed.

## Reference + Index Layout (nf-core convention)

Post-Issue #63, the pipeline keeps three distinct top-level directories:

| Dir | Contents | Gitignored | Examples |
|-----|----------|------------|----------|
| `references/` | User-provided + downloaded reference data | yes | `GRCh38.primary_assembly.genome.fa`, `gencode.v47.annotation.gtf.gz`, `vdjdb/<release>/`, `imgt_germlines/` |
| `indices/` | Pipeline-built alignment indices | yes | `indices/hisat2/`, `indices/star/` |
| `resources/test/` | Small committed test fixtures + chr22 local-dev cache | partial (test fixtures may be committed; chr22 data is gitignored) | `chr22.fa`, `chr22.gtf.gz`, `hisat2_index/` |
| `models/` | Committed ML model artifacts (fitted, small, reproducible) | no | `models/calibrator_v1.joblib` |

`models/` is the stable production home for committed model artifacts (e.g. the immunogenicity `calibrator_v1.joblib` consumed by the default DAG via `config.yaml` `calibrator.artifact`). It is committed by default (not gitignored), unlike `references/`/`resources/`, so a production artifact and any future re-fit live in the tree rather than a gitignored path. Chosen over `resources/calibrator/` in [Issue #908](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/908) to avoid a `.gitignore` negation and to keep the "`resources/` is retired" convention clean.

The production `resources/` directory is no longer used. `setup_vm.sh` performs a one-shot idempotent migration on existing VMs (moves old `resources/*` → `references/` and `resources/{hisat2,star}_index/` → `indices/{hisat2,star}/`). MHCflurry sentinel moved from `resources/mhcflurry_models.done` to `<mhcflurry.models_cache_dir>/.download_done` (default `~/.mhcflurry/`); the actual model cache still lives in MHCflurry's platformdirs default.

## HISAT2 Index Cache Invalidation

Snakemake skips the index download if `indices/hisat2/` already exists (it checks for an `index.done` sentinel file). Changing `hisat2_prebuilt_url` in `config/config.yaml` does **not** invalidate this cache — the old index silently persists and will be used on the next run.

**When changing `hisat2_prebuilt_url`:** delete the index directory on the VM before running:

```bash
gcloud compute ssh neoepitope-pipeline --zone=europe-west4-a --tunnel-through-iap \
  --command="rm -rf ~/splice-neoepitope-pipeline/indices/hisat2/"
```

(The chromosome naming mismatch in Issue #148 was caused by this exact scenario, when the index lived at `resources/hisat2_index/` pre-#63.)

## Config Migration Notes

### `assembly:` section removed (PR #99)
The old `config.yaml` had an `assembly:` block (`upstream_nt`, `downstream_nt`, `contig_length`). These keys were replaced by `translation.peptide_lengths` in PR #99; flank size is now derived automatically as `3 * (max_length - 1)`. If you have a local override `config.yaml` still containing `assembly:`, Snakemake will silently ignore those keys — remove them to avoid confusion.

## Known Dependency Issues (Fixed)

### `hisat2.yaml` — samtools omitted (libdeflate conflict)
`regtools >= 1.0.0` requires `libdeflate >= 1.26`. No bioconda linux-64 `samtools`/`htslib` build is compiled against `libdeflate >= 1.26`, so conda cannot satisfy both in the same env. (`hisat2` itself does NOT have this constraint — only regtools does.)
**Fix:** `samtools` is omitted from `hisat2.yaml` entirely. The pipeline uses system samtools installed via `apt-get` in `setup_cloud.sh` (currently 1.13 on Ubuntu 22.04). System PATH is visible inside activated conda envs, so no path wiring is needed.
**Workarounds that were tried and rejected:**
- `samtools >= 1.20` pin — made things worse; solver produced an env without samtools (exit 127 on linux-64)
- Fully unpinned samtools — solver falls back to samtools 1.3.1 (2016, 10 versions behind); not acceptable long-term
**TODO(#237):** revisit once bioconda ships an htslib/samtools build against libdeflate >= 1.26. Re-tested 2026-05-06: bioconda's 2026-03 samtools/htslib refactor (samtools 1.23.1) did NOT lift the cap — solver still reports `libdeflate >=1.20,<1.26.0a0` for modern htslib builds, conflicting with regtools 1.0.0's `libdeflate >=1.26`. Workaround remains in place.

### `run_mhcflurry.py` — Class1PresentationPredictor genotype API
`Class1PresentationPredictor.predict()` is a genotype-level call: pass all patient HLA alleles at once (≤6 as a list), get one best-allele prediction per peptide back. Do NOT repeat a single allele N times (that was the `Class1AffinityPredictor` convention and raises `ValueError`).
`predict_to_dataframe()` does not exist on `Class1PresentationPredictor` — use `predict()` which returns a DataFrame directly.

### `python.yaml` — per-platform torch selection (PEP 508 markers)
`python.yaml` selects torch by host platform via PEP 508 markers in the `pip:` block, so the env builds everywhere the project runs:
- **macOS arm64 (local CPU keep-alive, the active target):** `torch >=2.12` — a plain PyPI wheel (CPU/MPS). The `+cu126` build has **no arm64 wheel**, so without the marker the env build fails outright (`No matching distribution found for torch==2.12.0+cu126`).
- **Linux (GPU revival path):** `torch ==2.12.0+cu126` via the cu126 channel (`--extra-index-url https://download.pytorch.org/whl/cu126`), gated to `platform_system == "Linux"`. The `+cu126` local-version tag forces the **Pascal-compatible** build (PyPI's same-numbered wheel lacks the tag and sorts lower per PEP 440); preserved live-in-file for a P100 revival. The extra-index-url is harmless on macOS (unused).

The full cu126/Pascal rationale + the empirical P100 verification ([Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352)) is archived in [`docs/legacy/gcp_p100_setup.md`](docs/legacy/gcp_p100_setup.md). When a second **live** GPU target is funded, promote the markers to explicit per-target env files (`conda: config["conda_envs"][...]`) — markers can't distinguish Pascal-Linux from modern-CUDA-Linux. Verified 2026-06-26: the marker form builds + runs the chr22 CPU core green on macOS arm64 (9/9 steps).

`_has_gpu()` in `run_mhcflurry.py` uses a PyTorch smoke-test kernel (not TensorFlow) to catch the SM-incompat case: `torch.nn.functional.relu(torch.zeros(2, device="cuda"))`. TF reported GPU available even when PyTorch kernels would fail — both must work because MHCflurry 2.2.x uses PyTorch for inference.

### `run_mhcflurry.py` — no ProcessPoolExecutor with GPU
Running MHCflurry alleles in parallel with `ProcessPoolExecutor` crashes on GPU: each worker process initialises its own CUDA context, competing for the same device. Even on CPU, 6 workers × ~8 GB model = ~48 GB RAM → OOM on a 52 GB VM.
**Fix:** single `predict()` call in the main process with all alleles as a genotype. GPU parallelism applies within that call (all peptides batched at once by TF/PyTorch).

## sra-tools Note
Use version `3.1.1` on GCP VMs — newer versions (3.4.x) have a segfault bug.
On macOS arm64, sra-tools conda installation is unreliable due to libcurl/openssl conflicts. Use ENA HTTPS download instead (see `scripts/prepare_test_data.sh`).

## regtools Argument Order
`regtools junctions extract` requires all options before the positional BAM argument. Placing `-o` after the BAM causes `Error parsing inputs!(2)`.
```bash
regtools junctions extract -s XS -a 8 -m 50 -M 500000 -o out.bed input.bam
```

## regtools BED12 — `chromStart`/`chromEnd` are anchor outers, NOT donor/acceptor
`regtools junctions extract` emits BED12. Columns 2-3 (`chromStart`, `chromEnd`) are the **anchor outer boundaries** of the spliced read pile-up, not the intron donor/acceptor. The actual intron coords are recovered from `blockSizes` (col 11) and `blockStarts` (col 12):
```
donor    (0-based)            = chromStart + blockSizes[0]
acceptor (0-based, exclusive) = chromStart + blockStarts[1]
```
Treating cols 2-3 as donor/acceptor shifts every junction by the anchor lengths (typically 100–150 bp on each side) — and silently misses every GENCODE-annotated junction in downstream filtering. Issue #370 replaced the buggy inline `awk` in `alignment.smk` with `workflow/scripts/bed12_to_junctions.py`, which does the blockSizes math correctly.

**Safety net (CI canary, [Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377)):** unit tests of `bed12_to_junctions.py` catch *that* specific off-by-anchor bug but not the *next* one of the same family (a swapped BED12 column, a future off-by-one, a regtools version bump that changes BED12 semantics). The `annotate-flag-canary` CI job (`.github/workflows/annotate-canary.yml`) adds a second, independent source of truth: it runs `regtools junctions annotate` on a hermetic synthetic fixture (`resources/test/annotate_canary/`, regenerable via its `regenerate.py`) and asserts our home-rolled `annotated` flag agrees with regtools' `known_junction` flag for ≥99% of junctions (agreement **and** coverage, so a total coordinate desync can't pass vacuously). A coordinate-semantics divergence fails the job loudly. Comparison logic: `workflow/scripts/crosscheck_annotate_flag.py` (pure functions unit-tested in `workflow/tests/test_crosscheck_annotate_flag.py`; the end-to-end regtools run is the path-filtered CI job, not pytest). The job is path-filtered to PRs touching `alignment.smk`, the junction-extraction/annotation scripts, or the fixture. **regtools gotcha baked into the script: a gzipped GTF makes `regtools junctions annotate` silently emit zero rows — the script gunzips a `.gz` GTF to a temp file first.**

The STAR path is unaffected by the anchor-outer issue — `SJ.out.tab` cols 2-3 are 1-based intron donor/acceptor directly. STAR has its own silent-contamination bug class on **col 4** instead (next section).

## STAR `SJ.out.tab` — `strand=0` rescue from intron motif

`SJ.out.tab` col 4 encodes strand: `0=undefined, 1=+, 2=-`. STAR sets col 4 = 0 when it cannot infer strand from the intron motif (non-canonical splice site, or insufficient evidence in a non-stranded library). The original inline awk in `alignment.smk` silently emitted these records with strand `.`, which flowed through to `assemble_contigs.py` — `bedtools getfasta` without `-s` takes forward-orientation flanking sequence for strand-`.` records, so true minus-strand junctions got sequence-reversed contigs and `translate_peptides.py` read them in the wrong frame.

**Fix:** `workflow/scripts/star_sj_to_junctions.py` rescues strand from **col 5** (intron motif) when col 4 = 0:

| col 5 | motif    | strand |
|-------|----------|--------|
| 0     | non-canonical | dropped (cannot infer) |
| 1     | GT/AG    | `+` |
| 2     | CT/AC    | `-` |
| 3     | GC/AG    | `+` |
| 4     | CT/GC    | `-` |
| 5     | AT/AC    | `+` |
| 6     | GT/AT    | `-` |

Truly non-canonical (col 5 = 0) and unknown codes are **dropped** rather than emitted as strand `.`, on the theory that contaminating the candidate set is worse than slightly under-recalling. Per-run breakdown (direct vs rescued vs dropped) is logged at INFO. Issue #374 replaced the inline awk with this script, mirroring the HISAT2 path's [bed12_to_junctions.py](workflow/scripts/bed12_to_junctions.py) extraction.

The HISAT2 path is unaffected by this — `regtools junctions extract` derives strand from the XS BAM tag set during alignment and never emits `.` strand.

## UCSC vs ENSEMBL Chromosome Naming
Both naming conventions use "GRCh38" in filenames, making it easy to mix them silently.
- `hg38_*` (UCSC) — chromosomes have `chr` prefix: `chr1`, `chr2`, ..., `chrM`
- `grch38_*` (ENSEMBL) — no prefix: `1`, `2`, ..., `MT`

The GENCODE primary assembly FASTA (`GRCh38.primary_assembly.genome.fa`) uses **UCSC naming** (`chr` prefix) despite being distributed by GENCODE/ENSEMBL. All prebuilt HISAT2 indices and BED files in this pipeline must use `hg38_*` (UCSC), not `grch38_*` (ENSEMBL). A mismatch causes `bedtools getfasta` to return empty sequences, silently skipping all junctions in `assemble_contigs.py` (Issue #148).

## Local Test Dataset (chr22)
For development and testing on macOS (M1, 8 GB RAM):
```bash
bash scripts/prepare_test_data.sh   # one-time: downloads reference + FASTQs
snakemake --cores 4 --use-conda --configfile config/test_config.yaml
```
- Reference: chr22 FASTA (UCSC hg38) + GENCODE v47 GTF filtered to chr22
- FASTQs: 500K reads each from a matched gastric cancer pair via ENA HTTPS (no sra-tools):
  - **SRR9143066** — Primary Tumor (gastric cancer surgical section)
  - **SRR9143065** — Solid Tissue Normal (adjacent stomach tissue)
- Both samples are single-end Illumina HiSeq 3000; HISAT2 handles this via `-U` mode
- HISAT2 index stored in `resources/test/hisat2_index/` (separate from production)
- All test outputs go to `results/test/` and `logs/test/`
- **STAR is not usable for local development** — its genome index build requires >8 GB RAM, exceeding the M1 8 GB limit. HISAT2 was chosen for local testing specifically because its index fits within available memory.

## Board status governance — late-commitment Kanban (left side)

The GitHub project board (user project #9, "JH M Lee Lab") runs **late-commitment Kanban**: the commitment point — where an option becomes committed work — sits at the **`Backlog → Ready` boundary**, not at intake triage. The three left-side statuses:

| Status | Meaning | Milestone / Target |
|--------|---------|--------------------|
| **No Status** | Untriaged intake inbox — newly filed, not yet categorized. | none |
| **Backlog** | Triaged but **uncommitted** options: stage S#, priority, size, role label(s), priority rationale set - but not yet committed. | **none** - a Backlog item legitimately has no milestone/Target (not a drift finding) |
| **Ready** | **Committed** + Definition-of-Ready pull queue: refined enough to start, committed on its `arc:` label. | milestone only for milestone-tracked lifecycle work; Target inherited from that milestone, or set by hand if a flow item has a real deadline (else empty) |

**The commitment act (`Backlog → Ready`)** confirms the issue meets the **Definition of Ready** (scope clear, blockers cleared) - that *is* the commitment, signalled by the Status move itself. It does **not** require a milestone. Most board work is **reactive/improvement flow** (PM governance, Dev tooling, MM memory curation, exploratory research) and commits **milestone-free** on its `arc:` label; give a flow item a board Target date by hand only if it has a real deadline, else leave Target empty. Only **milestone-tracked lifecycle deliverables** (the pipeline S1-S7 research track) *also* take an `i<N> - S<N>` milestone here via the capacity decision tree - that `gh issue edit N --milestone "<full name>"` edit drives **Target-date sync** (`Target = milestone.due_on`) + a capacity recheck via the `recheck_dispatch.py` hook. Commitment is **PM-coordinated** (the capacity call needs the cross-cutting portfolio view); DoR-readiness is usually confirmed by the implementing role. (Why milestone-free flow: [Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) facet 2, 2026-06-30 - we are an all-Kanban shop; see "Work types" under the three-axis model below.)

**The Status move is the commitment signal, NOT the milestone.** A freshly-triaged issue - including a new sub-issue - enters Backlog uncommitted and is committed by the move to Ready at its own `Backlog → Ready` boundary; a milestone, when one applies (lifecycle work only), is assigned there too but is no longer the signal. Sub-issues do **not** inherit a parent's milestone (parents are themselves un-milestoned - see "Parent/epic Status" below and `.agents/memory/shared/feedback_parent_sub_issues.md` §Inheritance).

**Phases-as-sub-issues smell (heuristic, not a gate):** at triage / `Backlog → Ready`, an Issue whose internal phases map to distinct PRs / roles / commit-points is under-decomposed — file a sub-issue per phase and convert the parent to a structural epic. Trigger = *independent shippability*, not "has phases". Full rule + [Issue #569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569) precedent: `.agents/memory/shared/feedback_parent_sub_issues.md` → "Phases-as-sub-issues".

Left-side transitions: **intake-triage** (`No Status → Backlog`) → **commitment** (`Backlog → Ready`, PM-coordinated) → **pull** (`Ready → In progress`). The **right side** (`In progress → Ready for review → In review → Done`) is unchanged. Giving `Ready` a distinct job structurally kills the JIT-Ready anti-pattern: you cannot reach In progress without first crossing the commitment act.

Full rules: `.agents/memory/shared/feedback_board_hygiene.md` (sweep cadence, DoR, commitment-act mechanics) and `.agents/memory/feedback_milestones.md` (milestone naming, the capacity decision tree).

### Parent/epic Status — the "Epic" park (Pattern A2)

Parents/epics do **not** flow through the leaf workflow columns. They sit in a dedicated **`Epic`** Status option, and their real progress is read off GitHub's **native sub-issue progress bar** (completion %, auto-derived from children). Decided in [Issue #776](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/776) (2026-06-19) — **Pattern A** (eliminate the drift class, don't police it) + implementation **A2** (park + native bar, not a custom "Epic status" field).

**Why:** a parent given a single draggable Status drifts silently whenever its children move without it — forward-drift (parent ahead of children, e.g. [Issue #680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680)) or stale-terminal (closed/done parent stranded mid-column, e.g. [Issue #232](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/232)). A2 removes the draggable Status from parents entirely, so parent state **cannot disagree** with children. The progress bar is *computed*, never stored — drift-proof by construction — and native (no custom field to maintain). Tradeoff accepted: parents lose To-Do/In-Progress/Done granularity (completion-% only), fine at our low parent count. (Rejected: Pattern B derived-mirror + re-mirror sweep — keeps policing drift and inherits the [Issue #406](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/406) read-back lag; A1 dedicated Epic-status field — a *stored* field that can re-drift unless auto-derived, plus custom-field maintenance.)

**Operating rules:**
- A parent (`subIssuesSummary.total > 0`) belongs in Status **`Epic`** while open; on completion it closes → **Done** (closed parents are terminal, never re-parked).
- Read parent progress from the **native sub-issue bar**, not the Status field. A full bar on an open parent (e.g. all children done) signals ready-to-close.
- The `recheck_parent_status` check (one of the checks dispatched by the `recheck_dispatch.py` hook) **drops its child→parent Status mirror for parents** (they no longer mirror a leaf state); it narrows to leaf-only. **Landed in [Issue #794](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/794):** `classify_drift` suppresses the FORWARD/BACKWARD ladder-mirror for a parent parked in the off-ladder `Epic` Status — without the guard `Epic` reads as rank-0 (Backlog), so any active child spuriously flagged BACKWARD DRIFT, fighting the park. The completion-close signal (all children closed → parent should close to Done) and the [Issue #632](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/632) NOT_PLANNED scope-review flag are **preserved** (neither is a leaf-status mirror).
- Coupling: with parents parked in `Epic` and roadmap visibility carried by the `arc:` label, the milestone-pin-on-parent anchor was no longer load-bearing → **[Issue #690](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/690) sub-question A retired it** (2026-06-19): **parents/epics go un-milestoned.** Milestones are dated stage/time slices on **leaves** (set at `Backlog → Ready`); an epic's cross-iteration roadmap visibility rides its `arc:` label (+ an arc-grouped Project view), never a milestone pin. Matches standard epic practice — epics span sprints and never enter one; their leaves do. Migration: stripped the pins from the 3 then-anchored parents (#416, #547, #680; all arc-covered) + cleared their board Target date. **Sweep rule:** an open parent carrying a milestone is now drift — strip it (after confirming an `arc:` label carries visibility).

Full rule: `.agents/memory/shared/feedback_board_hygiene.md` (parent-status governance section).

## Three-axis work model (stage / arc / due-date)

Work on board #9 is structured along **three orthogonal axes**, each carried by the object that fits it (de-overloading landed by [Issue #693](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/693)):

| Axis | Question | Carried by |
|------|----------|-----------|
| **Stage** | where in the data-science lifecycle? | the **milestone** (`S<N>` in its name) - **lifecycle-track work only**; reactive/improvement flow work carries no stage |
| **Arc** | which long-running narrative throughline? | a never-closing **`arc:<slug>` label** on the issue (+ `arc-phase:active\|next\|later` focus slate, **~3 active, advisory** per [Issue #1102](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1102) — full spec: `docs/superpowers/specs/2026-06-05-arc-work-structuring-design.md`) |
| **Due-date** | by when? | **lifecycle work:** the milestone's `due_on` → synced to the board **Target date** field. **Flow work:** Target date set by hand only if there's a real deadline, else left empty |

**Work types - flow / milestone / sprint (we are an all-Kanban shop, [Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) facet 2, 2026-06-30).** Three kinds of work, committed three ways:
- **Reactive / improvement flow** (the default) - PM governance, Dev tooling + hooks, MM memory curation, exploratory research. **Kanban flow**: commit Backlog → Ready on the `arc:` label, **no milestone**; Target only if there's a real deadline. This is most of the board.
- **Forecastable, batchable, deadline-driven delivery** - committed via a **Scrum sprint** (the GitHub Iteration field). **None exists today**, so the sprint track is a **reserved tool**, stood up only when such a push appears (e.g. a feature-freeze hardening run toward a paper submission). Not instantiated now - no empty track.
- **Lifecycle deliverables** (the pipeline S1-S7 research track) - committed with a dated **`i<N> - S<N>` milestone** (a *goal*, not a time-box; the `i<N>` is variable-length). This is the milestone column above.

Web-grounded canonical direction: **flow the reactive/improvement/tech-debt work, sprint the forecastable feature delivery.** The inverse (sprinting meta-work via short milestone boxes) is a double anti-pattern - it has no real-world precedent and overloads milestones as time-boxes. Role-meta milestones (`pm-i*` / `dev-i*`) were exactly that fiction and are **retired** (see milestone-naming below).

**Rendered active-arc slate ([Issue #759](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/759)).** `scripts/pm/arc_taxonomy.tsv` is the source of truth for the arc **roster + slate phase** (membership is label-authoritative - the `arc:<slug>` label set at triage - per [Issue #973](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/973)) but is invisible day-to-day, so the slate is mirrored, human-readably, into a **pinned `📌 Active arc slate` draft card** at the top of board #9 — every `arc:<slug>` grouped by `active` / `next` / `later` phase, glanceable without leaving the board view. (The board README is *not* the surface — in the new Projects experience the README renders only on the Settings page, so it's repointed at the card.) The card is **hand-maintained**: refresh it at each arc review (PM-coordinated) whenever an arc's phase changes (membership follows the live `arc:<slug>` labels, not a hand-kept list). It is a *surfaced view*, never the source of truth - the TSV is canonical for the roster + phase, the labels for membership; drift between them (a stale `arc-phase`, a multi-arc issue, an unknown slug) is caught by `bash scripts/pm/apply_arc_labels.sh --check` and fixed by the reconcile (no `--check`). (Decided 2026-06-26 over a live renderer / generated docs artifact: arc-phase changes only happen at periodic PM-coordinated arc reviews, so the manual refresh folds into a ritual already touching these labels. The draft card is skipped by `scripts/board_open_items.py` so it never pollutes intake/hygiene sweeps.)

**Milestone naming (terse, post-#693):** `i<N> - S<N> - <Stage Name>` - e.g. `i5 - S3 - Data Preparation`. The trailing ` - <Arc>` content-suffix was **dropped**: a cross-cutting theme belongs on a label, not a milestone - [GitHub allows one milestone per issue](https://docs.github.com/en/issues/using-labels-and-milestones-to-track-work/about-milestones) and milestones *close*, whereas an arc never does (themes/epics → labels is documented GitHub practice). The arc now lives only as the `arc:<slug>` label (taxonomy in `scripts/pm/arc_taxonomy.tsv`). Excluded from the rename: closed milestones (historical record). The `pm-i*` / `dev-i*` **role-meta milestone format is retired** ([Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) facet 2, 2026-06-30) - that work now flows milestone-free (see "Work types" above); no successors are ever opened. The closed `pm-i*` / `dev-i*` milestones stay as historical record.

**Four distinct senses of "iteration" — do not conflate:**
- **`i<N>` (ours)** = a *pass through the DS lifecycle* (S1→S7); **variable-length, work/learning-driven**. Lives in the milestone name. Not a calendar time-box.
- **Scrum sprint** = a *fixed-length calendar time-box*. **We run no sprints today** - we are an all-Kanban shop; the sprint is a **reserved tool** for a future forecastable batch-delivery push only (see "Work types" above).
- **GitHub Iteration field** = GitHub's implementation of Scrum sprints (auto-rolling fixed windows, `@current`). **Held in reserve, not adopted** - see below.
- The lifecycle **stage** `S<N>` is occasionally called an "iteration" loosely — it is a *phase*, not a time-box.

**Time-boxing decision ([Issue #693](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/693), 2026-06-10; reaffirmed [Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) facet 2, 2026-06-30): Target date, not the Iteration field.** The chronological clock lives in milestone `due_on` → the board **Target date** field (auto-synced by `recheck_dispatch.py`), or a hand-set Target on a milestone-free flow item with a real deadline. The native Projects **Iteration field is held in reserve, not adopted** - it encodes fixed Scrum-cadence sprints, which fit only forecastable batch delivery; with none today we run pure **Kanban/flow** (late-commitment Ready queue + WIP limits) *and* variable-length `i<N>` lifecycle passes. The **Roadmap view** plots off **Target date**; **Start date is left unused** (empty on all items). Pointing the Roadmap's date field at Target date is a **manual UI step** - ProjectV2's API can't mutate view config. **Stand the Iteration field up** only if/when a forecastable-delivery push appears (then sprint swimlanes / `@current` auto-roll / velocity Insights become useful) - none are used today. Optional future nicety: populate **Start date** for roadmap *duration bars* (start→deadline) - a deliberate date-field choice, still not the Iteration field.

## Priority semantics - coarse class-of-service, not a fine rank ([Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) facet 3, 2026-06-30)

**Priority (`P0`-`P3`) is a coarse class-of-service signal, deliberately NOT a fine ranking.** The Backlog is ~98% `P2`/`P3` (P2 alone ~64% of 114 items, 2026-06-30), so a band that holds the majority was never a ranking - it is a binary. We make that explicit and stop pretending otherwise:

- **Read the band coarsely.** `P0`/`P1` = expedite / act-soon class; `P2`/`P3` = the standard option pool. P2-vs-P3 is **not** a meaningful sequence - do not agonize over it, and do not read a P2-heavy Backlog as a drift finding (it is the expected shape of an idea-rich, low-stakes-intake research project).
- **The real commitment-time selector is dynamic, not the priority band:** `arc:` focus slate (`active` phase) + freshness + Definition-of-Ready, decided when capacity signals a pull - not pre-baked into a stored rank. This is already how `shared/feedback_best_next_issue.md` selects; this just makes it the *documented* selector and demotes priority to a coarse class hint feeding into it.
- **Prune, don't rank.** An unordered option pool must stay bounded, so the refinement pass (widening [Issue #814](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/814)) actively **declines/closes** stale low-stakes options (closing comment + reason) rather than letting them accrue forever. Pruning is the hygiene that makes an unordered backlog workable.
- The required `**Priority rationale:**` line **stays** - it keeps us honest about which coarse class an item is in (why P0/P1 vs the pool). We simply no longer over-read the band as a rank.

**Why** (web-grounded 2026-06-30): this is the canonical Kanban answer. David J. Anderson - **["Banish Priority and Prioritization"](https://djaa.com/ban-priority-and-prioritization/)** - argues priority is a *proxy variable* masking real risk (cost of delay, skills, dependencies), and that a backlog clustering at one band is the documented failure mode; the prescription is **risk profiles / classes of service** in place of priority, **selection / scheduling / replenishment** (dynamic pull) in place of prioritization, and to **"maintain backlogs as unordered lists."** **What we rejected:** *forced-ranking a per-role top-N* (the explicit anti-pattern - a manual total-order on a deep backlog "requires 100% manual intervention" for no flow benefit at our scale), and a *separate Class-of-Service field* (our coarse `P0`-`P3`, read as classes, already carries the same signal at our scale - revisit only if classes ever need distinct *policies*, e.g. a fixed-date SLA class).

## WIP limits — starvation, swarm, and review-debt guards

The board has three WIP guards running in `check_ready_queue.sh`, each addressing a different constraint:

| Guard | Section | What it prevents | Default | Constraint |
|---|---|---|---|---|
| **Ready starvation** | § Replenishment trigger | Empty Ready queue (no curated shortlist when a role sits down) | floor 5/role, cap 23 total | Agent throughput |
| **In-progress swarm** | § Per-role In-progress cap | Duplicate-PR swarms (two roles accidentally racing the same work) | cap 3/role | Agent throughput |
| **Review debt** | § Review-column WIP limit | Unreviewed PRs piling up past the human's review capacity | limit 10 | Human review bandwidth |

The first two (Ready starvation + In-progress swarm) constrain agent throughput. The third (review debt) constrains the real binding constraint for a 1-human/N-persona team: **human review bandwidth** ([Issue #928](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/928), grounded in [`research/agent_team_governance_research_2026-07.md`](research/agent_team_governance_research_2026-07.md) §7). All three are **advisory** (house style: guards stay advisory until a defect recurs >=2x).

**These three are the *only* WIP limits, and they all count work in flight.** The `arc-phase:active` focus slate is **not** a fourth one: it groups *categories*, not items, and a WIP limit on a category has no external analogue (SAFe's portfolio WIP limit counts epics moving through Kanban states, never strategic themes - [agility-at-scale](https://agility-at-scale.com/safe/lpm/work-in-progress/)). Its former hard cap of 3 was **softened to an advisory `~3` guideline** in [Issue #1102](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1102), after measurement showed `active` does not predict throughput: over 2026-06-19..07-10 a `later` arc out-delivered two of the three `active` arcs, roughly a quarter of closed work carried no `arc:` label at all, and the board's only In-progress item sat in a `next` arc. The cap's sole observed effect was to force a demotion of the highest-throughput arc during the [Discussion #949](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/discussions/949) arc call, for zero behavior change. `apply_arc_labels.sh` now prints a `NOTE` above the guideline and never fails on it (tunable via `ACTIVE_SLATE_GUIDELINE`).

**Replenishment trigger ([Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) facet 1, 2026-06-30) - a proactive shelf with interpretable shortfall.** The per-role Ready floor (default 5) is a **target shelf depth**: keep ~5 DoR-ready items committed per role so that whenever a role sits down there is always a curated shortlist to pull, without paying the commitment-decision tax mid-session. It is **proactive** (kept stocked ahead of demand), **not** gated on current consumption - the buffer's whole value is being stocked *before* demand arrives. A shortfall (`ready < floor`) is routed by whether committable work exists: **`[REPLENISH role]`** when the role has triaged Backlog candidates (commit the highest-priority DoR-ready ones; if none meet DoR, that's a grooming gap - refine/file, **never stuff** low-value work to hit the number), vs **`[GROOMING-GAP role]`** when the role has no Backlog candidates at all (the remedy is intake/grooming, not commitment). In progress is reported as demand *context*, never a gate. The total cap (default 23 = floor x 4 roles + 3 headroom) is a soft WIP limit so Ready doesn't over-deepen. **All four roles - PM / Sci / Dev / MM - are subject to the floor** ([Issue #1006](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1006), 2026-07-04, retiring the [#705](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/705) MM exemption): the exemption's premise (a fixed floor would nag a bursty/blocked MM lane into stuffing junk) was obsoleted by facet-1's `[GROOMING-GAP]` routing - a thin MM lane reads GROOMING-GAP (groom, don't stuff), so a genuinely-empty MM lane is honest, not a violation. Web-canonical calibration ([Businessmap](https://businessmap.io/kanban-resources/getting-started/what-is-wip), [Kanban WIP guides](https://teachingagile.com/kanban/introduction/kanban-wip-limits)): start every workstream at a **uniform** limit and only scale by capacity *after* flow data shows a mismatch (the weekly Service Delivery Review flow-metrics plan, #902 follow-up, is where that tuning lives — this doc names the starting guard; the SDR owns the data-driven retune).

**Rule (In-progress swarm guard): per-role advisory WIP cap of 3 on `In progress`** (PM / Sci / Dev / **MM** = Memory Manager each). MM was previously excluded here "because it implements no board-tracked work, same reason it's out of the #754 Ready count" - but [Issue #1006](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1006) (2026-07-04) retired that premise (MM commits `role:memory_manager` board work and is now a full floor role), so the exclusion's justification collapsed and MM is folded in here for symmetry. More than 3 In-progress items carrying a single role's `role:` label = a flag, **not** a hard block — surfaced in the **Daily Stand-up WIP-awareness beat** (this just gives that existing beat a concrete trigger number). Classify by the item's `role:` *implementer* label, not by who filed it (same convention as the Ready floor).

**Why In-progress cap 3** (matches standard Kanban: "2–3 items per person, per-person WIP is the right starting point"):
- **per-role, not per-column** — roles are the throughput unit here; consistent with the #754 per-role Ready gate and the arc `active`-slate guideline (~3, advisory).
- **cap 3** — top of the standard 2–3/person band; one actively-worked + slack for items blocked pre-PR (work moves *out* of In-progress into the separate `Ready for review` / `In review` columns, so In-progress holds only pre-PR active work). **Tunable** — tighten toward 2 if In-progress fragmentation shows up; the weekly Service Delivery Review (#902 follow-up) is the forum for data-driven retune.
- **advisory, not blocking** — the major tools (Azure Boards, Jira) implement WIP as a *visual warning when exceeded*, not a pull-stop; advisory matches our house style of keeping guards advisory until a defect recurs >=2x (then escalate per the mechanism-over-memory ladder).

**Rule (review-column WIP limit): advisory cap of 10 on `Ready for review` + `In review` combined ([Issue #928](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/928), 2026-07-07).** For a 1-human/N-persona team, human review bandwidth is the binding throughput constraint — "the more you parallelize your code generation, the more review debt you create" ([`research/agent_team_governance_research_2026-07.md`](research/agent_team_governance_research_2026-07.md) §7, hand-verified verbatim from MindStudio + InfoWorld). Default 10 (top of the 5–10 cards/reviewer band from practitioner guidance). Counts PRs only (a PR and its linked Issue both sit in review columns, so counting all cards would double-count). Tunable via `--review-limit` or `REVIEW_WIP_LIMIT` env. The `[REVIEW-DEBT]` surfacing fires in `check_ready_queue.sh`'s status line; it's a Daily Stand-up WIP-awareness signal (same tier as the In-progress cap), not a dispatch gate.

## Coordinator autonomy envelope - reversibility ladder gating dispatch

The pre-authorized policy for what a coordinating agent may commit/dispatch **without** per-action human approval, and where it must stop. It makes the two existing autonomy rules (`autonomy-merge-gate-cadence` + `decision-ratification != action-authorization`) explicit and act-indexed, and is the keystone gating [epic #1072](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1072)'s dispatch rungs 3-5. Full design + rationale: [`docs/superpowers/specs/2026-07-08-coordinator-autonomy-envelope-design.md`](docs/superpowers/specs/2026-07-08-coordinator-autonomy-envelope-design.md); operational memory: `shared/feedback_coordinator_autonomy_envelope.md`.

**The boundary (unambiguous):** reversible coordination is pre-authorized; the irreversible outward act and the novel/governance decision are the human's.

**The 3-tier ladder** (keyed on reversibility / blast-radius; linear, not a matrix, for action-time read speed):
- **Tier A - autonomous** - reversible, in-envelope. Do it.
- **Tier B - autonomous up to the gate** - irreversible-routine: run the whole chain autonomously, the human performs the single irreversible act.
- **Tier C - always stop** - novel / precedent / governance / cross-role-conflict / irreversible-outward beyond merge. Propose first.

**Act -> tier mapping:**

| Act | Tier | Note |
|---|---|---|
| **Commit** (Backlog -> Ready) | A | up to the per-role floor, on the active arc, DoR met |
| **Dispatch** (route ping, trigger pull) | B now -> A on maturation | promotes once the review-axis rung ships + `hook_fires.jsonl` shows safe behavior; the one act whose tier climbs (the Supervised -> Delegated trajectory in miniature) |
| **Merge** (`gh pr merge`) | B | run the chain, stop at the merge command; flips to C when the PR itself trips a Tier-C trigger |
| **Escalate** (flag / ask) | A | never need approval to raise a concern |

**Two refinements:** (1) **Tier C is an override condition, not a fourth act** - its triggers (novel scope / precedent-or-governance change / cross-role conflict / irreversible-outward beyond merge) bump *any* act to always-stop (editing memory or CLAUDE.md is itself Tier C). (2) **Reversible-but-novel splits**: a novel *instance* is Tier A + surface (do it, explain); a *precedent-setting* act is Tier C (stop).

**Altitude:** these are **act tiers** (per-act, reversibility). Distinct from the paper's **working modes** (posture over a workstream, maturing over time = the epic's altitude). Same word root, different altitude.

## GitHub Safety Wrappers

Mechanisms that fire automatically to enforce GitHub-related discipline rules that have broken repeatedly despite being documented in memory. Per the mechanism-over-memory ladder (memory → inline Always-in-effect → mechanism), these are the rung-3 escalation when a rule has slipped ≥2× on the same shape.

> **⚠️ Hook loading — `.agents/settings.json` edits are hot-reloaded (as of Claude Code 2.1.199).** Editing `.agents/settings.json` (the canonical file; Claude Code discovers it through the `.claude` -> `.agents` symlink, so the legacy `.claude/settings.json` path resolves to the same file) mid-session is **re-read live**: the config is reloaded on each subsequent tool call, so a newly-wired hook fires immediately and a removed/changed hook takes effect at once - **no session restart needed**. Verified 2026-07-03 (Claude Code 2.1.199, [Issue #935](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/935)) with a controlled 3-probe test in a fresh no-branch-switch session: the em-dash guard denied a probe *before* a no-op settings edit, *still* denied after it, and *stopped* denying once its hook was deleted from the config (the deleted-hook file was valid JSON, ruling out a parse-failure fallback rather than a genuine live re-read). **The behavior has flip-flopped across builds - which is exactly why the version pin matters.** Three data points: (1) an earlier note ([`research/lab_notebook/pm.md`](research/lab_notebook/pm.md) under `## 2026-05-18`) recorded a mid-session *hot-reload* of a PostToolUse hook (in `settings.local.json`); (2) then on the Claude Code build current 2026-05-29 a mid-session edit silently *deactivated* ALL hooks (every PreToolUse + PostToolUse) until a session restart - the footgun that silently caused [PR #560](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/560) + [PR #562](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/562) to miss auto-boarding; (3) now 2026-07-03 it hot-reloads again. (The 2026-05-18 point is a different file and 6 weeks stale, so corroborating not conclusive - but it independently supports the PostToolUse half below.) Because the behavior is demonstrably **Claude-Code-version-dependent**, on an unknown or older build still verify what's live with `/hooks` after editing. The 3-probe test exercised a **PreToolUse** guard; a per-tool-call re-read would also cover **PostToolUse** hooks (independently corroborated by the 2026-05-18 observation), but **Stop** hooks fire at turn end, *outside* the per-tool-call cycle, and were not tested - don't assume they hot-reload. Also: hooks fire only on **Claude's tool calls**, not on commands a human runs in a terminal - so a firing test must be driven by a Claude session.

> **⚠️ Newline blindness is in the HOOK's own matcher, NOT in the `if:` gate ([Issue #1142](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1142), 2026-07-13).** **A newline is a command separator in shell but NOT in `shlex`.** So a hook that tokenizes the raw command is **silently wrong on every multi-line command** (a heredoc `--body-file`, or simply a second line). It bites in two shapes, and the second is the one that hides:
>
> - **command-start walk** (`at_command_start`): the trigger lands after an ordinary word instead of a separator, so the match never fires and the hook **allows**. Killed `post_gh_pr_create` on the dominant PR shape ([Issue #1130](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1130)), plus `check_board_query_pagination` and `check_gh_issue_develop_parent`.
> - **subcommand split** (`split_subcommands`): two commands on two lines **MERGE into one token block**, and the merge defeats the guard *without* any matcher failing. `git commit --amend` ⏎ `git push --force` resolved to `commit`, so the force push was never examined; `cd /tmp` ⏎ `cd /etc/secret` returned only the first target; `git commit` ⏎ `git push` collapsed to one index so the chain-guard passed. It also **false-positives**: `git push` ⏎ `git branch -f` merged into something that read as a forced push.
>
> **The fix is `_shell_parse.normalize_command(cmd)` before tokenizing** (strips heredoc bodies, turns *unquoted* newlines into separators). It **strengthens** the anti-false-positive property rather than weakening it: a heredoc body is *removed* instead of tokenized, so a comment body that merely *discusses* `gh pr create` still cannot match. **All five `gh`/`git`-matching hooks now call it, and any new one MUST** - shipping without it is shipping a hook that is born dead on the shape we actually use.
>
> **Probe with a second REAL command, not an `echo` prefix.** The merge shape only bites when both lines are commands the hook cares about, so a probe like `echo x` ⏎ `git push --force` passes while the real `git commit` ⏎ `git push --force` fails. The first #1142 sweep used the `echo` prefix, declared four hooks clean, and missed three live ones; the PR's bot review caught them.
>
> **The `if:` gate is NOT the culprit, and de-gating is NOT the fix.** #1142 was originally filed blaming the `if:` filter, on a probe that could not separate the two layers (it used a hook whose *inner* matcher was blind too, so the gate looked guilty). An isolating matched-pair probe (a throwaway hook with **no** inner matcher, substring-matching only) showed `if: "Bash(gh *)"` **fires correctly** on both a heredoc-led and a plain multi-line command: the gate is subcommand-aware. Dropping the gates would add ~180 ms x 9 hooks to *every* Bash call and fix nothing. **Keep the gates; normalize inside the hook.** The gate remains vendor-documented **best-effort**, so it must never be a hook's *only* matcher - every hook self-matches internally and fails open, which is what makes the gate a safe performance optimization rather than a correctness boundary.
>
> **Method note (the reason this took two attempts):** the first fix repaired the matcher one layer *below* the one under suspicion and shipped believing the mechanism was alive. When a hook does not fire, **isolate the layer** - drive the real trigger, read `.agents/hook_fires.jsonl`, and control with a probe that varies exactly one thing. A payload piped into `main()` bypasses both the gate and the harness, and can only ever confirm.

### `@claude` mention guard (user-level PreToolUse hook in `~/.claude/settings.json`)

Refuses any `gh (pr|issue) (comment|create|edit)` whose `--body` contains a literal `@claude` substring, except the exact canonical review-trigger `--body "@claude review"` (and the single-quoted variant). The hook runs `check_at_claude.py` on stdin-piped PreToolUse JSON and emits a `permissionDecision: deny` when the guard fires.

**Registration - user-level, single canonical source ([Issue #799](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/799) / [#665](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/665) gap 1, 2026-07-04):** the guard is registered **only** in `~/.claude/settings.json` (matcher `Bash`, `if: "Bash(gh *)"`, timeout 5), **not** in the project `.agents/settings.json` (the project-level entry was retired in #799). It must fire regardless of cwd: a project-repo comment can be drafted from the project clones *or* the personas-repo cwd (where the Memory Manager operates), and only a user-level (global) registration fires cross-cwd - a project-scoped hook fires only when its own repo is cwd. Keeping both registrations meant a double-fire in the project cwd, and that redundant pair was the exact drift risk #799 closed. The script stays single-sourced in the project repo (`.claude/hooks/check_at_claude.py`, via the `.claude` -> `.agents` symlink); the user-level command references it by **absolute path** (unavoidable - `${CLAUDE_PROJECT_DIR}` is undefined at user scope). Tradeoff: the user-level registration is **machine-local, not version-controlled**, so a fresh machine is unguarded until `~/.claude/settings.json` is set up (acceptable - the guard only matters where an agent runs `gh` locally; CI never drafts these). On a **missing script path the user-level command currently fails open** (a moved/renamed clone means the hook simply never runs); making it fail **closed** (loud breakage, not a silent bypass - the failure mode that bit three times) is the decided fix (Developer, #799), **pending** the Memory Manager wiring it into the user-level entry. That fail-closed hardening plus VCS-capture of the canonical guard block on the personas side are tracked in [#978](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/978); until they land, the user-level registration is the sole live source but its missing-path behavior is still fail-open.

**Why:** the `Claude Code` GitHub Action subscribes to `issues` events and triggers on ANY literal `@claude` — including inside parens, ACs, code spans, quoted historical refs, or markdown link descriptions. The rule lives in `.agents/memory/shared/feedback_no_at_claude_mention.md` and is inlined into shared Always-in-effect, yet broke on PR #359 (2026-05-13) and originally on Issue #272 (2026-05-06).

**Workaround for non-trigger references:** use `@-claude` (zero-width hyphen between `@` and `claude`). The literal substring `@claude` must not appear.

### `gh issue develop` parent guard ([PreToolUse hook](.agents/settings.json))

Refuses any `gh issue develop <N>` (number or issue-URL form) when Issue `N` is a parent/epic — i.e. its `subIssuesSummary.total > 0`. The hook runs `.agents/hooks/check_gh_issue_develop_parent.py` on stdin-piped PreToolUse JSON and emits a `permissionDecision: deny`. It **fails open** on every uncertain path (unparseable command, no issue number, repo unresolvable, `gh`/network error) — only a *confirmed* parent denies — so a hiccup never blocks legitimate `gh issue develop`. Each deny appends one line to `.agents/hook_fires.jsonl` (gitignored, [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453) fire-log infra).

**Why:** branching off a parent creates a PR↔parent `closingIssuesReferences` edge that auto-closes the epic on merge, silently orphaning its sub-issues. That bit [PR #543](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/543) → [parent Issue #538](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/538) (2026-05-28). The "parents have no branches/PRs" rule was Reference-tier in `memory/shared/feedback_parent_sub_issues.md` and didn't load at the gh-develop target-pick moment — rung-3 mechanism escalation ([Issue #549](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/549)).

**Workaround:** branch off a leaf sub-issue instead, or file a closure sub-issue (cf. [Issue #548](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/548)) under the epic and develop off that.

### Board-query pagination guard ([PreToolUse hook](.agents/settings.json))

Refuses any command-start `gh api …` invocation whose args contain `projectV2` **and** an `items(first: …)` connection but **no** pagination token (`hasNextPage` / `endCursor` / `after:`). The hook runs `.agents/hooks/check_board_query_pagination.py` on stdin-piped PreToolUse JSON (pure string inspection — no `gh`/network I/O) and emits a `permissionDecision: deny`. It **fails open** on any parse miss or untokenizable command. The two safe patterns pass untouched: the per-issue `issue(number:N){ projectItems … }` lookup (no `projectV2.items`) and `scripts/board_open_items.py` (carries `pageInfo`/`after`; also a python invocation, never inspected). A `gh pr comment`/`gh issue create` whose *body* merely discusses the pattern does not match — only a `gh api` at a command start counts (shlex command-start tokenization). Each deny appends one line to `.agents/hook_fires.jsonl` (gitignored, [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453) fire-log infra).

**Why:** the board (#9) holds 700+ items sorted **Done-first**, so a single-page `projectV2 { items(first: N) }` query silently returns mostly-closed items and *hides the open Ready/In-progress work past position N* — with no error. That truncation produced a wrong "Ready is empty" read on 2026-06-12 (a `first: 100` query showed only the first 14% of 711 items). The correct paginating helper (`scripts/board_open_items.py`) and the rule to use it live in `memory/shared/feedback_board_queries.md`, but the rule keeps slipping at the moment someone types a quick ad-hoc query — rung-3 mechanism escalation ([Issue #717](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/717)).

**Workaround:** use `scripts/board_open_items.py` (`--role`/`--status`/`--json`), or for a single-issue lookup query `issue(number:N){ projectItems … }`. To roll your own board scan, add a `pageInfo { hasNextPage endCursor }` + `after:` cursor loop. (A hook only fires on agent tool-calls, not a human's terminal query; the softer `gh issue list --limit N` truncation is handled by convention — default `--limit 1000`.)

### No-em-dash guard ([PreToolUse hook](.agents/settings.json))

Refuses any `Edit` / `MultiEdit` / `Write` that **net-adds** an em-dash (`U+2014`) or en-dash (`U+2013`) to Claude's newly-written content. The hook runs `.agents/hooks/check_no_emdash.py` on stdin-piped PreToolUse JSON and emits a `permissionDecision: deny` pointing at the plain-hyphen rule. **Character scope = em + en** (both have slipped; horizontal ellipsis and other smart punctuation are deliberately out of scope). It scans only the **delta**: for `Edit`/`MultiEdit` each `new_string` is compared against its `old_string`; for `Write` the new `content` is compared against the file's current on-disk content (empty for a brand-new file). So editing a file that legitimately already contains em-dashes (most historical memory files, including this one) never false-positives unless the edit introduces a *new* one. (`Write` is coarser than `Edit`/`MultiEdit`: it compares whole-file counts, so a single `Write` that both removes a pre-existing dash *and* adds a new one nets to zero and slips through - `Edit`/`MultiEdit` are span-local and exact. `NotebookEdit` is deliberately not guarded - notebook cells are out of scope for now.) It **fails open** on any parse miss, unreadable file, or unexpected shape, and appends one line per deny to `.agents/hook_fires.jsonl` (gitignored, [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453) fire-log infra).

**Why:** Jin-Ho's global "never an em-dash, use a plain hyphen" rule lives in `~/.claude/CLAUDE.md` (memory tier) but slipped in Claude's own additions twice inside two sessions ([PR #916](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/916) em+en, [PR #917](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/917) em), crossing the mechanism-over-memory threshold (>= 2x on the same shape). Rung-3 escalation ([Issue #920](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/920)), sibling to the `@claude` and `gh issue develop` parent guards.

**Escape hatch:** for a genuine verbatim need (quoting external text, a fixture that must contain the character), set `CLAUDE_ALLOW_EMDASH=1` in the environment, or write to an allowlisted path (substring match on `check_no_emdash` and `.agents/hook_fires.jsonl`). The guard only sees Claude's tool-calls, never a human editing in an external editor. **Wiring this hook edited `.agents/settings.json`; as of Claude Code 2.1.199 such an edit is hot-reloaded, so wiring a hook now takes effect live** (see the Hook loading note at the top of this section - mid-session settings edits are re-read live, no restart needed; on the older Claude Code build current 2026-05-29 a settings edit instead deactivated every hook until restart).

### Review-request board-advance (PostToolUse hook, `.agents/hooks/post_gh_pr_review_request.py`)

Not a guard (it never denies) but a **board-automation** sibling of `post_gh_pr_create.py`. When a bot review is requested on a PR - `gh pr comment <ref> --body "@-claude review"` (the canonical trigger, hyphenated here only to dodge the mention guard) or a direct `gh pr review <ref>` - it resolves the PR's linked Issue(s) via `closingIssuesReferences` and sets each card's Status on project #9 to **In review** (`singleSelectOptionId df73e18b`), closing the manual `Ready for review` -> `In review` hop that kept being forgotten and stranding cards (e.g. [#406](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/406), [#234](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/234)). It reuses the sibling's `PROJECT_ID` / `STATUS_FIELD_ID` plumbing (no duplicated field-ID drift), is **idempotent** (already-In-review is a no-op), does **not** overwrite a parent parked in `Epic` or a terminal `Done` card, and **fails open** on any parse miss / untracked repo / `gh` error. A real flip logs one line to `.agents/hook_fires.jsonl` (gitignored, [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453) fire-log infra).

**Why:** the `Ready for review` -> `In review` transition is board discipline that memory did not fix - it slipped repeatedly across roles, meeting the mechanism-over-memory threshold (>= 2x on the same shape). Rung-3 escalation ([Issue #996](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/996)).

**Residual (documented, not solved):** a human reviewing directly on github.com emits no local command, so a local PostToolUse hook cannot see it - that path still needs a manual flip (or a future GitHub Action / webhook). The bot-review path is the dominant one here, so the hook covers the common case.

### Auto-requested first-pass bot review (PostToolUse hook, `.agents/hooks/post_gh_pr_create.py`)

Since [Issue #1073](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1073) the PR-open hook **auto-posts the bot-review trigger** on a non-trivial PR, so the first-pass review no longer waits on a human remembering to offer it. **Non-trivial = non-draft AND no `skip-bot-review` HTML-comment marker in the PR body** (drafts are excluded for free - a draft is not up for review, which is already the hook's own Status semantic; the marker mirrors the existing `<!-- skip-lab-notebook: routine -->` convention). The human `gh pr merge` gate is **unchanged** - the bot is first-pass, never a merge authority.

**The opt-out marker must be on its own line.** This diverges deliberately from the unanchored `_SKIP_LAB_NOTEBOOK` regex, and the reason is empirical: the first unanchored draft matched a *backtick-quoted, mid-sentence mention*, so [PR #1124](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1124) - whose body documents the opt-out - silently opted **itself** out of the review it exists to request. Anchoring to a whole line is what separates *using* the directive from *talking about* it. (The lab-notebook marker carries the same latent hazard, unfixed - a PR body quoting it inline would skip that gate too.)

**Why:** human review bandwidth is the binding throughput constraint for a 1-human/N-persona team ([`research/agent_team_governance_research_2026-07.md`](research/agent_team_governance_research_2026-07.md) §7), so making the bot's first pass automatic widens the gate that caps everything else. It is the highest-leverage rung of [epic #1072](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1072), and the rung the **Dispatch** Tier-B -> Tier-A promotion in the coordinator autonomy envelope is gated on. Chosen over folding the request into the PR-open *checklist* because that checklist is a memory rule (ladder rung 2) and this exact rule already slipped twice on it ([PR #441](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/441) + [PR #442](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/442), one morning) - which is why `bot_review_offer.py` exists at all. Answering a memory slip with more memory is the wrong rung.

**Two interactions, both load-bearing:**
- **No double-request against the merge gate.** `tools/ci/bot_review_offer.py` is a *detector*, not a prompter: `audit_and_merge.sh` prompts only when the literal review trigger is absent from the PR's comments. Once the hook posts that string the gate reads `OFFERED` and skips its prompt, so the two mechanisms compose **with no code change to the gate**. A cross-module test pins the agreement (swap the hook's trigger for the non-firing `@-claude review` reference form and it goes red).
- **Status is delegated, not duplicated.** When the hook auto-requests, it calls `post_gh_pr_review_request.apply_review_request()` (extracted for exactly this second caller) to set **In review** on the PR card + every linked Issue, rather than its usual `Ready for review`. The review-request hook's PostToolUse matcher only fires on *Claude's Bash calls*, so it structurally cannot see this hook's **subprocess** `gh pr comment` - leaving the card at `Ready for review` for the duration of a live review would re-create the exact stranding [Issue #996](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/996) fixed.

**Residual (documented, not solved):** same as its sibling - a hook fires only on Claude's tool-calls, so a PR opened by a human in a terminal gets no auto-request and still falls back to the merge-time offer prompt. That path degrades to *today's* behavior, never worse. Agent-authored PRs are the dominant path here.

### Closure-ritual gate (`scripts/audit_and_merge.sh`)

Refuses to run `gh pr merge` if any `- [ ]` remains on the PR body Test plan OR on any Issue in `closingIssuesReferences` (under that Issue's Acceptance criteria), OR if any linked Issue lacks a `**Priority rationale:**` line, OR — added via [Issue #559](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/559) — if the **would-be squash body** (PR title + body + every commit message) contains a closing keyword (`close|fix|resolve` + `#N`) targeting an Issue **outside** `closingIssuesReferences`, OR — added via [Issue #409](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/409) — if a **role-tagged linked Issue** lacks a `## <merge-date>` lab-notebook entry in `research/lab_notebook/<role>.md` referencing the PR/Issue (routine-ship bypass: `<!-- skip-lab-notebook: routine -->` in the PR body), OR — added via [Issue #665](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/665) — if a **cross-repo closing forward-link** in the PR body (`owner/repo#N` or a full issue URL — the personas-PR→project-Issue case that native `closingIssuesReferences` can't express) points at an Issue whose **Acceptance criteria** still have unticked boxes. Operational usage lives under [Merge workflow](#merge-workflow); this section captures the mechanism-over-memory rationale.

**Why (gates 1+2):** the closure-ritual rule has broken 4× in 10 days despite being inlined into shared MEMORY.md — [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280), [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299), [PR #328](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/328), [Issue #347](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/347). Memory of a declarative rule cannot reliably survive the action-distance from session start to `gh pr merge`. The script enforces at the exact moment of action — the merge invocation itself — so the gate cannot be forgotten regardless of how much downstream context has accumulated.

**Why (stray-closer gate, `tools/ci/stray_closers.py`):** a closing keyword + `#N` in a *commit-message body* auto-closes Issue N on merge (the squash commit inherits it), but the PR's `closingIssuesReferences` API surfaces only PR-**body** link edges — so a pre-merge `closingIssuesReferences` check passes clean and misses it. [PR #543](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/543) auto-closed parent epic Issue #538 this way (2026-05-28; commit body `would auto-close #538`). The gate assembles the full squash text and flags any closer pointing outside the intended set. It **fails open** (non-blocking) on a missing interpreter or `gh` error.

**Why (bot-review-offer gate, `tools/ci/bot_review_offer.py`):** the rule to proactively offer a `@-claude review` after a non-trivial PR (`shared/feedback_github_workflow.md`) broke twice on one morning — [PR #441](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/441) + [PR #442](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/442) merged without it (2026-05-21) — crossing the memory→mechanism threshold ([Issue #443](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/443)). After the closure-ritual checks pass, the gate detects whether the real trigger (`@claude review`, the literal string the GitHub Action fires on — *not* the hyphenated reference form) already appears in the PR's comments. If not: **interactive** runs prompt (a) offer now / (b) skip-trivial / (c) cancel; **non-interactive** runs (Claude/CI — the actor whose slips this targets) **block with exit 1** rather than silently merging, with guidance to offer-and-re-run or pass `--skip-review-offer`. Detection **fails open** to "offered" on a missing interpreter or `gh` error.

**Why (lab-notebook gate, `tools/ci/lab_notebook_gate.py`):** the lab-notebook entry was the one closure-ritual element checked *only post-merge* — `tools/ci/closure_audit.py` (run by the closure-audit workflow) posts a marker comment naming a missing `## <date>` entry *after* the PR has already merged, a cleanup loop, not prevention ([PR #403](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/403) merged 2026-05-19 without a scientist entry; cleaned up in [PR #408](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/408)). [Issue #409](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/409) moves the same check to merge time. The gate **single-sources** `closure_audit.audit_pr_pre_merge` (reusing `is_exempt` / `skip_lab_notebook` / `resolve_roles` / `collect_notebook_gaps`), reads the entry from the **working tree** (so run `audit_and_merge.sh` from the PR branch where the entry was written), and **fails open** (exit 0 + warning) on a missing interpreter or `gh` error. It enforces no *new* policy — it blocks exactly what the post-hoc bot would have flagged, just earlier. This supersedes the "NOT gated by `scripts/audit_and_merge.sh`" note in `shared/feedback_lab_notebook.md`: the routine-vs-non-routine brittleness that previously blocked a deterministic gate is handled by opt-out (the `<!-- skip-lab-notebook: routine -->` marker), not auto-classification.

**Why (stray-AC-box lint, `tools/ci/ac_section_lint.py`):** gate 2 blocks only on unticked boxes *under* a `## Acceptance criteria` heading, so an Issue whose gating boxes live under a different heading has 0 boxes "under Acceptance criteria" and **passes silently** — the AC enforcement is bypassed by a heading-name deviation. [PR #724](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/724) merged [Issue #569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569) with its gating boxes under `## Plan (phased)` this way ([Issue #730](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/730)). The fix is a **convention, not a parser change**: broadening the gate to also read `Plan`/`Tasks`/`Checklist` headings would false-block the *non-gating* boxes those sections routinely carry (exploratory steps, nice-to-haves), making the gate non-deterministic. So the gate stays keyed to **one canonical `## Acceptance criteria` heading**, and a **non-blocking lint** warns when a linked Issue has unticked boxes but no AC section, naming the count + non-AC heading(s). It **never blocks** (the right fix is moving the boxes, not refusing the merge) and **fails open** on a missing interpreter or `gh` error. The pure scan (`closure_audit.scan_ac_boxes`) is the shared helper that the post-merge bot's `check_ac` is being scoped onto in the sequenced follow-up [Issue #726](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/726).

**Why (cross-repo AC gate, `tools/ci/cross_repo_ac_gate.py`):** gate 2 audits ACs only for Issues in the PR's native `closingIssuesReferences`, which GitHub populates **only within one repo**. So a **cross-repo close** — the Memory Manager's personas-repo PR closing a project-repo Issue — slips the AC gate **entirely** (the project Issue's ACs are never checked before the personas PR merges); the cross-repo "close" is expressed only as text in the PR body ([Issue #665](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/665), gap 2). The gate (2b) parses the PR body for cross-repo closing forward-links — a closing keyword (`close`/`fix`/`resolve`) co-occurring on a line with an `owner/repo#N` ref or a full issue URL, in the GitHub form *or* the project's `[Issue #N](url) (closes)` link+keyword form — and audits each target Issue's ACs from **its own** repo (`fetch_issue(..., repo=owner/repo)`), blocking on any unticked box. It **single-sources** `closure_audit.collect_cross_repo_ac_gaps` (reusing the same `check_ac` the same-repo gate's bash mirrors), honors the `REPO` override so it composes with the sibling gates (#607), excludes the PR's own repo (native references cover same-repo), and **fails open** (exit 0 + warning) on a missing interpreter or `gh` error. The keyword↔ref association is line-scoped (a deliberate over-audit: an extra target whose ACs are ticked passes anyway) — to reference a cross-repo Issue *without* gating on it, keep the closing keyword off that line. Companion gap 1 of #665 (bot-mention guard cross-repo canonicalization) is MM-coordinated and ships separately.

**Issue-authoring convention (the rule the lint backstops):** put **deliverable-gating** checkboxes under a canonical `## Acceptance criteria` heading. A phased plan is fine — phases may live *inside* that section, or as a separate narrative `## Plan` section — but the boxes the closure ritual gates on are the **AC** ones. Non-AC checklists (plan steps, exploratory tasks) are intentionally exempt from the gate. (See `feedback_issue_creation_canonical_sections.md`.)

**Workaround for genuine deferrals:** comment-defer per `.agents/memory/shared/feedback_closure_ritual.md` by ticking `- [x]` with a link to the carrier Issue (or remove the line entirely). The script does NOT parse "deferred" inline tags — keep the tick/remove convention explicit.

**Workaround for a non-closing `#N` reference:** break the keyword→`#N` adjacency — neutral phrasing like "related to Issue #N", or move the keyword away from the token — in the PR body **and** commit messages. See `shared/feedback_hash_numbers.md` (Companion foot-gun section). If the close is genuinely intended, add the Issue to `closingIssuesReferences` (link it in the PR body) so the gate recognizes it.

## Branch naming

Branches follow the canonical `<type>/<role>/issue-<N>-<slug>` pattern (Conventional-Branch format + a `role` segment the board/closure automation consume). Create them with **`scripts/new_branch.sh`** ([Issue #578](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/578)) rather than by hand — raw `gh issue develop` auto-slugs the full Issue title into long/mangled names, and bare `git checkout -b` drops the Development-panel link.

```bash
scripts/new_branch.sh <issue#> <short-slug> [--type T] [--role R] [--dry-run]
#   new_branch.sh 578 branch-helper              → feat/pm/issue-578-branch-helper
#   new_branch.sh 578 branch-helper --type spike → spike/pm/issue-578-branch-helper
scripts/new_branch.sh --no-issue <type> <role> <short-slug> [--dry-run]   # issueless fallback
```

- **type** derives from the Issue title's Conventional-Commit prefix (`feat(scripts): …` → `feat`); `--type` overrides; no prefix and no `--type` → error (never guessed).
- **role** comes from the Issue's `role:<x>` label; exactly one → used, multiple → `--role` required, none → `--role` required.
- **slug** is a **short, human-supplied** kebab descriptor (positional arg), sanitized (lowercase, non-`[a-z0-9-]` dropped); omitting it errors with the Issue title + a suggested invocation — it is never auto-slugged from the title.
- **Wrap, not replace `gh issue develop`** (recorded decision): the helper feeds `gh issue develop` our canonical `--name`, preserving the Issue↔branch Development-panel link (there is no retroactive link). `--no-issue` is the *only* sanctioned `git checkout -b` path.
- **Parent guard:** the script refuses an Issue with `subIssuesSummary.total > 0` (epic) — branching off a parent creates a `closingIssuesReferences` edge that auto-closes it on merge. It replicates the `gh issue develop` PreToolUse hook, which can't see the nested `gh` call. Branch off a leaf sub-issue instead.
- `--dry-run` prints the computed name and exits without any git/gh mutation.

This is the interim convention until the helper is habitual; the underlying "Issue-tracked branches use `gh issue develop`, never `git checkout -b`" rule is enforced by the parent-guard hook above.

## Merge workflow

For all `gh pr merge` invocations, prefer the closure-ritual gate:

```bash
bash scripts/audit_and_merge.sh <PR_NUMBER> [--squash|--merge|--rebase] [--delete-branch|--no-delete-branch] [--skip-review-offer]
```

Defaults: `--squash --delete-branch`. The script audits the PR body Test plan + every linked Issue's Acceptance criteria for unticked `- [ ]` boxes, prints any gaps to stderr, exits 1 without merging if any remain. On a clean audit it then runs the bot-review-offer gate (offer a `@-claude review` first, or pass `--skip-review-offer` for a trivial PR) before forwarding to `gh pr merge` with the chosen flags.

This is the operational path for shipping; bare `gh pr merge` bypasses the closure-ritual gate and should not be used.
