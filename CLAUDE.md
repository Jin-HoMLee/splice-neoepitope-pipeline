# Splice Neoepitope Pipeline — Project Notes

## Instructions for Claude
Update this file only when something is **not derivable from the code or git history** — non-obvious decisions, known gotchas, workarounds, and infrastructure facts. Do not use it as a change log; that's what `git log` is for.

## Project Overview
Modernised reimplementation of a 2015 cancer neoepitope prediction pipeline (Jin-Ho Lee, Seoul National University). Identifies tumor-specific splice junctions from RNA-Seq data and predicts MHC-binding neoepitopes.

## Infrastructure
- Running on GCP Compute Engine VMs — see `docs/google_cloud_guide.md` for full setup
- Current production VMs: `neoepitope-pipeline` (n1-highmem-8 + P100, Phase 1), `pipeline-spot-gpu` (n1-standard-4 + P100, Phase 3); zone `europe-west1-b` (us-central1 has been exhausted in the past)
- `neoepitope-orchestrator` (e2-micro) — lightweight companion VM that starts and manages the pipeline VM in detached mode; stays running cheaply between pipeline runs
- GCS bucket: `gs://splice-neoepitope-project` — results at `.../results/<patient_id>/`, logs at `.../logs/`
- `run_cloud_gpu.sh` defaults to the current local branch; the VM git-pulls it automatically — no `--branch` flag needed unless deliberately running a different branch on the VM
- **NVIDIA driver pinned to `nvidia-headless-570-server` (DKMS)** — do not upgrade. Driver ≥575 dropped P100 Pascal (SM 6.0) support. Image family `common-cu129-ubuntu-2204-nvidia-580` is used but the driver is overridden to 570 in the setup script.
- Pipeline is run with `snakemake --cores $(nproc) --use-conda --rerun-triggers mtime` inside a `tmux` session
- GitHub project board: user project #9 ("JH M Lee Lab") under user `Jin-HoMLee` — query via `gh api graphql` with `user(login: "Jin-HoMLee") { projectV2(number: 9) { ... } }` (it's a user project, not org)
- `main` branch protection: required CI checks (`pipeline-pytest`, `pipeline-snakemake-dry-run`) + squash-merge default. The "Require branches to be up to date before merging" rule was **removed 2026-05-09** — it fired on every PR cut from a worktree branch lagging `main`, even without real conflicts. Don't suggest `gh pr update-branch` or `--admin` workarounds for "branch is behind main" anymore (the rule is gone). Real file-level conflicts still need manual resolution.

### Expected unavailability — P100 in `europe-west1-b`

P100 capacity in `europe-west1-b` has gone sustained-exhausted multiple times — not a transient blip:

- 2026-05-06: 11 launch attempts over ~7h16m, all `ZONE_RESOURCE_POOL_EXHAUSTED` on `n1-highmem-8 + nvidia-tesla-p100` (10:32 → 17:48 BST)
- 2026-05-08: capacity probe still exhausted at ~16:00 UTC; the rolling outage spanned ~46h across 05-06 → 05-08

**Mitigation:** retry overnight or next morning — capacity has typically recovered after a sustained gap. Longer-term fix tracked in [Issue #285](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/285) (P100 contingency epic) and [Issue #310](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/310) (T4/L4 hybrid fallback in `run_cloud_gpu.sh`, blocked on a Google T4 quota grant).

**Last verified:** 2026-05-08.

## Pipeline Design Decisions

### Junction origin classification
Normal samples are used to filter tumor junctions at the junction level, not at the prediction level. The hierarchy:
```
all junctions
  └─ annotated        (in GENCODE)            → discard
  └─ unannotated      (not in GENCODE)
       ├─ normal_shared  (also in normal)  → discard (kept in TSV for reference)
       └─ tumor_exclusive    (absent in normal) → neoepitope prediction
```
This is the clinically correct approach: a junction present in matched normal tissue is not tumor-specific and should not be a neoepitope target. The Fisher's exact test (end-of-pipeline statistical comparison) was removed in favour of this upstream filtering step.

When no normal sample is present, all unannotated junctions are labeled `tumor_exclusive` with a warning — the pipeline still runs.

### TCRdock via Docker
TCRdock runs inside a Docker container (`docker/Dockerfile.pipeline`) rather than a conda env. The conda approach failed due to irreconcilable cuDNN/JAX/openmm version conflicts. The Docker image bundles CUDA 11.8, cuDNN 8, Python 3.10, JAX 0.3.25, AlphaFold params, and BLAST — the host only needs the NVIDIA Container Toolkit. Running CUDA 11.8 inside the container on a host with a newer driver (e.g. 12.8) is supported by NVIDIA's forward-compatibility guarantee.

### PDB chain relabelling
AlphaFold outputs all residues as a single chain (A). `relabel_pdb_chains()` in `run_tcrdock.py` reassigns chain IDs (A=MHC, B=peptide, C=TCR-α, D=TCR-β) using per-chain sequence lengths from TCRdock's `alphafold_setup/targets.tsv`. The report injects PDB COMPND records so Mol* displays meaningful chain names in the sequence panel instead of "Polymer 1/2/3/4".

## Experiment notebooks live under `research/experiments/`

Per-Issue experimental work (analysis notebooks + their cached outputs + a slide deck + a one-page README) lives at `research/experiments/issue_NNN_<short-content-desc>/`. The slide deck (`slides.qmd`) is co-located here rather than under `research/slides/` — see "Slide decks for experiment Issues" below.

Layout per experiment:

```
research/experiments/issue_NNN_<short>/
├── README.md          # one-page: goal, parent issue link, status, outputs index, cross-experiment deps
├── notebook.ipynb     # the analysis
└── outputs/           # cached artifacts (parquet, tsv, png)
```

**Distinguish from `research/notebooks/`:** that folder holds stable per-patient analyses (`patient_001_results.ipynb`, `patient_002_results.ipynb`) — long-lived, manuscript-supporting. The experiments/ folder holds scoped per-Issue work that lands once and then becomes a frozen reference. A separate migration Issue will move existing per-Issue work (`issue_224_*`, `issue_299_*`) from `research/notebooks/` into `research/experiments/`; the convention is established by #225.

### Cross-experiment data sharing

1. **Default:** each experiment owns its outputs in `<experiment>/outputs/`.
2. **Shared between ≥2 experiments → `research/experiments/_shared/`.** Promote an artifact here only when a 2nd consumer materializes (YAGNI; don't pre-share). Filenames carry provenance (`gtex_panel_chr22_snaptron_v1.parquet`, not `gtex_panel.parquet`).
3. **Cross-experiment read, single consumer → explicit path reference.** Document in the consumer's README under "Cross-experiment deps".
4. **Promoted to production → `resources/`.** When an artifact becomes a stable pipeline input.

**Bad practices (any size):** symlinks across experiments, copying artifacts, one experiment writing into another's `outputs/`.

### Size guidance

- **< 10 MB:** check into git.
- **10–100 MB:** check into git, commit a regenerator script alongside.
- **> 100 MB:** keep out of git; store in `gs://splice-neoepitope-project/experiments/<issue>/`; commit a `data_manifest.yaml` in the experiment folder listing artifact paths + checksums + fetch commands. Pin the manifest schema when the first artifact crosses 100 MB.

## Slide decks for experiment Issues

Every experiment-tier Issue ships a Quarto slide deck alongside the per-patient notebook + manuscript work. Decks live **co-located with the notebook** at `research/experiments/issue_NNN_<short-content-desc>/slides.qmd`. Shared scaffolding (`_template.qmd`, `nature.csl`) stays centralized at `research/slides/` and is referenced via `csl: ../../slides/nature.csl` from each deck. See [`research/slides/README.md`](research/slides/README.md) for the full convention (Quarto rationale, render commands, figure-source pattern, Zotero linkage, install). Scope: **one deck per experiment Issue**, not per sub-issue; no decks for closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Figures regenerate from `research/experiments/issue_NNN_<short>/outputs/*.parquet` (notebook outputs) and from a local `figures/_regenerate_figures.py` (deck-only figures) so the **notebook stays canonical**. Rationale for co-location: notebook + outputs + deck rename / archive / migrate as a unit; figure paths shorten from `../../experiments/issue_NNN/outputs/...` to `outputs/...`. (`research/slides/issue_393_alphagenome_chr22_poc/` predates this rule and is migrated under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).)

**Render tooling:** Quarto via `brew install --cask quarto` (macOS dev). PDF render (beamer) is currently disabled in decks — Quarto/pandoc's auto-emitted preamble hits a `\makesavenoteenv{longtable}` interaction with `footnotehyper.sty` that fails LaTeX compile. HTML reveal.js is the primary delivery; for a PDF handout, "Print → Save as PDF" from Chrome works on `slides.html`. A proper LaTeX fix is tracked as a follow-up.

## Slide decks for eval Issues

Every tool-evaluation Issue (e.g. [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES, [Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) ImmSET, [Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) Boltz-2) also ships a Quarto deck — a **tool primer** distinct from the consolidating publication deck. Decks live at `research/evals/issue_NNN_<tool>/slides.qmd`, sibling to `research/experiments/`. Shared scaffolding (`_template.qmd`, `nature.csl`) at `research/slides/` is referenced via `csl: ../../slides/nature.csl`.

**Why required even though evals are XS/S-sized:** visual condensation per tool is a different artifact from the consolidating decision deck (e.g. [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432) TCR-pMHC scorer landscape). The per-tool deck answers *"what is this tool, and how does it work?"* — a comprehension precondition that lab notebook prose alone doesn't satisfy. The consolidating deck answers *"why this (a/b/c) decision across all 5 evals?"*. Both have value.

**Format (~8-10 slides; title + 7-9 content):** (i) the question (should we integrate?), (ii) tool primer (architecture, citation, code/weights status), (iii) why-it-works mechanism, (iv) integration map (block diagram showing where it'd plug in), (v) reasons-to-(a)/decline-(b)/skip-(c), (vi) decision + sub-issue ref, (vii) open scientific questions for the integration (if (a)), (viii) references. [`research/evals/issue_218_hermes/`](research/evals/issue_218_hermes/) is the reference example (9 slides total).

**Scope:** one deck per eval Issue. Established 2026-05-26 alongside the HERMES + ImmSET evals. Audience: lab seminar / external talk + manuscript-figure rehearsal for the consolidating publication Issue.

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
| Pytest test suite | `workflow/tests/.venv` | `pyenv local 3.13.5 && python -m venv workflow/tests/.venv && workflow/tests/.venv/bin/pip install -r workflow/tests/requirements-test.txt` | [workflow/tests/README.md](workflow/tests/README.md) |
| Research / Jupyter notebooks | `research/.venv` | `cd research/ && pyenv local 3.14.4 && python -m venv .venv && .venv/bin/pip install -r requirements.txt` | [research/README.md](research/README.md) |
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
- **Conda solver failure** — e.g. `IMGTgeneDL>=0.7.0` pinned when PyPI max is `0.6.1`. Envs are built lazily on first execute.
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

The production `resources/` directory is no longer used. `setup_vm.sh` performs a one-shot idempotent migration on existing VMs (moves old `resources/*` → `references/` and `resources/{hisat2,star}_index/` → `indices/{hisat2,star}/`). MHCflurry sentinel moved from `resources/mhcflurry_models.done` to `<mhcflurry.models_cache_dir>/.download_done` (default `~/.mhcflurry/`); the actual model cache still lives in MHCflurry's platformdirs default.

## HISAT2 Index Cache Invalidation

Snakemake skips the index download if `indices/hisat2/` already exists (it checks for an `index.done` sentinel file). Changing `hisat2_prebuilt_url` in `config/config.yaml` does **not** invalidate this cache — the old index silently persists and will be used on the next run.

**When changing `hisat2_prebuilt_url`:** delete the index directory on the VM before running:

```bash
gcloud compute ssh neoepitope-pipeline --zone=europe-west1-b --tunnel-through-iap \
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

### `python.yaml` — PyTorch SM 6.0 / P100 compatibility
The conda-forge `pytorch` build pulled into `python.yaml` ships against CUDA 12.8/12.9, which dropped SM 6.0 (Pascal). On a P100, `torch.cuda.is_available()` still returns `True` but kernel dispatch fails silently or with a cryptic error.
**Current workaround:** pin `torch>=2.0,<2.5` in `python.yaml` (installs 2.4.1, which is built against an older CUDA that still includes SM 6.0 kernels).
**The pin can be lifted via the cu126 pip channel — verified empirically [Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352) 2026-05-27.** Pascal removal happened only in CUDA 12.8/12.9 binaries (per the [PyTorch dev-discuss thread](https://dev-discuss.pytorch.org/t/cuda-toolkit-version-and-architecture-support-update-maxwell-and-pascal-architecture-support-removed-in-cuda-12-8-and-12-9-builds/3128)); CUDA 12.6 wheels remain Pascal-compatible through `torch 2.12` (build matrix in [pytorch/pytorch#178665](https://github.com/pytorch/pytorch/issues/178665) — `Pascal(6.0)` listed for cu126/2.12). Verified on `neoepitope-pipeline` 2026-05-27: `torch==2.12.0+cu126` dispatches ReLU on a P100 cleanly (`tensor([0., 0.], device='cuda:0')`). To lift the pin: refactor `python.yaml` to pip-install torch from `https://download.pytorch.org/whl/cu126` instead of conda-forge, and verify MHCflurry 2.2.x still works with the bumped version. Maxwell/Pascal/Volta remain "feature-complete with no further enhancements planned" so this isn't free forever — but is good for the cu126 channel's support window.
`_has_gpu()` in `run_mhcflurry.py` uses a PyTorch smoke-test kernel (not TensorFlow) to catch this case: `torch.nn.functional.relu(torch.zeros(2, device="cuda"))`. TF reported GPU available even when PyTorch kernels would fail — both must work because MHCflurry 2.2.x uses PyTorch for inference.

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

## GitHub Safety Wrappers

Mechanisms that fire automatically to enforce GitHub-related discipline rules that have broken repeatedly despite being documented in memory. Per the mechanism-over-memory ladder (memory → inline Always-in-effect → mechanism), these are the rung-3 escalation when a rule has slipped ≥2× on the same shape.

### `@claude` mention guard ([PreToolUse hook](.claude/settings.json))

Refuses any `gh (pr|issue) (comment|create)` whose `--body` contains a literal `@claude` substring, except the exact canonical review-trigger `--body "@claude review"` (and the single-quoted variant). The hook runs `.claude/hooks/check_at_claude.py` on stdin-piped PreToolUse JSON and emits a `permissionDecision: deny` when the guard fires.

**Why:** the `Claude Code` GitHub Action subscribes to `issues` events and triggers on ANY literal `@claude` — including inside parens, ACs, code spans, quoted historical refs, or markdown link descriptions. The rule lives in `.claude/memory/shared/feedback_no_at_claude_mention.md` and is inlined into shared Always-in-effect, yet broke on PR #359 (2026-05-13) and originally on Issue #272 (2026-05-06).

**Workaround for non-trigger references:** use `@-claude` (zero-width hyphen between `@` and `claude`). The literal substring `@claude` must not appear.

### Closure-ritual gate (`scripts/audit_and_merge.sh`)

Refuses to run `gh pr merge` if any `- [ ]` remains on the PR body Test plan OR on any Issue in `closingIssuesReferences` (under that Issue's Acceptance criteria). Operational usage lives under [Merge workflow](#merge-workflow); this section captures the mechanism-over-memory rationale.

**Why:** the closure-ritual rule has broken 4× in 10 days despite being inlined into shared MEMORY.md — [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280), [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299), [PR #328](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/328), [Issue #347](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/347). Memory of a declarative rule cannot reliably survive the action-distance from session start to `gh pr merge`. The script enforces at the exact moment of action — the merge invocation itself — so the gate cannot be forgotten regardless of how much downstream context has accumulated.

**Workaround for genuine deferrals:** comment-defer per `.claude/memory/shared/feedback_closure_ritual.md` by ticking `- [x]` with a link to the carrier Issue (or remove the line entirely). The script does NOT parse "deferred" inline tags — keep the tick/remove convention explicit.

## Merge workflow

For all `gh pr merge` invocations, prefer the closure-ritual gate:

```bash
bash scripts/audit_and_merge.sh <PR_NUMBER> [--squash|--merge|--rebase] [--delete-branch|--no-delete-branch]
```

Defaults: `--squash --delete-branch`. The script audits the PR body Test plan + every linked Issue's Acceptance criteria for unticked `- [ ]` boxes, prints any gaps to stderr, exits 1 without merging if any remain. On a clean audit it forwards to `gh pr merge` with the chosen flags.

This is the operational path for shipping; bare `gh pr merge` bypasses the closure-ritual gate and should not be used.
