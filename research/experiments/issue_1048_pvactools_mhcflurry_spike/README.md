# Issue #1048 - Can pVACtools drive its source modules on MHCflurry-only (class I, zero NetMHC)?

**Feasibility spike.** Positioning + 3-layer architecture: [#1045](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1045).
Sibling audit (open splice callers): [`issue_965_open_caller_audit`](../issue_965_open_caller_audit/README.md).

**Verdict: GREEN.** pVACtools 7.0.1 runs its source modules end-to-end for class I on **MHCflurry only**, with **no NetMHC binary, no IEDB standalone install, no `mhcnuggets`, and no `tensorflow` present anywhere in the environment**.
"Absorb the framework" is real: Layer A is a **wrapper job**, not a partial rebuild.

Run on macOS arm64 (Apple M1), CPU-only, 2026-07-09.

## What was tested

| AC | Check | Result |
|----|-------|--------|
| 1 | Open-only env; NetMHC binaries confirmed absent | **PASS** - `netMHCpan` / `netMHC` / `netMHCIIpan` / `netMHCcons` / `netMHCstabpan` all absent from `PATH`; no `~/mhc_i`, `~/mhc_ii`, `/opt/iedb` |
| 2 | `pvacseq` runs to completion on a VCF, MHCflurry-only class I | **PASS** - 9,132 epitopes scored, 27 pass filters |
| 3 | At least one further source module runs MHCflurry-only | **PASS** - `pvacbind` on a peptide FASTA: 37 epitopes scored, 3 pass filters |
| 4 | GREEN/RED documented; license + dependency gotchas recorded | this file |
| 5 | If GREEN: wrapper surface noted | see "Wrapper surface" below |

`pvacsplice` was **not** run (see "Not tested" below) - it is the module we care about most, and it is the one this spike did *not* settle.

### Biological sanity check (not just plumbing)

The `pvacbind` fixture ([`fixtures/peptides.fa`](fixtures/peptides.fa)) embeds three known HLA-A\*02:01 epitopes inside longer sequences.
Of the 37 scored 9-mers, exactly the three canonical **binders (by IC50)** survive filtering.
"Binder" is the right word here and not a vocabulary slip: pVACtools filters on **affinity** (`ic50`), which is precisely the quantity the house presentation vocabulary distinguishes itself from (see "Wrapper surface" below).

| Epitope | Source |
|---------|--------|
| `GILGFVFTL` | influenza A M1 58-66 |
| `SLYNTVATL` | HIV-1 gag p17 77-85 |
| `KLQEEIPVL` | third seeded binder |

So the pipeline is scoring real presentation, not merely completing without error.
`pvacseq`'s top hit on the bundled HCC1395 example is `TLYRVPLLV` (OSTC, IC50 17.7 nM), attributed to the `MHCflurry MT IC50 Score` column.

## The load-bearing finding: `mhcnuggets` is a hard dependency, and it is the only tensorflow vector

`pvactools==7.0.1` declares **both** of these as hard install requirements (not extras):

```
mhcnuggets==2.4.1      # JHU, NON-COMMERCIAL license
mhcflurry==2.0.6       # Apache-2.0, but PINNED and tensorflow-backed
```

So a plain `pip install pvactools` **cannot** produce an open-only environment: it drags in `mhcnuggets`, which drags in `tensorflow`.

**This is dodgeable, and dodging it is not a hack.** `mhcnuggets` is only *needed* if you select `MHCnuggetsI` / `MHCnuggetsII` as prediction algorithms. Installing pVACtools with `--no-deps` and supplying the dependency list minus `mhcnuggets` yields a working install where `MHCflurry` and `MHCflurryEL` remain selectable. Verified: `mhcnuggets`, `tensorflow`, and `keras` are all absent from [`env_lock.txt`](env_lock.txt), and both modules still run.

Two consequences worth carrying into #1045:

1. **Licensing.** We never install or redistribute the non-commercial `mhcnuggets`. The open-only claim survives, but only because we build the env deliberately. Any naive `pip install pvactools` in a Dockerfile silently reintroduces a non-commercial dependency.
2. **The `mhcflurry==2.0.6` pin is wrong for us and must be overridden.** MHCflurry 2.0.6 lazily `import tensorflow` inside `Class1NeuralNetwork.merge()` at *prediction* time (an undeclared runtime dep - the install succeeds, the prediction crashes). Our pipeline runs **MHCflurry 2.2.x on PyTorch**. Upgrading to `mhcflurry==2.2.1` fixes the crash, removes the tensorflow requirement entirely, and **still works with pVACtools**, because pVACtools invokes MHCflurry as a *subprocess CLI* (`mhcflurry-predict --alleles ... --peptides ... --out ...`), not as a Python API. That CLI contract is stable across 2.0.6 and 2.2.1, and pVACtools already renames both the old (`mhcflurry_prediction`) and new (`mhcflurry_affinity`) output columns.

The subprocess-CLI coupling is the reason this spike is GREEN. pVACtools is not entangled with MHCflurry's internals; it shells out. That also means **our production MHCflurry version is the one that scores**, which is exactly the property the 3-layer architecture wants.

## Gotchas (all environment, none fatal)

| Gotcha | Symptom | Fix |
|---|---|---|
| `mhcflurry` 2.0.6 pin | `ModuleNotFoundError: No module named 'tensorflow'` at predict time | override to `mhcflurry>=2.1` (we used 2.2.1, torch) |
| Bare `mhcflurry-predict` lookup | `Exception: An error occurred while calling MHCflurry:` with an **empty** message | pVACtools calls the bare CLI name; the venv `bin` must be on `PATH`. Invoking `/abs/path/to/pvacbind` is not enough. The empty message is a swallowed `FileNotFoundError` |
| `h5py` undeclared | `ModuleNotFoundError: No module named 'h5py'` importing `pvactools.lib.output_parser` | normally arrives transitively via mhcnuggets/tensorflow; install `h5py` explicitly once mhcnuggets is dropped |
| `setuptools>=81` | `ModuleNotFoundError: No module named 'pkg_resources'` | pin `setuptools<81` (mhcflurry still imports `pkg_resources`) |
| `tkinter` required at import | `ModuleNotFoundError: No module named '_tkinter'` | `pvactools/lib/vector_visualization.py` imports `turtle` at package import. Use a Python built with Tk (our `pyenv 3.10.13`; `pyenv 3.11.5` / homebrew 3.11 / 3.13 lack it) |
| macOS: `shellinford` C++ build | `fatal error: 'vector' file not found` building `vaxrank` -> `shellinford` | **local toolchain defect, not pVACtools.** `/Library/Developer/CommandLineTools/usr/include/c++/v1/` exists but has no `vector` header, so `clang++` cannot compile *any* C++ here. Build with `CPPFLAGS="-isysroot $(xcrun --show-sdk-path) -I$(xcrun --show-sdk-path)/usr/include/c++/v1"` |

## Wrapper surface (AC 5)

If we absorb pVACtools as Layer A, these are the CLI entry points to wrap. Each takes an input, a sample name, an allele list, an algorithm list, and an output dir:

```
pvacseq   run <vep_annotated.vcf.gz> <sample> <alleles> MHCflurry <outdir> -e1 <lengths>
pvacbind  run <peptides.fa>          <sample> <alleles> MHCflurry <outdir> -e1 <lengths>
pvacfuse  run <agfusion_dir>         <sample> <alleles> MHCflurry <outdir> -e1 <lengths>
pvacsplice run <regtools_junctions>  <sample> <alleles> MHCflurry <outdir> ...   # not tested here
```

Notes for the wrapper:
- Algorithm string is literally `MHCflurry` (affinity) or `MHCflurryEL` (eluted-ligand). `pvactools valid_algorithms` lists all 22; the NetMHC family is *listed* but unusable without the binaries, so the wrapper must pin the algorithm rather than expose the full list.
- Output lands at `<outdir>/MHC_Class_I/<sample>.MHC_I.{all_epitopes,filtered,all_epitopes.aggregated}.tsv`.
- `--iedb-install-directory` is **not** required when only MHCflurry is selected. Despite the internal method being named `call_iedb()`, no IEDB call occurs.
- pVACtools emits `ic50` / `percentile`, i.e. **affinity** semantics. Our house vocabulary is presentation (`presentation_score`, `presentation_percentile`). A wrapper that adopts pVACtools' columns must not silently rename them into our presentation vocabulary; they are the affinity-only quantities. Prefer `MHCflurryEL` plus our own `Class1PresentationPredictor` pass for the ranked output.
- This spike ran `-e1 9` (9-mers only), which is sufficient for a feasibility proof but is **not** the production setting. The wrapper must pass the full range from `config.yaml` `translation.peptide_lengths`, currently `[8, 9, 10]`, as `-e1 8 9 10`.

## Not tested (scope boundary, and the honest gap)

- **`pvacsplice`** - the module most relevant to this project. It needs a regtools junction file plus a reference FASTA/GTF, and its published path is mutation-driven (it takes a DNA VCF alongside junctions). Whether it can be fed junction-only, without a VCF, is the open question flagged in the Issue's Risks section and it remains open. This spike deliberately settled the *class-I MHCflurry* question (the primary) and left the pvacsplice input-contract question (the secondary) unanswered.
- **`pvacfuse`** - needs an AGFusion/Arriba directory; not built.
- **Class II** - unchanged and still blocked. Every class-II algorithm pVACtools offers is NetMHCIIpan / MHCnuggetsII / MixMHC2pred / NNalign / SMMalign, i.e. NetMHC-family or non-commercial. Consistent with [#1049](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1049): no license-clean class-II presentation predictor exists.
- **Linux / CI** - only macOS arm64 was exercised. The tkinter and shellinford gotchas are macOS-flavored; a Linux env will hit different ones.

## Reproduce

```bash
bash reproduce.sh              # builds the venv, runs both modules, diffs against outputs/
STRICT=1 bash reproduce.sh     # additionally FAILS if a fresh output diverges from outputs/
```

A plain `bash reproduce.sh` exiting 0 means **"both modules ran open-only"**, not "the outputs were byte-identical" - the committed-vs-fresh diff is advisory by default, because float jitter across MHCflurry/torch builds and platforms is expected and is not a defect. The byte-identical claim below was substantiated with `STRICT=1` on the authoring machine (macOS arm64, `env_lock.txt` versions); reproduce it that way if you want the same guarantee.

Environment is pinned in [`env_lock.txt`](env_lock.txt) (92 packages). Note the lock is a *record*, not a constraint file consumed by `reproduce.sh`; the script pins only the load-bearing versions so the resolver stays free elsewhere.

Committed artifacts:
- [`fixtures/peptides.fa`](fixtures/peptides.fa) - 3 sequences seeding known A\*02:01 epitopes
- [`outputs/pvacbind.MHC_I.filtered.tsv`](outputs/pvacbind.MHC_I.filtered.tsv) - 3 rows
- [`outputs/pvacseq.MHC_I.filtered.tsv`](outputs/pvacseq.MHC_I.filtered.tsv) - 27 rows (pVACtools' bundled HCC1395 example VCF)

**Created by:** Scientist
