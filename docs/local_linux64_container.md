# Running the pipeline locally on macOS arm64 (linux-64 container)

Implements [ADR-0003](adr/0003-converge-on-star.md).
This is a **local-dev affordance only**. On a native linux-64 host (CI, any funded VM) you do not need any of this: run `snakemake --use-conda` directly.

## The one-time setup

```bash
brew install colima docker
colima start --vm-type=vz --vz-rosetta --cpu 4 --memory 4 --disk 30
```

Colima is Apache-2.0 and costs nothing, which is why it was chosen over Docker Desktop under the $0 posture.

**Do not pass `--arch x86_64`.**
That emulates the *entire VM* in qemu and is far slower.
The form above keeps the VM native aarch64 and lets Rosetta translate only the `linux/amd64` container binaries.

## Running the chr22 test through the STAR path

```bash
bash scripts/run_local_linux64.sh --cores 4 --use-conda \
  --configfile config/test_config.yaml config/test_star_config.yaml
```

The wrapper builds `docker/Dockerfile.local` on first use, starts Colima if the daemon is down, and forwards everything else to `snakemake`.
Both config files go in a **single** `--configfile` invocation; Snakemake deep-merges the later file over the earlier one, so the STAR overlay changes only `alignment.aligner` and inherits the rest.
Passing `--configfile` twice does not work (see CLAUDE.md, "configfile flag collapsing").

To run the HISAT2 path instead, drop the overlay and go native; no container is involved.

## Why a container, and not simply "Linux"

The load-bearing property is **`x86_64`**, not Linux. This is the part that is easy to get wrong:

- bioconda ships **no `osx-arm64` STAR** at our pinned version, and *both* macOS builds (native arm64 and Rosetta osx-64) build the chr22 index correctly but then read **zero reads** from every FASTQ, including a synthetic read cut from the chr22 reference itself. The index generator works; the read parser returns instant EOF.
- bioconda also ships **no `linux-aarch64` STAR**. A dry-solve fails with `star =2.7.10b does not exist`. So a *native* Linux VM on an M1 (which is what Colima gives you by default) hits the same wall one platform over.

A `linux/amd64` container, translated by Rosetta, is therefore the cheapest x86_64 Linux available on Apple Silicon.
Verified end-to-end through the full two-pass DAG on 2026-07-15 (Issue #1162): both chr22 samples aligned (tumor 3,379 / normal 3,105 raw junctions, 500,000 input reads each), 388 classified records, index 407 MB, about 3 minutes per sample on 4 threads.
The two-pass run needs the overcommit fix in the gotchas below; without it STAR dies at 2nd-pass start.

## Why Snakemake runs *inside* the container, rather than calling it per-rule

Snakemake has **no Docker or Podman execution backend**. Its `container:` directive is executed by Apptainer/Singularity, and `docker://` is only an image address that Apptainer pulls and converts. Apptainer has no native macOS build, so the documented macOS route is Apptainer-inside-a-Lima-VM, which means Snakemake ends up running inside a Linux VM anyway.

Given that, running Snakemake itself inside one `linux/amd64` container is the smaller moving-part count for the same result, and it leaves the Linux path completely untouched. Shelling out to `docker run` from inside a rule was considered and rejected: it is a hand-rolled backend sitting outside the tool's own model.

**Follow-up worth taking later:** `snakemake --containerize` generates a Dockerfile from the workflow's existing conda envs, which would pin the whole software stack into one image and give local, CI, and production a byte-identical environment. That is a strict improvement on this setup, not a replacement for it.

## Gotchas that cost real time

- **STAR's two-pass mode vforks `zcat` at 2nd-pass start, and a swapless VM refuses the fork.** The colima `vz` VM has no swap, so once STAR's genome and suffix array are resident, `--readFilesCommand zcat` fails with `Failed vforking readFilesCommand / 12: Cannot allocate memory` under Linux' default heuristic overcommit, even though a `vfork` shares the parent's address space and immediately execs a tiny process. It is **not** the `--limitSjdbInsertNsj` cap: the two-pass junction insertion completes *before* this fails. `scripts/run_local_linux64.sh` sets `vm.overcommit_memory=1` in the VM on every run to grant the reservation. If you start Colima and call `snakemake` by hand instead of through the wrapper, set it yourself: `colima ssh -- sudo sysctl -w vm.overcommit_memory=1`.
- **Colima only mounts `$HOME` and `/tmp/colima` into the VM.** A `-v` bind of any other host path silently produces an **empty directory** inside the container. The tool then writes into the container's ephemeral filesystem, exits 0, and the host sees nothing. The wrapper refuses to run if the repo is outside `$HOME` for exactly this reason.
- **Never send the tool's log to `/dev/null` while probing.** The failure above looked like a clean success precisely because the evidence had been discarded.
- Conda envs built inside the container live in `.snakemake/conda-linux64`, deliberately separate from the default `.snakemake/conda` prefix a native macOS run would use. The two hold incompatible binaries for the same env hashes; keeping them apart stops a host run and a container run from corrupting each other.
- **The image is built once and never auto-rebuilt.** The wrapper skips the build whenever the tag exists (`docker image inspect`), so if `docker/Dockerfile.local` later changes, an existing clone keeps the **stale** image. Force a rebuild with `docker rmi splice-neoepitope-local:linux64` (or set a fresh `SNP_LOCAL_IMAGE` tag) after editing the Dockerfile.
