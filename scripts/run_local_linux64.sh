#!/usr/bin/env bash
# Run Snakemake inside a linux-64 container on a macOS arm64 dev box (Issue #1162,
# ADR-0003). Everything after `--` is passed straight to snakemake.
#
#   bash scripts/run_local_linux64.sh --cores 4 --use-conda \
#        --configfile config/test_config.yaml config/test_star_config.yaml
#
# On a native linux-64 host you do NOT need this: call snakemake directly.
#
# Why a container and not a Linux VM: bioconda has no linux-aarch64 STAR, so a
# native aarch64 VM fails exactly like macOS arm64 does. The load-bearing
# property is x86_64. See docker/Dockerfile.local for the full rationale.
set -euo pipefail

IMAGE="${SNP_LOCAL_IMAGE:-splice-neoepitope-local:linux64}"
DOCKERFILE="docker/Dockerfile.local"
# Keep linux-64 conda envs out of the default .snakemake/conda prefix: a macOS
# host run and a container run would otherwise race for the same hashed env dirs
# with incompatible binaries in them.
CONDA_PREFIX_DIR="${SNP_CONDA_PREFIX:-.snakemake/conda-linux64}"

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

die() { printf 'error: %s\n' "$1" >&2; exit 1; }

# A leading `--` is accepted but not forwarded: snakemake reads a bare `--` as the
# end of its configfile list, which would silently change how its args parse.
if [ "${1:-}" = "--" ]; then
  shift
fi

command -v docker >/dev/null 2>&1 \
  || die "docker CLI not found. Install with: brew install colima docker"

# Colima is our runtime (Apache-2.0, no licensing question). Start it if the
# daemon is not reachable. A user on Docker Desktop just needs the app running;
# in that case the colima calls are skipped.
if ! docker info >/dev/null 2>&1; then
  if command -v colima >/dev/null 2>&1; then
    printf '>> docker daemon not reachable; starting colima\n' >&2
    # vz + Rosetta: the VM stays native aarch64 and Rosetta translates only the
    # linux/amd64 container binaries. Do NOT use `--arch x86_64`, which emulates
    # the entire VM in qemu and is far slower.
    colima start --vm-type=vz --vz-rosetta --cpu 4 --memory 4 --disk 30
  else
    die "docker daemon not reachable and colima is not installed."
  fi
fi

# STAR's --twopassMode Basic vforks `zcat` (--readFilesCommand) at the start of
# 2nd-pass mapping, once the genome + suffix array are resident. The colima vz VM
# has NO swap, so under Linux' default heuristic overcommit the fork's memory
# reservation is refused and STAR dies "Failed vforking readFilesCommand /
# 12: Cannot allocate memory" - even though the child immediately execs a tiny
# process and never copies the parent (a vfork shares the address space). Grant
# overcommit so the reservation always succeeds. Idempotent; applied every run so
# it holds regardless of who started colima. Issue #1162.
if command -v colima >/dev/null 2>&1 && colima status >/dev/null 2>&1; then
  colima ssh -- sudo sysctl -w vm.overcommit_memory=1 >/dev/null 2>&1 \
    || printf '>> warning: could not set vm.overcommit_memory=1 in the VM; STAR 2-pass may OOM at fork\n' >&2
else
  # Docker Desktop (or any non-colima daemon): we cannot reach the VM to set the
  # sysctl. Warn rather than silently skip - if that VM is swapless too, STAR's
  # 2-pass would OOM exactly as it did on colima before this fix.
  printf '>> note: not on colima; vm.overcommit_memory left unset. If STAR 2-pass dies "Cannot allocate memory" at fork, set vm.overcommit_memory=1 in your Docker VM.\n' >&2
fi

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  printf '>> building %s (first run only)\n' "$IMAGE" >&2
  docker build --platform linux/amd64 -t "$IMAGE" -f "$DOCKERFILE" .
fi

# The bind mount must live under $HOME: colima only mounts $HOME and /tmp/colima
# into the VM, and a -v of any other path silently yields an EMPTY directory
# inside the container (the tool then "succeeds" and writes nothing to the host).
case "$REPO_ROOT" in
  "$HOME"/*) : ;;
  *) die "repo is outside \$HOME ($REPO_ROOT); colima will not mount it. Move the clone under \$HOME." ;;
esac

CONTAINER_HOME="$REPO_ROOT/.snakemake/container-home"
mkdir -p "$CONDA_PREFIX_DIR" "$CONTAINER_HOME"

# Interactive only when there is a real TTY: `-t` aborts under a non-TTY caller
# (CI, an agent shell) with "the input device is not a TTY".
TTY_FLAGS=""
if [ -t 0 ] && [ -t 1 ]; then
  TTY_FLAGS="-it"
fi

# We run as the host uid/gid so files written into the mounted repo are not
# root-owned. That makes the image's own /opt/conda/pkgs unwritable, so conda's
# package cache is redirected into the repo alongside the env prefix.
exec docker run --rm $TTY_FLAGS \
  --platform linux/amd64 \
  -v "$REPO_ROOT":"$REPO_ROOT" \
  -w "$REPO_ROOT" \
  -u "$(id -u):$(id -g)" \
  -e HOME="$CONTAINER_HOME" \
  -e CONDA_PKGS_DIRS="$CONTAINER_HOME/conda-pkgs" \
  -e XDG_CACHE_HOME="$CONTAINER_HOME/cache" \
  "$IMAGE" \
  snakemake --conda-prefix "$CONDA_PREFIX_DIR" "$@"
