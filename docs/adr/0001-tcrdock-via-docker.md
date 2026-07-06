# ADR-0001: Run TCRdock in a Docker container, not a conda env

- **Status:** Accepted
- **Date:** 2026-07-06 (migrated from `CLAUDE.md`; the decision itself predates this record)
- **Deciders:** Developer

## Context

TCRdock needs CUDA, cuDNN, JAX, OpenMM, AlphaFold params, and BLAST at mutually compatible versions.
Provisioning that stack as a conda env failed due to irreconcilable cuDNN/JAX/openmm version conflicts.

## Decision

TCRdock runs inside a Docker container (`docker/Dockerfile.pipeline`) rather than a conda env.
The image bundles CUDA 11.8, cuDNN 8, Python 3.10, JAX 0.3.25, AlphaFold params, and BLAST; the host only needs the NVIDIA Container Toolkit.

## Consequences

- The host needs only the NVIDIA Container Toolkit, not a matched CUDA/cuDNN/JAX toolchain.
- Running CUDA 11.8 inside the container on a host with a newer driver (e.g. 12.8) is supported by NVIDIA's forward-compatibility guarantee.
- TCRdock is decoupled from the conda-managed per-rule envs: it is invoked as a container, not a Snakemake `--use-conda` rule.

## Alternatives considered

- **conda env for the whole TCRdock stack** - rejected: irreconcilable cuDNN/JAX/openmm version conflicts.
