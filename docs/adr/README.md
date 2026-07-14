# Architecture Decision Records (`docs/adr/`)

This directory holds **technical decision records**: point-in-time "we chose X because Y"
decisions about the pipeline's architecture, infrastructure, and tooling.

An ADR captures *why* a choice was made so a future reader (or agent) does not silently
undo it. Each ADR is a numbered, mostly-immutable file: once Accepted, a decision is
changed by writing a *new* ADR that supersedes the old one (mark the old one
`Superseded by [ADR-XXXX]`), not by editing history.

## Format

Lightweight [MADR](https://adr.github.io/madr/)-style: **Context / Decision / Consequences**
(+ optional *Alternatives considered*). Copy [`0000-adr-template.md`](0000-adr-template.md)
and give it the next free number: `NNNN-short-slug.md`.

## Index

| ADR | Title | Status |
|-----|-------|--------|
| [0001](0001-tcrdock-via-docker.md) | Run TCRdock in a Docker container, not a conda env | Accepted |
| [0002](0002-cloud-storage-vs-compute.md) | Treat cloud storage and compute differently after the GCP decommission | Accepted |
| [0003](0003-converge-on-star.md) | Converge on STAR as the single aligner, run locally via a linux-64 container | Accepted |

## Why this home exists (and the routing rule)

This directory and its sibling [`docs/design/`](../design/) exist to unload the project-root
`CLAUDE.md`, which had grown to mix three [Diátaxis](https://diataxis.fr/) content types in
one file. The cut is **by audience and content type**:

| Content type | Home | Examples |
|--------------|------|----------|
| **Decision** ("we chose X because Y") | `docs/adr/` | TCRdock-via-Docker, torch `cu126`/P100 pin, GCP zone migration, regtools/libdeflate `samtools` workaround, `assembly:` config removal |
| **Explanation** ("here is how/why this works") | [`docs/design/`](../design/) | junction-origin classification, PDB chain relabelling, Snakemake-8 gotchas (conceptual), BED12 / STAR `SJ.out.tab` semantics |
| **Agent reference** (commands, env matrix, hook behavior, "do X not Y" gotchas) | **stays in `CLAUDE.md`** | the env matrix, GitHub safety-wrapper hooks, merge workflow |

When a decision block moves here, `CLAUDE.md` keeps a **terse one-line pointer** back to
the ADR so an agent reading the reference file still finds the rationale.

This is the technical sibling of the board-governance-prose extraction tracked in
[Issue #769](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/769) (PM-owned,
`arc:board-governance`). To avoid the two extractions re-overlapping: **board / project /
Kanban governance prose is #769's territory (its own canonical home); technical design
rationale is this directory's.** A block about *how the board works* goes to #769's home,
not here.
