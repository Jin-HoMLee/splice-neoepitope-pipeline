# Design explanation docs (`docs/design/`)

This directory holds **design explanation**: narrative "here is how/why this part of the
pipeline works" documents. Unlike an [ADR](../adr/) (a point-in-time decision, mostly
immutable), a design doc is **editable** and evolves with the code it describes.

These are the [Diátaxis](https://diataxis.fr/) *explanation* type: understanding-oriented
prose, not a command reference and not a single dated decision.

## Index

| Doc | Topic |
|-----|-------|
| [junction-origin-classification.md](junction-origin-classification.md) | How tumor-vs-normal junction filtering classifies each junction |

## Relationship to `docs/adr/` and `CLAUDE.md`

See [`docs/adr/README.md`](../adr/README.md) for the full routing rule. In short: a crisp
*decision* ("we chose X because Y") is an ADR; ongoing *explanation* of how something works
is a design doc; operational *agent reference* (commands, env matrix, hooks) stays in the
project-root `CLAUDE.md`.
