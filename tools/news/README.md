# Morning-news release-feed poller (Issue #639)

Deterministic engine for the morning-routine **News phase**. Instead of asking
the web "what's new with Snakemake?" (which re-surfaces the same evergreen
listicles every day), this polls our actual dependency feeds and reports a tool
**only when its upstream version moved since we last reported it**.

This is the project-repo half of [Issue #639](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/639).
The companion News-phase memory-rule changes (suppress-tracked-instead-of-annotate;
recency-shaped discovery search) live in `shared/feedback_morning_routine.md` and
are landed separately by MM.

## Why it kills the repetition

The morning briefing used to re-surface Snakemake 9, the Astral stack, etc. for
weeks because the News phase had **no cross-day memory** after `news_log.md` was
retired ([Issue #484](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/484)),
and the dedup step only *annotated* tracked items rather than suppressing them.

The watermark restores cross-day memory in a **minimal, single-writer-per-role**
shape that sidesteps both reasons #484 retired the old log (3-role write
contention; a brittle 7-day time window). It stores only the **last version
reported per tool** — a high-water mark, not a log.

## Usage

```bash
conda activate snakemake          # provides python + PyYAML; gh is on PATH
python tools/news/poll_releases.py --role developer
```

Typical morning output (only genuine deltas appear):

```
## Dependency deltas (since last briefing)
- **snakemake** — 9.21.0 → 9.22.0  (tracked: Issue #200)
- **samtools** — 1.22.1 → 1.23.1

_Watch-list (candidate-to-adopt, muted): uv, ruff_
```

Flags: `--no-write` (don't persist the watermark), `--no-check-tracked` (skip the
open-Issue annotation lookup), `--json` (machine-readable), `--watermark PATH`,
`--basket PATH`.

## How it works

1. **`tools.yaml`** — the basket: each dependency with its `feed_type`
   (`pypi` / `github` / `bioconda` / `pytorch_cu126` / `manual`) and any guards.
2. **`watermark_<role>.yaml`** (gitignored) — `{tool: last_version}` + `last_briefing`.
3. **`poll_releases.py`** — fetch latest per feed → diff against the watermark →
   surface only changed tools → write the watermark back.

### Pure-delta semantics

- A tool surfaces **once** when its version moves, then goes quiet until it moves
  again. The horizon is *"since I last reported it"* (unbounded) — **not** a time
  window. A multi-day gap yields one clean `9.21 → 9.25` delta, not a stale repeat.
- A genuine bump surfaces **even if a board Issue already tracks it** (it gets
  annotated with the Issue number); the watermark mutes the *repeat*, not the
  first sighting. (Decision recorded in #639: "pure delta", not suppress-tracked.)
- **First run is silent**: it seeds every tool's baseline at the current upstream
  version. So this answers *"what's new since we adopted the poller"*, **not**
  *"are we behind our pins"* — an upstream-vs-pin **drift audit** is a separate,
  out-of-scope concern (e.g. OptiType upstream `1.5.0` vs our `=1.3.5` pin is not
  reported here).

### Guards (avoid false "you're behind" noise)

| Guard | Tools | Behavior |
|-------|-------|----------|
| `frozen: true` | `jax`, `TCRdock` | Never flagged; never enters the watermark. Pinned together for the AlphaFold-2.3.2-era Docker compat. |
| `max_version: X` | `IMGTgeneDL` (`0.6.1`) | Never flagged above the cap — `>=0.7.0` breaks the conda solve. |
| `watch_only: true` | `gcp_dlvm_image` | Surfaced as a *watch* note, never an actionable bump — a newer DLVM image is a P100 regression risk ([Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522)), not a win. |
| `pytorch_cu126` feed | `torch` | Polls the cu126 wheel channel (not plain PyPI) and preserves the `+cu126` tag — a bare PyPI poll always shows a newer non-Pascal build. |

## Tests

```bash
workflow/tests/.venv/bin/python -m pytest tools/news/ -v
```

All unit tests are **fully offline** — network and `gh` calls are
dependency-injected. CI runs them in the `ci-tools-pytest` job. Live-feed checks
would carry the `live` marker (opt-in via `-m live`); none ship by default.
