# research/

Scientific work associated with the splice neoepitope pipeline project.

| Path | Contents |
|------|----------|
| `manuscript/` | Manuscript sections (INTRODUCTION, METHODS, RESULTS, DISCUSSIONS, CONCLUSIONS) |
| `lab_notebook.md` | Chronological session log — all decisions, experiments, and findings |
| `scripts/` | Research support scripts (literature management, etc.) |

## Jupyter notebooks (`notebooks/`)

Per-patient result interpretation notebooks. Each patient run gets its own notebook (`patient_NNN_results.ipynb`). Cross-patient comparisons go in `results_comparison.ipynb`.

**Setup (one-time):** requires [`uv`](https://docs.astral.sh/uv/) (`brew install uv` or `curl -LsSf https://astral.sh/uv/install.sh | sh`).
```bash
cd research/
pyenv local 3.14.4            # sets .python-version (gitignored, per-clone)
uv venv --python 3.14.4 .venv                              # fast venv creation
uv pip install --python .venv/bin/python -r requirements.txt  # 10-100x faster than pip
```
Only the two per-clone pyenv venvs use `uv`; the `snakemake` conda env and per-rule `--use-conda` envs are unchanged (`uv` is pip-side only, never touches conda's solver).

**Usage:** Open any `.ipynb` in VSCode and select `research/.venv` as the kernel (Python 3.14.4).

**Headless execution (one-time, per-clone):** register this clone's venv as the `splice-neoepitope-research` Jupyter kernelspec so `jupyter nbconvert --execute` and agent/CI contexts resolve to a real interpreter:
```bash
research/.venv/bin/python -m ipykernel install --user \
  --name splice-neoepitope-research \
  --display-name "Splice Neoepitope (research)"
```
This is a machine-local registration (`~/Library/Jupyter/kernels/` on macOS, `~/.local/share/jupyter/kernels/` on Linux for the GPU-revival path), not tracked in git - like the `.venv` itself, each clone runs it once. The kernelspec name is global, so the last clone to run it becomes the canonical execution clone; run it from whichever clone you execute notebooks in. After registration, run notebooks from this venv's jupyter (`research/.venv/bin/jupyter nbconvert --to notebook --execute <nb>`) with no `--ExecutePreprocessor.kernel_name` override needed.

The notebooks read pipeline results from **Cloudflare R2** (`results/<patient_id>/` on the project bucket) via the [`r2_io.py`](notebooks/r2_io.py) helper, which replaced the decommissioned GCS/`gsutil` loader at the GCP exit ([#854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)). Credentials come from the project-root `.env` (gitignored): `R2_ENDPOINT`, `R2_BUCKET`, `R2_ACCESS_KEY_ID`, `R2_SECRET_ACCESS_KEY`. Reads are locally cached (`$R2_CACHE_DIR`, default `~/.cache/splice-neoepitope-r2`) and ETag-validated, so a re-run is cache-backed and offline-tolerant:

```python
from r2_io import r2_read_tsv_cached
report = r2_read_tsv_cached("results/patient_001/reports/report.tsv")
```

The bare `import r2_io` resolves because notebooks run with `research/notebooks/` as the cwd (Jupyter's default); adjust the import path if you run the snippet from elsewhere.

---

## scripts/zotero_add.py

Adds a paper to the [Splice Neoepitope Pipeline Zotero collection](https://www.zotero.org/) (collection key `Z38GTJNW`) by DOI, and optionally attaches a structured note.

**Setup:** create a `.env` file in the project root with your credentials:
```
ZOTERO_USER_ID=<your user ID>
ZOTERO_API_KEY=<your API key>
```

**Add a new paper:**
```bash
python research/scripts/zotero_add.py "10.1038/s41586-024-08552-0" \
  --tags splice-neoepitope public-neoantigen \
  --note "<p><strong>Results:</strong><br>• ...</p><p><strong>Methods:</strong><br>• ...</p><p><strong>Limitations:</strong><br>• ...</p>"
```

**Update the note on an existing item:**
```bash
python research/scripts/zotero_add.py --update-note ITEM_KEY \
  --note "<p><strong>Results:</strong><br>• ...</p>..."
```

Note format (three sections, max 2 bullets each):
- **Results:** `• [finding] — [why it matters generally] / [why it matters for us]`
- **Methods:** `• [tool/approach] — [reuse potential]`
- **Limitations:** `• [gap] — [how we address it / open question]`

---

## scripts/aggregate_cohort.py

Combines per-patient `report.tsv` files into a single cohort table for cross-patient analysis (Issue #84). Pipeline runs are single-patient by design, so this is a post-pipeline research-time step rather than a Snakemake rule.

Each per-patient `report.tsv` carries the schema `patient_id | stage | metric | value | notes`; aggregation is a sorted concatenation.

**By explicit input paths:**
```bash
python research/scripts/aggregate_cohort.py \
  --inputs results/patient_001/reports/report.tsv \
           results/patient_002/reports/report.tsv \
  --output research/cohort_summary.tsv
```

**By patient IDs (paths resolved as `{results-root}/{ID}/reports/report.tsv`):**
```bash
python research/scripts/aggregate_cohort.py \
  --patients patient_001 patient_002 \
  --results-root results \
  --output research/cohort_summary.tsv
```

Output is sorted by `(patient_id, stage, metric)` for stable diffs.
