# research/

Scientific work associated with the splice neoepitope pipeline project.

| Path | Contents |
|------|----------|
| `manuscript/` | Manuscript sections (INTRODUCTION, METHODS, RESULTS, DISCUSSIONS, CONCLUSIONS) |
| `lab_notebook.md` | Chronological session log — all decisions, experiments, and findings |
| `scripts/` | Research support scripts (literature management, etc.) |

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
