## Junction/repeat overlap - normal (chr22)

| set | junctions | splice site in a repeat |
|---|---|---|
| lost to the filter | 198 | 182/198 (91.9%) |
| retained | 1516 | 1266/1516 (83.5%) |
| gained (should be 0) | 0 | - |

Unstratified enrichment among lost junctions: **1.10x** (CONFOUNDED - see below)

### Stratified by annotation status (this is the one to read)

| stratum | junctions | splice site in a repeat |
|---|---|---|
| lost, annotated | 14 | 0/14 (0.0%) |
| lost, unannotated | 184 | 182/184 (98.9%) |
| retained, annotated | 245 | 8/245 (3.3%) |
| retained, unannotated | 1271 | 1258/1271 (99.0%) |

**Enrichment within the unannotated pool: 1.00x** - the unannotated pool is what the filter actually draws from, so this is the number that says whether it targets repeats.

Single-read junctions: 183/198 (92.4%) of lost vs 1349/1516 (89.0%) of retained (a large gap would mean the filter is really just removing low-coverage calls).

| repeat class at splice site | lost | retained |
|---|---|---|
| SINE | 170 | 1200 |
| (none) | 16 | 250 |
| LINE | 10 | 40 |
| Retroposon | 2 | 10 |

