## Junction/repeat overlap - tumor (chr22)

| set | junctions | splice site in a repeat |
|---|---|---|
| lost to the filter | 267 | 248/267 (92.9%) |
| retained | 1605 | 1325/1605 (82.6%) |
| gained (should be 0) | 0 | - |

Unstratified enrichment among lost junctions: **1.13x** (CONFOUNDED - see below)

### Stratified by annotation status (this is the one to read)

| stratum | junctions | splice site in a repeat |
|---|---|---|
| lost, annotated | 11 | 0/11 (0.0%) |
| lost, unannotated | 256 | 248/256 (96.9%) |
| retained, annotated | 276 | 9/276 (3.3%) |
| retained, unannotated | 1329 | 1316/1329 (99.0%) |

**Enrichment within the unannotated pool: 0.98x** - the unannotated pool is what the filter actually draws from, so this is the number that says whether it targets repeats.

Single-read junctions: 229/267 (85.8%) of lost vs 1409/1605 (87.8%) of retained (a large gap would mean the filter is really just removing low-coverage calls).

| repeat class at splice site | lost | retained |
|---|---|---|
| SINE | 234 | 1258 |
| (none) | 19 | 280 |
| Retroposon | 12 | 11 |
| LINE | 2 | 45 |

