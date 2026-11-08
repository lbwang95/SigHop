# Publication checks

The separate public copy was checked on Linux with GCC 13.3.0 on 2026-11-08.

| Check | Result |
|---|---|
| Six current C++ source files | Match the frozen V4 SHA-256 values |
| Eight fitted models and 61 query files | All 69 bundled-input checksums match |
| `make` and original correctness tests | `ALL TESTS PASSED` |
| Point-only construction tests | `ALL POINTSONLY TESTS PASSED` |
| Free-flow reference tests | `ALL FFREF TESTS PASSED` |
| `make model-tools` | All four model/reference programs compile |
| README example | Five methods complete 20 queries each; 90 answer checks, zero mismatches |
| DIMACS archive extraction | Archive and all eight input checksums match |
| Public case-analysis script | 40 profile, 20,000 pair, and 3,802 map rows are byte-identical to the original V4 outputs |
| Recorded CSVs | 35 summary files; no nonzero reported mismatches |

The example is a functionality check, not a replacement benchmark. The recorded V2-V4 timing results are preserved separately. The large road-network experiments and model fitting were not repeated for publication.
