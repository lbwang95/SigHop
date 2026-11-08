# Data sources and format

## Included inputs

`data/manifest.json` records the size, SHA-256, and header or query count of each bundled input. Run `python3 scripts/verify_data.py` from `code/` to verify them. These are copies of the frozen experimental inputs.

CD and CHI include ranks 2, 4, 6, and 8. Rank 6 was fitted separately with the same pipeline. The models store edge loadings and snapshot factors, not raw vehicle trajectories. Their weights are measured in seconds.

| Dataset | Source | Period | Conditions |
|---|---|---|---|
| CD | [Chengdu traffic data, Figshare v4](https://doi.org/10.6084/m9.figshare.7140209.v4) | 2015-06-01 to 2015-07-15 | The source is CC0. Its published two-minute speeds were aggregated to complete quarter-hours. |
| CHI | [Chicago Traffic Tracker historical congestion estimates](https://data.cityofchicago.org/Transportation/Chicago-Traffic-Tracker-Historical-Congestion-Estima/sxs8-h27x) | 2019-01-01 to 2019-02-14 | Source use is governed by the [Chicago Data Portal terms](https://www.chicago.gov/city/en/narr/foia/data_disclaimer.html). These are published congestion estimates, not measured trip ground truth. |
| NY, BAY, COL, FLA | [9th DIMACS Implementation Challenge](http://www.diag.uniroma1.it/challenge9/download.shtml) | Static road networks | Preserve the source attribution and applicable source terms. Travel-time files use deciseconds. |

The supplied real-model provenance is also recorded in `data/sources.json`. The CD model retains 1,780 snapshots; CHI retains 3,268. Both have 45 calendar days. Only snapshots with no positive published observations are removed after fitting. The same retained snapshots are used at every rank.

## DIMACS download

The [data-v1 attachment](https://github.com/lbwang95/SigHop/releases/download/data-v1/sighop-dimacs-v1.tar.gz) contains the eight original `.gr` and `.co` inputs used in the experiments. They are stored outside Git to keep cloning small. `data/dimacs_manifest.json` contains every checksum. The official direct links returned 404 during preparation, so the verified original files are archived here. No modified-weight replacement is used.

```bash
cd code
python3 scripts/download_dimacs.py
# Or use a previously downloaded copy:
python3 scripts/download_dimacs.py --archive /path/to/sighop-dimacs-v1.tar.gz
```

The files are installed under `data/dimacs/`, which is ignored by Git. NY, BAY, and COL are also available in the author's [RSP data directory](https://github.com/lbwang95/RSP/tree/main/data), but checksum verification is needed before substituting a mirror.

Generate one semi-synthetic input using the bundled CD factors:

```bash
mkdir -p ../data/generated
bin/gen dimacs --gr ../data/dimacs/USA-road-t.NY.gr \
  --co ../data/dimacs/USA-road-d.NY.co --unit 0.1 \
  --r 4 --z-from ../data/models/CD_r4.bin --days 90 --seed 1 \
  --out ../data/generated/NY_r4_d90.bin
```

This generator uses DIMACS free-flow times, spatially varying nonnegative loadings, and day-level sampling of CD factors. These are semi-synthetic traffic snapshots, not observed traffic on the four DIMACS networks.

## Query files

Each text line is `source target`, using zero-based vertex IDs of the corresponding processed model. NY/BAY/COL/FLA have ten distance groups, D1 through D10, with 1,000 pairs per group. The real-data banks retain the original files and counts; their exact counts are in the manifest. They were created using the rank-8 model and then reused across ranks. Do not replace them with newly generated queries when comparing the recorded experiments.

`data/queries/CD/MAP.txt` contains 1,901 pairs `(1291, target)`, one for every other reachable vertex. The MAP benchmark runs independent point-to-point queries; it is not an optimized single-source computation.

## HTQDATA1 binary format

All files use little-endian numbers. Arrays are stored consecutively, without padding.

| Field | Type and length |
|---|---|
| Magic | 8 bytes: `HTQDATA1` |
| Header | Seven int64 values: `n, m, r, N, days, slots_per_day, has_coordinates` |
| Slot duration | One float64, in seconds |
| Arc sources, arc targets | Two int32 arrays of length `m` |
| Edge loadings U | `m * r` float32 values, edge-major |
| Snapshot factors Z | `N * r` float32 values, snapshot-major |
| Day, slot, weekday tags | Three int32 arrays of length `N` |
| Optional x, y coordinates | Two float32 arrays of length `n` |

The travel time of edge `e` in snapshot `i` is the dot product of `U[e]` and `Z[i]`. Coordinate zero of every factor is 1. Weekday 0 is Monday. The implementation is in `code/src/common.hpp`.

## Model fitting and reference tools

The portable repository includes the original fitting and reference code under `code/model/`, separated from the current query implementation:

```bash
cd code
make model-tools
python3 -m pip install -r requirements.txt
```

`prepare_chengdu.py` consumes the source catalog, link table, and losslessly packed speed files produced by `pack_chengdu.py` (the downloader additionally needs PyArrow). `aggregate_fixed_slots.py` performs overlap-weighted pace aggregation into complete 15-minute slots. `prepare_chicago.py` consumes checked daily speed and geometry CSVs and their `download_manifest.jsonl`; it handles the irregular reports directly. Provider downloads remain separate from the small repository.

For a prepared city grid containing `meta.txt`, `arcs.bin`, `coords.bin`, `days.bin`, `slots.bin`, `obs.bin`, and the frozen holdout masks:

```bash
mkdir -p ../data/work/fitted
OMP_NUM_THREADS=4 bin/nmf --in /path/to/prepared-grid --out ../data/work/fitted \
  --ranks 2,4,6,8 --iters 60 --seed 1 \
  --holdout-mask /path/to/prepared-grid/holdout.bin --min-speed 0.000001 --name CITY
```

The holdout unit is an edge-day block, with seed 20261105 and probability 0.05. Free-flow speeds are estimated from training observations only. Copy the grid's `preparation.json` and `edge_mapping.csv` into the fitted directory, then use `model/drop_empty_worlds.py --input ... --out ...` to remove wholly unobserved worlds while preserving fitted values. Exact bundled model hashes, rather than a fresh fit on a different machine, define the released experimental inputs.
