# Data sources and format

## Datasets

`data/manifest.json` lists the file sizes, checksums, model dimensions, and query counts. Run `python3 scripts/verify_data.py` from `code/` to check the files.

CD and CHI include fitted models at ranks 2, 4, 6, and 8. Each model stores edge loadings and snapshot factors, with travel times measured in seconds.

| Dataset | Source | Period | Notes |
|---|---|---|---|
| CD | [Chengdu traffic data, Figshare v4](https://doi.org/10.6084/m9.figshare.7140209.v4) | 2015-06-01 to 2015-07-15 | The source is CC0. Its published two-minute speeds were aggregated to complete quarter-hours. |
| CHI | [Chicago Traffic Tracker historical congestion estimates](https://data.cityofchicago.org/Transportation/Chicago-Traffic-Tracker-Historical-Congestion-Estima/sxs8-h27x) | 2019-01-01 to 2019-02-14 | Published congestion estimates; see the [Chicago Data Portal terms](https://www.chicago.gov/city/en/narr/foia/data_disclaimer.html). |
| NY, BAY, COL, FLA | [9th DIMACS Implementation Challenge](http://www.diag.uniroma1.it/challenge9/download.shtml) | Static road networks | Travel-time files use deciseconds. |

Preprocessing details are recorded in `data/sources.json`. CD has 1,780 snapshots and CHI has 3,268, both spanning 45 calendar days. Snapshots with no positive speed observations are removed after fitting. All ranks use the same snapshots.

## DIMACS download

The [download archive](https://github.com/lbwang95/SigHop/releases/download/data-v1/sighop-dimacs-v1.tar.gz) contains a travel-time graph (`.gr`) and a coordinate file (`.co`) for each of NY, BAY, COL, and FLA. Checksums are listed in `data/dimacs_manifest.json`.

```bash
cd code
python3 scripts/download_dimacs.py
# Or use a previously downloaded copy:
python3 scripts/download_dimacs.py --archive /path/to/sighop-dimacs-v1.tar.gz
```

The script checks the archive and installs the files under `data/dimacs/`. NY, BAY, and COL are also available in [RSP](https://github.com/lbwang95/RSP/tree/main/data).

Generate one semi-synthetic input using the bundled CD factors:

```bash
mkdir -p ../data/generated
bin/gen dimacs --gr ../data/dimacs/USA-road-t.NY.gr \
  --co ../data/dimacs/USA-road-d.NY.co --unit 0.1 \
  --r 4 --z-from ../data/models/CD_r4.bin --days 90 --seed 1 \
  --out ../data/generated/NY_r4_d90.bin
```

The semi-synthetic traffic is generated from DIMACS free-flow times, spatially varying nonnegative loadings, and day-level sampling of CD factors.

## Query files

Each line contains `source target`, using zero-based vertex IDs of the corresponding model. NY/BAY/COL/FLA have ten distance groups, D1 through D10, with 1,000 pairs per group. Counts for CD and CHI are listed in the manifest. Their query sets were generated using the rank-8 model and reused across ranks.

`data/queries/CD/MAP.txt` contains 1,901 pairs `(1291, target)`, one for every other reachable vertex. The MAP benchmark evaluates each pair as a separate point-to-point query.

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

The fitting and reference-query programs are in `code/model/`:

```bash
cd code
make model-tools
python3 -m pip install -r requirements.txt
```

`prepare_chengdu.py` reads the source catalog, link table, and speed files produced by `pack_chengdu.py` (which requires PyArrow). `aggregate_fixed_slots.py` aggregates them into 15-minute slots using overlap-weighted pace. `prepare_chicago.py` reads daily speed and geometry CSVs with a `download_manifest.jsonl` file and aggregates the irregular reports directly.

For a prepared city grid containing `meta.txt`, `arcs.bin`, `coords.bin`, `days.bin`, `slots.bin`, `obs.bin`, and the holdout masks:

```bash
mkdir -p ../data/work/fitted
OMP_NUM_THREADS=4 bin/nmf --in /path/to/prepared-grid --out ../data/work/fitted \
  --ranks 2,4,6,8 --iters 60 --seed 1 \
  --holdout-mask /path/to/prepared-grid/holdout.bin --min-speed 0.000001 --name CITY
```

The holdout unit is an edge-day block, with seed 20261105 and probability 0.05. Free-flow speeds are estimated from training observations only. Copy the grid's `preparation.json` and `edge_mapping.csv` into the fitted directory, then use `model/drop_empty_worlds.py --input ... --out ...` to remove snapshots with no observations. The models in `data/models/` are those used for the reported experiments.
