# SigHop: Efficiently Answering Travel-Time Quantile Queries over Historical Traffic

This repository contains the source code and data for SigHop.

SigHop answers travel-time quantile queries over historical traffic using a low-rank model of edge travel times. The code also includes CCH-All, CCH-Batch, Dijkstra, and TTP-adapt.

## Usage

### Environment

- Linux, or Windows with WSL.
- A C++17 compiler with OpenMP, GNU Make, and Bash.
- Python 3.11 or later. NumPy is needed for case analysis and model preparation.

### Compilation

```bash
git clone https://github.com/lbwang95/SigHop.git
cd SigHop/code
make
bash tests/run_tests.sh
```

The tests check query answers against shortest-path searches and print `ALL TESTS PASSED` on success. Tests for the two ablations are also available:

```bash
bash tests/run_pointsonly_tests.sh
bash tests/run_ffref_tests.sh
```

### Running Example

To run all five methods on the Chengdu model, use the following command in `code/`:

```bash
bash scripts/run_example.sh
```

The script builds the index and runs 20 queries per method. It writes the results to `code/results/example/`. An output directory can also be specified: `bash scripts/run_example.sh results/example2`.

### Execution

Build the index:

```bash
mkdir -p idx results
bin/htq build --data ../data/models/CD_r4.bin --name CD_r4 \
  --prune points --threads 8 --index idx/CD_r4_points.idx --csv results/build.csv
```

Run SigHop queries for 17:30 on the weekdays in the final 30 days:

```bash
bin/htq query --data ../data/models/CD_r4.bin --name CD_r4 \
  --method sighop --index idx/CD_r4_points.idx \
  --qdir ../data/queries/CD --sets D8 --per 1000 --rounds 3 --warmup 16 \
  --phi 'cd1730:dow=0-4;days=15-45;slots=70-71' --csv results/query.csv
```

`--method` accepts `sighop`, `cch`, `cchbatch`, `dijkstra`, and `ttp`. The last four methods do not need `--index`. The query ablations are `sighop-nohub`, `sighop-countsep`, and `sighop-ffref`. Index construction accepts `--prune points`, `pointsonly`, or `box`.

`days` and `slots` use half-open ranges; `dow=0-4` means Monday through Friday. There are 96 quarter-hour slots per day, so slot 70 starts at 17:30. Quantile levels default to `0.5,0.8,0.95` and can be changed with `--alphas`.

## Data Description

The fitted traffic models and query sets are included in `data/`. The four DIMACS networks are available as a separate [download](https://github.com/lbwang95/SigHop/releases/download/data-v1/sighop-dimacs-v1.tar.gz).

| Data | Location | Description |
|---|---|---|
| Chengdu (CD), Chicago (CHI) | `data/models/` | Fitted models at ranks 2, 4, 6, and 8; `HTQDATA1` format |
| Query sets | `data/queries/` | The original source-target pairs used in the experiments |
| NY, BAY, COL, FLA | Download attachment | Original DIMACS travel-time graphs and coordinates |

Download and verify the DIMACS files from `code/`:

```bash
python3 scripts/download_dimacs.py
python3 scripts/verify_data.py
```

See [Data](docs/DATA.md) for the sources and file formats, and [Experiments](docs/EXPERIMENTS.md) for the settings and scripts used in the paper.

## Code Structure

```text
code/
  src/          SigHop, graph decomposition, baselines, and command-line drivers
  model/        Traffic-model fitting and reference queries
  scripts/      Data checks, experiment job generation, scheduling, and case analysis
  tests/        Correctness tests
data/
  models/       Fitted CD and CHI models
  queries/      Fixed query sets
docs/           Data format and experiment instructions
```
