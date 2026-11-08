# SigHop: Efficiently Answering Travel-Time Quantile Queries over Historical Traffic

This repository contains the implementation, fitted traffic models, query sets, and main experimental results for SigHop.

Given a source, a target, a time description, and quantile levels, SigHop returns quantiles of the shortest travel times over the selected historical snapshots. Edge travel times are represented by nonnegative edge loadings and snapshot factors. The implementation includes SigHop, CCH-All, CCH-Batch, Dijkstra, and TTP-adapt.

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

The final test line should be `ALL TESTS PASSED`. Additional tests cover the point-only construction and the free-flow reference query variant:

```bash
bash tests/run_pointsonly_tests.sh
bash tests/run_ffref_tests.sh
```

### Running Example

The repository includes the small Chengdu model and its original query sets. From `code/`, run:

```bash
bash scripts/run_example.sh
```

This builds an index and runs 20 queries with each of the five methods. Results are written to `code/results/example/`. Use a different output directory for another run, for example `bash scripts/run_example.sh results/example2`.

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

The Git repository contains only small fitted models, fixed query sets, and summary results. The four large DIMACS inputs are a separate [download attachment](https://github.com/lbwang95/SigHop/releases/download/data-v1/sighop-dimacs-v1.tar.gz). Generated traffic and indexes are rebuilt locally.

| Data | Included or downloaded | Description |
|---|---|---|
| Chengdu (CD), Chicago (CHI) | `data/models/` | Fitted models at ranks 2, 4, 6, and 8; `HTQDATA1` format |
| Query sets | `data/queries/` | The original source-target pairs used in the experiments |
| NY, BAY, COL, FLA | Download attachment | Original DIMACS travel-time graphs and coordinates |
| Recorded results | `results/v2/`, `results/v3/`, `results/v4/` | Main numerical results, including timeout counts |

Download and verify the DIMACS files from `code/`:

```bash
python3 scripts/download_dimacs.py
python3 scripts/verify_data.py
```

The downloader checks both the archive and the individual files. [Data sources and format](docs/DATA.md) describe provenance, units, licenses, and model preparation. [Experiment instructions](docs/EXPERIMENTS.md) explain generation, scheduling, and the result files.

## Code Structure

```text
code/
  src/          SigHop, graph decomposition, baselines, and command-line drivers
  model/        Original traffic-model fitting and reference-analysis tools
  scripts/      Data checks, experiment job generation, scheduling, and case analysis
  tests/        Correctness tests
data/
  models/       Small fitted CD and CHI models
  queries/      Fixed query sets
results/        Recorded numerical results
docs/           Data format, experimental settings, and release notes
```

The six files in `code/src/` are unchanged from the final V4 experiments; their hashes are recorded in [source_hashes.json](docs/source_hashes.json). The public copy adapts paths and removes machine-specific controllers. It does not change the algorithms or the experiment parameter matrix.
