# Experiments

Run the commands below from `code/`. Output files are written to `code/results/`.

## Tests

```bash
make && bash tests/run_tests.sh
bash tests/run_pointsonly_tests.sh
bash tests/run_ffref_tests.sh
bash scripts/run_example.sh
```

The tests compare query answers against shortest-path searches. The example runs all five methods on CD. The `verified` column records the number of queries checked; CCH-All uses a smaller sample checked directly with Dijkstra. Query verification uses a relative tolerance of 1e-3, with an absolute tolerance of 1e-3 near zero.

## Main road-network experiments

Download the DIMACS networks and use the query sets in `data/queries/`:

```bash
python3 scripts/download_dimacs.py
python3 scripts/verify_data.py
mkdir -p idx results logs
export GRAPHS="NY BAY COL FLA"
export BUILD_THREADS=32
for phase in gen build mem time; do
  python3 scripts/make_jobs.py "$phase" > "jobs_${phase}.tsv"
  python3 scripts/sched.py "jobs_${phase}.tsv" \
    --max-mem-gb 900 --max-time-jobs 12 --numa-policy preferred
done
python3 scripts/collect.py results
```

`gen` creates traffic from the CD factors. The query sets are already included, so query generation is omitted above. For a new dataset, `make_jobs.py genq` generates query sets in `data/queries/`.

The experiments used a server with 192 physical cores and a total scheduling budget of 900 GiB. Each timing task uses one physical core and a separate L3 cache domain, with at most 12 tasks running together. Memory is allocated on the local NUMA node first and may extend to other nodes; tasks whose reservation exceeds a single node run alone. Set `GRAPHS="NY"` to run one network at a time when disk space is limited.

The matrix uses rank 4, 90 days, 08:00 weekdays, and quantiles 0.5/0.8/0.95 as defaults. Rank, history, quantile count, time description, noise, and ablation settings are defined in `make_jobs.py`. Dijkstra and TTP use smaller query counts than the indexed methods, as recorded in each job and CSV. TTP uses a 60-second per-query timeout and stops early if its first five queries all time out.

Thread-count builds must run alone:

```bash
python3 scripts/make_jobs.py threads > jobs_threads.tsv
python3 scripts/sched.py jobs_threads.tsv --max-jobs 1 --max-mem-gb 900
```

Additional experiment jobs are available through `--extra X1` through `--extra X7`. Their parameters are defined in `scripts/make_jobs.py`.

## Additional baselines and ablations

`python3 scripts/make_jobs.py --v3 time` generates TTP jobs for the selected-day experiment and the TTP parameter sweep. It also generates point-only queries when the corresponding build result is available.

The `--v4` tasks cover the remaining parameter baselines, query ablations, serial construction experiments, and CD case study. They use the datasets, points indexes, and local `results/build.csv` and `results/mem.csv` produced by the main stages above. Run the query check first, then the six construction jobs in order. The seventh job builds the CD index if needed:

```bash
python3 scripts/make_jobs.py --v4 default-check > jobs_v4_check.tsv
python3 scripts/sched.py jobs_v4_check.tsv --max-mem-gb 900 --max-time-jobs 12
for i in 1 2 3 4 5 6 7; do
  python3 scripts/make_jobs.py --v4 "build-$i" > "jobs_v4_build_$i.tsv"
  python3 scripts/sched.py "jobs_v4_build_$i.tsv" --max-jobs 1 --max-mem-gb 900
done
python3 scripts/make_jobs.py --v4 case > jobs_v4_case.tsv
python3 scripts/sched.py jobs_v4_case.tsv --max-jobs 1 --max-mem-gb 900
python3 scripts/make_jobs.py --v4 time > jobs_v4_time.tsv
python3 scripts/sched.py jobs_v4_time.tsv --max-mem-gb 900 --max-time-jobs 12 --numa-policy preferred
```

Run these builds without other experiment tasks. Temporary indexes are named `idx/tmp_v4_*` and can be removed after each build, except for the FLA box index, which is needed for the D1/D10 query checks. Keep the CD index for MAP timing. Build timeouts apply separately to the shortcut and label stages.

## Chengdu case analysis

Install NumPy and run `scripts/case_v4.py` through the case task above. The task builds `bin/model_case` to calculate a reference profile from the model, checks the Python profile against it, and then computes the pair and map outputs. No precomputed result files are required. To invoke the Python script directly, first run `make bin/model_case`.

The profile is for 1291 to 1127 over 22 weekdays in days [15,45). Shortest-path ties use sequential float64 component sums, a heap ordered by (distance, vertex ID), original outgoing edge order, and predecessor updates only for strictly shorter paths. Routes are distinguished by their directed edge sequences.

`mean_eta` is the shortest travel time under the mean factor vector. `late_mean` counts days whose daily optimum exceeds that estimate; it is a lower bound on late days when following the fixed mean-condition route. `fixed_p50` and `fixed_p95` evaluate that fixed route across the snapshots. MAP timings run 1,901 independent source-target queries for each method.
