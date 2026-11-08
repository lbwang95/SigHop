# Experiments

All commands below run from `code/`. Recorded CSVs are under the repository's top-level `results/`; new executions write to `code/results/`. They are kept separate so a replay does not overwrite the recorded results.

## Small example and correctness

```bash
make && bash tests/run_tests.sh
bash tests/run_pointsonly_tests.sh
bash tests/run_ffref_tests.sh
bash scripts/run_example.sh
```

The first test ends with `ALL TESTS PASSED`. The example uses all five methods on the bundled CD model. Query CSVs record `completed`, `timeouts`, `verified`, and `mismatches`. A reported nonzero mismatch stops the scheduler. Verification uses the original implementation's floating-point tolerance; CCH-All uses a smaller direct-Dijkstra check sample, recorded in `verified`.

## Main road-network experiments

Download the DIMACS attachment and use the included frozen query sets:

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

`gen` creates traffic from the frozen CD factors. Query generation is omitted above because the original queries are already included. `make_jobs.py genq` is available for a new dataset; do not run it over the frozen published query directories.

The full matrix needs substantial RAM and disk space. The original server had 192 physical cores; the scheduler reserved at most 900 GiB in total. Timing tasks use one physical core and one private L3 cache domain each, with at most 12 tasks at once. The NUMA policy prefers local memory and permits allocation on other nodes. Large jobs reserve the entire scheduling budget. There is no uniform 64 GiB process limit. Use `GRAPHS="NY"`, then another graph, when indexes do not fit together. Increasing or reducing concurrency changes the timing conditions and should be reported.

The matrix uses rank 4, 90 days, 08:00 weekdays, and quantiles 0.5/0.8/0.95 as defaults. Rank, history, quantile count, time description, noise, and ablation settings are defined in `make_jobs.py`. Dijkstra and TTP use smaller query counts than the indexed methods, as recorded in each job and CSV. TTP uses a 60-second per-query timeout and stops early if its first five queries all time out.

Thread-count builds must run alone:

```bash
python3 scripts/make_jobs.py threads > jobs_threads.tsv
python3 scripts/sched.py jobs_threads.tsv --max-jobs 1 --max-mem-gb 900
```

Additional matrices are available through `--extra X1` through `--extra X7`. The recorded V2 results include X1-X4. X5 was stopped; X7 was not run. The V3 `x6.csv` contains the later TTP parameter checks. Missing measurements are not filled with synthetic numbers.

## V3 and V4 additions

The portable adapters preserve the original additional task definitions. V3 timing jobs are generated with `python3 scripts/make_jobs.py --v3 time`; their main addition is TTP for the selected-day sweep. Its point-only query task is generated only when the matching local build result exists.

For V4, first have the main generated datasets and points indexes available. Run the default query regression, then the six construction stages strictly one at a time. The seventh stage builds the small CD index if needed:

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

No other construction, timing, or case-analysis task should run during these V4 builds. Temporary construction indexes are named `idx/tmp_v4_*`; they can be removed after recording each build. Keep the FLA box index until its D1/D10 query checks have completed. The retained CD index is needed for MAP timing. Build timeout flags retain the original stage-level timing semantics.

The NY default regression should match the recorded D1/D10 `hubs_kept` and `cands_kept` values (1 / 1.032 and 12.627 / 139.169). Parallel construction may produce the small recorded `boxkept_mean` ordering differences. Do not select the best run or replace failed or timed-out rows.

## Chengdu case analysis

Install NumPy and run `scripts/case_v4.py` through the case task above. It first checks the 40-slot source-target profile against `results/v2/real/case_cd.csv`, then generates 20,000 pair-slot rows and 3,802 map rows. The larger pair/map CSVs are regenerated from the bundled model and queries instead of being stored in Git.

The profile is for 1291 to 1127 over 22 weekdays in days [15,45). Shortest-path ties use sequential float64 component sums, a heap ordered by (distance, vertex ID), original outgoing edge order, and predecessor updates only for strictly shorter paths. Routes are distinguished by their directed edge sequences.

`mean_eta` is the shortest travel time under the mean factor vector. `late_mean` counts days whose daily optimum exceeds that estimate; it is a lower bound on late days when following the fixed mean-condition route. `fixed_p50` and `fixed_p95` evaluate that fixed route across the snapshots. MAP timings run 1,901 independent source-target queries for each method and do not measure an optimized single-source algorithm.

## Recorded results

| Directory | Contents |
|---|---|
| `results/v2/` | Main distance/parameter experiments, index sizes/build costs, X1-X4, and real-data model/reference summaries |
| `results/v3/` | TTP day-count/parameter additions and the first point-only experiment |
| `results/v4/` | Completed parameter baselines, expanded query ablations, isolated construction results, and CD case summaries |

The V4 run completed 56 scheduled tasks, with 93,004 answer checks and zero mismatches. It retains 30 TTP single-query timeouts. These totals include the default query regression. An `OK` process or CSV status can still contain individual TTP timeouts: always inspect `completed` and `timeouts` when interpreting its latency. Model/reference errors in the real-data results are relative to published speed estimates with the documented common fill, not physical trip truth.
