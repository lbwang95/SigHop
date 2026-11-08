#!/usr/bin/env python3
"""Runs the jobs of one phase in parallel without letting them disturb each other.

  python3 scripts/sched.py JOBS.tsv [--max-mem-gb 900] [--max-time-jobs 12] [--force]

JOBS.tsv has one job per line: id <TAB> kind <TAB> mem_gb <TAB> threads <TAB> deps <TAB> command
  kind    gen | build | time. "time" jobs measure query times; "gen"/"build" jobs do not.
  mem_gb  reservation. The sum over running jobs never exceeds --max-mem-gb, and each job gets an
          address-space limit of 2 x mem_gb (+8 GB), so a runaway job fails with exit code 4 (MO)
          instead of triggering the system OOM killer.
  deps    comma-separated files that must exist, or "-". A job with a missing dependency is
          recorded as SKIP_DEP.

Isolation rules:
  * Every job is pinned with taskset to physical cores (one hardware thread per core).
  * build/gen jobs get whole L3 domains (CCDs on EPYC) with ceil(threads / cores-per-L3) domains.
  * time jobs are single-threaded and each gets its own L3 domain. With the user-approved
    preferred policy, memory is allocated on the local NUMA node first and may spill to other
    nodes. A time job whose reservation exceeds the largest node reserves the full scheduler
    budget and runs alone. --numa-policy bind restores the original strict binding.
    At most --max-time-jobs run at once, spread over both sockets.
  * A JOBS file must contain only one kind of job: never mix timing with building.

Logs: logs/<id>.log (stdout+stderr), logs/<id>.status (OK, TO, MO, SKIP_DEP, or FAIL:<code>),
and logs/jobs.tsv (id, status, seconds, cores). Finished jobs (OK or TO) are skipped on rerun
unless --force is given.
"""
import os, sys, time, subprocess, shutil, resource, argparse, signal, fcntl, json
from collections import defaultdict
from pathlib import Path


def node_capacities_gib():
    capacities = {}
    for path in Path('/sys/devices/system/node').glob('node[0-9]*/meminfo'):
        for line in path.read_text().splitlines():
            if 'MemTotal:' in line:
                capacities[int(path.parent.name[4:])] = int(line.split()[-2]) / (1024 ** 2)
    return capacities


def effective_reservation(kind, requested, budget, largest_node, policy):
    if kind == 'time' and policy == 'preferred' and requested > largest_node:
        return budget  # Whole budget means no other job can run at the same time.
    return requested


def topology():
    out = subprocess.run(["lscpu", "-p=CPU,CORE,SOCKET,NODE,CACHE"], capture_output=True, text=True).stdout
    first = {}
    for line in out.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cpu, core, sock, node, cache = line.split(",")
        key = (int(sock), int(core))
        if key in first:
            continue  # keep one hardware thread per physical core
        l3 = cache.split(":")[-1] if cache else sock
        first[key] = (int(cpu), int(sock), int(node) if node else 0, l3)
    doms = defaultdict(list)
    for (sock, core), (cpu, s, node, l3) in first.items():
        doms[(s, l3)].append((cpu, node))
    out = []
    for (s, l3), lst in sorted(doms.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        lst.sort()
        out.append({"socket": s, "cpus": [c for c, _ in lst], "node": lst[0][1]})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs")
    ap.add_argument("--max-mem-gb", type=float, default=900)
    ap.add_argument("--max-time-jobs", type=int, default=12)
    ap.add_argument("--max-jobs", type=int, default=10**9, help="maximum concurrent jobs; use 1 for Exp-10")
    ap.add_argument("--numa-policy", choices=["preferred", "bind"], default="preferred")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--logs", default="logs")
    a = ap.parse_args()
    os.makedirs(a.logs, exist_ok=True)
    lock = open(os.path.join(a.logs, "scheduler.lock"), "a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit("another scheduler owns this experiment directory")
    doms = topology()
    if not doms:
        sys.exit("lscpu gave no topology")
    per = len(doms[0]["cpus"])
    numactl = shutil.which("numactl")
    jobs = []
    for line in open(a.jobs):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        jid, kind, mem, thr, deps, cmd = line.split("\t", 5)
        jobs.append(dict(id=jid, kind=kind, mem=float(mem), thr=int(thr), deps=deps, cmd=cmd))
    kinds = {j["kind"] for j in jobs}
    if "time" in kinds and len(kinds) > 1:
        sys.exit("a jobs file with time jobs must not contain other kinds")
    capacities = node_capacities_gib()
    if "time" in kinds and a.numa_policy == "preferred" and (not capacities or not numactl):
        sys.exit("preferred timing policy requires numactl and readable NUMA capacities")
    largest_node = max(capacities.values()) if capacities else a.max_mem_gb
    for j in jobs:
        j["requested_mem"] = j["mem"]
        j["mem"] = effective_reservation(j["kind"], j["mem"], a.max_mem_gb, largest_node, a.numa_policy)
    if "time" in kinds:
        print("[policy] numa=%s largest_node=%.3fGiB large_jobs_reserve=%.0fGiB" %
              (a.numa_policy, largest_node, a.max_mem_gb), flush=True)

    def status_path(j):
        return os.path.join(a.logs, j["id"] + ".status")

    pending = []
    for j in jobs:
        sp = status_path(j)
        if not a.force and os.path.exists(sp) and open(sp).read().strip() in ("OK", "TO"):
            continue
        if not a.force and os.path.exists(sp) and open(sp).read().strip()=='MO' and Path('logs/final_mo',j['id']+'.json').exists():
            print('[reuse] %s final MO after900GiB retry' % j['id'],flush=True)
            continue
        pending.append(j)
    # order: time jobs alternate sockets; otherwise keep file order
    free = list(range(len(doms)))
    if "time" in kinds:  # alternate sockets so that concurrent timing jobs spread out
        rank = {}
        seen = defaultdict(int)
        for i in free:
            rank[i] = seen[doms[i]["socket"]]
            seen[doms[i]["socket"]] += 1
        free.sort(key=lambda i: (rank[i], doms[i]["socket"]))
    running = {}
    used_mem = 0.0
    summary = open(os.path.join(a.logs, "jobs.tsv"), "a")
    cancel_path=Path('results/provenance/revision_cancellations.json')
    cancelled=json.loads(cancel_path.read_text()).get('cancelled_jobs',{}) if cancel_path.exists() else {}

    def preserve(j):
        old = [os.path.join(a.logs, j["id"] + ext) for ext in (".log", ".status")]
        if any(os.path.exists(path) for path in old):
            archive = os.path.join(a.logs, "attempts", j["id"], str(time.time_ns()))
            os.makedirs(archive)
            for path in old:
                if os.path.exists(path): shutil.copy2(path, archive)

    def finish(p):
        nonlocal used_mem
        j, t0, ds, logf = running.pop(p)
        logf.close()
        code = p.returncode
        st = "OK" if code == 0 else "TO" if code == 3 else "MO" if code in (4, -9, 137) else "FAIL:%d" % code
        open(status_path(j), "w").write(st + "\n")
        summary.write("%s\t%s\t%.1f\t%s\n" % (j["id"], st, time.time() - t0, ",".join(map(str, ds))))
        summary.flush()
        free.extend(ds)
        if a.max_jobs == 1 and "time" not in kinds:
            free.sort()  # Exp-10: nested core sets, same socket, without round-robin drift.
        used_mem -= j["mem"]
        print("[done] %-50s %s (%.0fs)" % (j["id"], st, time.time() - t0), flush=True)
        return st.startswith("FAIL")

    while pending or running:
        failed = False
        for p in list(running):
            if p.poll() is not None:
                failed = finish(p) or failed
        if failed:
            # Only stop process groups created by this scheduler; keep all evidence.
            for p, (j, t0, ds, logf) in list(running.items()):
                try:
                    os.killpg(p.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            for p, (j, t0, ds, logf) in list(running.items()):
                try:
                    p.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                    p.wait()
                logf.close()
                open(status_path(j), "w").write("ABORTED_FAILURE\n")
                summary.write("%s\tABORTED_FAILURE\t%.1f\t%s\n" % (j["id"], time.time()-t0, ",".join(map(str, ds))))
            for j in pending:
                open(status_path(j), "w").write("NOT_RUN_FAILURE\n")
            summary.close()
            sys.exit("FAIL detected; stopped this phase")
        started = False
        for j in list(pending):
            if cancelled.get(j['id'],{}).get('cancelled'):
                if not os.path.exists(status_path(j)) or not open(status_path(j)).read().startswith('SKIP_PLAN'):
                    preserve(j)
                    open(status_path(j),'w').write('SKIP_PLAN removed by authoritative EXPERIMENTS_V2.md\n')
                    summary.write('%s\tSKIP_PLAN\t0.0\t-\n' % j['id']);summary.flush()
                pending.remove(j)
                print('[skip] %s (removed by updated plan)' % j['id'],flush=True)
                continue
            missing = [d for d in j["deps"].split(",") if d not in ("-", "") and not os.path.exists(d)]
            if missing:
                preserve(j)
                open(status_path(j), "w").write("SKIP_DEP " + " ".join(missing) + "\n")
                open(os.path.join(a.logs, j["id"] + ".log"), "w").write("# missing dependencies: " + " ".join(missing) + "\n")
                summary.write("%s\tSKIP_DEP\t0\t-\n" % j["id"]); summary.flush()
                pending.remove(j)
                print("[skip] %s (missing %s)" % (j["id"], missing[0]), flush=True)
                continue
            need = 1 if j["kind"] == "time" else max(1, -(-j["thr"] // per))
            if j["kind"] == "time" and sum(1 for x in running.values() if x[0]["kind"] == "time") >= a.max_time_jobs:
                break
            if len(running) >= a.max_jobs:
                break
            if need > len(free) or used_mem + j["mem"] > a.max_mem_gb:
                if not running:
                    if need > len(doms) or j["mem"] > a.max_mem_gb:
                        sys.exit("job %s cannot fit the machine" % j["id"])
                break
            ds = [free.pop(0) for _ in range(need)]
            cpus = [c for d in ds for c in doms[d]["cpus"]]
            if j["kind"] == "time":
                cpus = cpus[:1]
            else:
                cpus = cpus[:max(1, j["thr"])]
            env = dict(os.environ, OMP_NUM_THREADS=str(len(cpus) if j["kind"] != "time" else 1),
                       OMP_PROC_BIND="close", OMP_PLACES="cores")
            prefix = ["taskset", "-c", ",".join(map(str, cpus))]
            if j["kind"] == "time" and numactl:
                option = "--preferred=%d" if a.numa_policy == "preferred" else "--membind=%d"
                prefix = [numactl, option % doms[ds[0]]["node"]] + prefix
            lim = int((2 * j["mem"] + 8) * (1 << 30))

            def setlim(lim=lim):
                try:
                    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))
                except (ValueError, OSError):
                    pass  # not supported on this platform

            # Reruns must not erase the first MO/FAIL/aborted attempt.
            preserve(j)
            logf = open(os.path.join(a.logs, j["id"] + ".log"), "w")
            logf.write("# %s\n# cpus %s\n" % (j["cmd"], cpus))
            logf.write("# requested_mem_gib %.3f effective_mem_gib %.3f address_space_limit_gib %.3f\n" %
                       (j["requested_mem"], j["mem"], lim / (1 << 30)))
            logf.write("# numa_policy %s local_node %d launcher %s\n" %
                       (a.numa_policy if j["kind"] == "time" else "default", doms[ds[0]]["node"], " ".join(prefix)))
            logf.flush()
            p = subprocess.Popen(prefix + [sys.executable, "scripts/check_job.py", j["cmd"]], stdout=logf, stderr=subprocess.STDOUT, env=env,
                                 preexec_fn=setlim, start_new_session=True)
            running[p] = (j, time.time(), ds, logf)
            used_mem += j["mem"]
            pending.remove(j)
            started = True
            print("[run ] %-50s cpus=%d mem=%.0fG" % (j["id"], len(cpus), j["mem"]), flush=True)
        if not started:
            time.sleep(2)
    summary.close()


if __name__ == "__main__":
    main()
