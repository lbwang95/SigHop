#!/usr/bin/env python3
"""Summarizes results/*.csv into results/RESULTS.md (one table per experiment of PLAN.md).

  python3 scripts/collect.py [results]

Times are microseconds per query (median over rounds). TO = timed out, MO = out of memory,
NF = not run because N copies of the per-snapshot index exceed 1 TB, -- = not applicable.
"""
import csv, glob, os, sys
from collections import defaultdict

R = sys.argv[1] if len(sys.argv) > 1 else "results"
GRAPHS = ["NY", "BAY", "COL", "FLA"]
SETS = ["D%d" % i for i in range(1, 11)]


def rows(path):
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as f:
        rd = csv.reader(f)
        hdr = None
        for r in rd:
            if not r:
                continue
            if hdr is None or r == hdr:
                hdr = r
                continue
            out.append(dict(zip(hdr, r)))
    return out


def fmt(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v == 0:
        return "0"
    if abs(v) >= 1e5 or abs(v) < 1e-2:
        return "%.3g" % v
    if abs(v) >= 100:
        return "%.0f" % v
    return "%.3g" % v


def table(head, body):
    s = "| " + " | ".join(head) + " |\n| " + " | ".join("---" for _ in head) + " |\n"
    for b in body:
        s += "| " + " | ".join(str(x) for x in b) + " |\n"
    return s + "\n"


def tval(r):
    if r is None:
        return "NF"
    if r.get("status") == "TO":
        return "TO"
    v = fmt(r["mean_us"])
    if int(float(r.get("timeouts", 0) or 0)) > 0:
        v += " (%s/%s TO)" % (r["timeouts"], r["queries"])
    return v


build = rows(os.path.join(R, "build.csv"))
mem = {r["name"]: r for r in rows(os.path.join(R, "mem.csv"))}
B = {(r["name"], r["prune"]): r for r in build}
E = {e: rows(os.path.join(R, e + ".csv")) for e in ["e1", "e2", "e3", "e4", "e5", "e7", "e8", "e9", "e11"]+["x%d"%i for i in range(1,8)]}
status = {}
if os.path.exists("logs"):
    for p in glob.glob("logs/*.status"):
        status[os.path.basename(p)[:-7]] = open(p).read().strip()


def find(exp, name, method, set_="D10", phi=None):
    for r in E[exp]:
        if r["name"] == name and r["method"] == method and r["set"] == set_ and (phi is None or r["phi"] == phi):
            return r
    return None


md = ["# Results\n", __doc__.split("\n\n", 1)[1] if "\n\n" in __doc__ else "", "\n"]
if os.path.exists(os.path.join(R, "env.txt")):
    md.append("## Environment\n\n```\n" + open(os.path.join(R, "env.txt")).read() + "\n```\n\n")

# audit
bad = []
for e, rs in E.items():
    for r in rs:
        if int(float(r.get("mismatches", 0) or 0)) > 0:
            bad.append((e, r["name"], r["method"], r["set"], r["phi"], r["mismatches"], r["verified"]))
md.append("## Correctness audit\n\n")
md.append("All verified answers match.\n\n" if not bad else table(["exp", "name", "method", "set", "phi", "mismatches", "verified"], bad))
fails = sorted((k, v) for k, v in status.items() if v not in ("OK",))
if fails:
    md.append("Jobs that did not finish with OK:\n\n" + table(["job", "status"], fails))

# datasets and index cost
md.append("## Datasets and index cost (Table: index cost)\n\n")
body = []
for r in build:
    m = mem.get(r["name"], {})
    body.append([r["name"], r["prune"], r["N"], r["n"], r["m"], r["tw"], r["height_levels"], r["threads"], fmt(r["build_s"]),
                 fmt(r["label_mb"]), fmt(r["set_mean"]), r["set_max"], fmt(r["boxkept_mean"]), r["boxkept_max"],
                 fmt(r["peak_rss_mb"]), r["status"], fmt(m.get("cch_world_mb", "")), fmt(m.get("cch_all_mb", "")),
                 fmt(m.get("cch_custom_ms", "")), fmt(m.get("label2hop_all_mb", "")), fmt(m.get("model_mb", ""))])
md.append(table(["name", "prune", "N", "n", "m", "tw", "levels", "threads", "build s", "label MB", "set mean", "set max",
                 "box-kept mean", "box-kept max", "peak RSS MB", "status", "CCH MB/snap", "CCH-All MB", "custom ms",
                 "2-hop-All MB", "model MB"], body))
for (k, v) in sorted(status.items()):
    if k.startswith("build_") and v != "OK":
        md.append("- %s: %s\n" % (k, v))
md.append("\n")

MAIN = ["sighop", "cchbatch", "cch", "dijkstra", "ttp"]
md.append("## Exp-1: query time (us) vs query distance, r=4, 3 months, weekdays at 08:00\n\n")
for g in GRAPHS:
    nm = "%s_r4_d90" % g
    body = [[s] + [tval(find("e1", nm, m, s)) for m in MAIN] for s in SETS]
    k = next((r["k"] for r in E["e1"] if r["name"] == nm), "?")
    md.append("### %s (k = %s selected snapshots)\n\n" % (g, k) + table(["set"] + MAIN, body))

md.append("## Exp-2: query time (us) vs selected days (weekdays of the last w days), D10\n\n")
M2 = ["sighop", "cchbatch", "cch", "dijkstra"]
for g in GRAPHS:
    nm = "%s_r4_d90" % g
    body = []
    for w in [7, 14, 30, 60, 90]:
        src = ("e1", None) if w == 90 else ("e2", "wk%d" % w)
        rs = [find(src[0], nm, m, "D10", src[1]) for m in M2]
        k = next((r["k"] for r in rs if r), "?")
        body.append([w, k] + [tval(r) for r in rs])
    md.append("### %s\n\n" % g + table(["w (days)", "k"] + M2, body))

md.append("## Exp-3: query time (us) vs rank r, NY, D10\n\n")
for g in ["NY"]:
    body = []
    for r in [2, 4, 6, 8]:
        nm = "%s_r%d_d90" % (g, r)
        ex = "e1" if r == 4 else "e3"
        body.append([r] + [tval(find(ex, nm, m)) for m in M2])
    md.append("### %s\n\n" % g + table(["r"] + M2, body))

md.append("## Exp-4: query time (us) vs number of quantile levels, NY, D10\n\n")
body = []
for na in [1, 2, 3, 5, 10]:
    line = [na]
    for m in ["sighop", "cchbatch"]:
        src = E["e1"] if na == 3 else E["e4"]
        hit = [r for r in src if r["name"] == "NY_r4_d90" and r["method"] == m and r["set"] == "D10" and
               int(r.get("levels", 3)) == na]
        line.append(tval(hit[0]) if hit else "NF")
    body.append(line)
md.append(table(["|A|", "sighop", "cchbatch"], body))

md.append("## Exp-5: query pruning ablation, NY, r=4\n\n")
body = []
for s in SETS:
    a = find("e1", "NY_r4_d90", "sighop", s)
    b = find("e5", "NY_r4_d90", "sighop-nohub", s)
    c = find("e5", "NY_r4_d90", "sighop-countsep", s)
    body.append([s, tval(a), tval(b), tval(c), fmt(a["hubs"]) if a else "", fmt(a["hubs_kept"]) if a else "",
                 fmt(c["hubs_kept"]) if c else "", fmt(a["cands"]) if a else "", fmt(a["cands_kept"]) if a else "",
                 fmt(c["cands_kept"]) if c else ""])
md.append(table(["set", "SigHop us", "no separator pruning us", "separator by count us", "|H|", "|H| kept",
                 "|H| kept (count)", "signatures", "signatures kept", "signatures kept (count)"], body))

md.append("## Exp-7: history length, NY, r=4 (all weekdays of the history at 08:00), D10\n\n")
body = []
for d in [7, 30, 90, 180, 365, 730, 1095]:
    nm = "NY_r4_d%d" % d
    b = B.get((nm, "points"))
    m = mem.get(nm, {})
    ex = "e1" if d == 90 else "e7"
    rs = [find(ex, nm, x) for x in M2]
    k = next((r["k"] for r in rs if r), "?")
    body.append([d, b["N"] if b else m.get("N", "?"), k, fmt(b["build_s"]) if b else status.get("build_" + nm, "?"),
                 fmt(b["label_mb"]) if b else "", fmt(b["set_mean"]) if b else "", fmt(m.get("cch_all_mb", "")),
                 fmt(m.get("label2hop_all_mb", ""))] + [tval(r) for r in rs])
md.append(table(["days", "N", "k", "SigHop build s", "SigHop MB", "set mean", "CCH-All MB", "2-hop-All MB"] +
                ["%s us" % x for x in M2], body))

md.append("## Exp-8: index pruning ablation, NY, r=4 (points = box filter + point test; box = box filter only)\n\n")
body = []
for pr in ["points", "box"]:
    b = B.get(("NY_r4_d90", pr))
    body.append([pr, fmt(b["build_s"]) if b else status.get("build_NY_r4_d90" + ("_box" if pr == "box" else ""), "?"),
                 fmt(b["label_mb"]) if b else "", fmt(b["set_mean"]) if b else "", b["set_max"] if b else "",
                 fmt(b["peak_rss_mb"]) if b else ""] +
                [tval(find("e1" if pr == "points" else "e8", "NY_r4_d90", "sighop", s)) for s in ["D1", "D5", "D10"]])
md.append(table(["prune", "build s", "label MB", "set mean", "set max", "peak RSS MB", "D1 us", "D5 us", "D10 us"], body))

md.append("## Exp-9: time description, NY, D10\n\n")
body=[]
for tag in ["wk03","wk08","wk12","wk17","wk21","we08","all08"]:
    rs=[find('e1' if tag=='wk08' else 'e9','NY_r4_d90',m,'D10',None if tag=='wk08' else tag) for m in ['sighop','cchbatch','cch']]
    k=next((r['k'] for r in rs if r),'?');sig=rs[0] or {}
    body.append([tag,k]+[tval(r) for r in rs]+[fmt(sig.get(x,'')) for x in ['hubs','hubs_kept','cands_kept']])
md.append(table(['condition','k','sighop us','cchbatch us','cch us','hubs','hubs_kept','cands_kept'],body))

md.append("## Exp-10: isolated build thread sweep, NY\n\n")
thread_rows=rows(os.path.join(R,'threads.csv'));baseline=next((r for r in thread_rows if r['threads']=='1' and r['status']=='OK'),None)
body=[]
for th in [1,8,16,32,64]:
    hit=next((r for r in thread_rows if int(r['threads'])==th),None)
    if hit:
        speed=float(baseline['build_s'])/float(hit['build_s']) if baseline and hit['status']=='OK' and float(hit['build_s'])>0 else 'NA'
        body.append([th]+[fmt(hit[k]) for k in ['build_s','shortcut_s','label_s']]+[fmt(speed),fmt(hit['peak_rss_mb']),hit['status']])
    else:body.append([th,'NA','NA','NA','NA','NA',status.get('threads_NY_r4_d90_t%d'%th,'NF')])
md.append(table(['threads','build s','shortcut s','label s','speedup vs 1','peak RSS MiB','status'],body))

md.append("## Exp-11: day-level noise, NY, D10\n\n")
body=[]
for noise in [0.05,0.1,0.2,0.4]:
    nm='NY_r4_d90'+('' if noise==0.1 else '_n%g'%noise);b=B.get((nm,'points'));ex='e1' if noise==0.1 else 'e11'
    body.append([noise]+([fmt(b[x]) for x in ['build_s','label_mb','set_mean']]+[b['set_max'],b['status']] if b else ['NA']*4+[status.get('build_'+nm,'NF')])+[tval(find(ex,nm,m)) for m in ['sighop','cchbatch']])
md.append(table(['day noise','build s','label MB','set mean','set max','build status','sighop us','cchbatch us'],body))

real = sorted(glob.glob(os.path.join(R, "real", "*.csv")))
if real:
    md.append("## Real-traffic results (CD, CHI)\n\n")
    for p in real:
        rs = rows(p)
        if not rs:
            continue
        md.append("### %s\n\n" % os.path.basename(p) + table(list(rs[0].keys()), [list(r.values()) for r in rs]))

def extra_started(exp):
    return bool(E[exp]) or bool(glob.glob('jobs_'+exp+'_*.tsv'))

if extra_started('x1'):
    md.append('## X1: evening 17:00 distance curves\n\n')
    for g in GRAPHS:
        nm='%s_r4_d90'%g
        md.append('### '+g+'\n\n'+table(['set','sighop us','cchbatch us','cch us','dijkstra us'],[[s]+[tval(find('x1',nm,m,s,'wk17')) for m in ['sighop','cchbatch','cch','dijkstra']] for s in SETS]))
if extra_started('x2'):
    md.append('## X2: one-year histories and selected windows\n\n')
    for g in ['BAY','COL','FLA']:
        nm='%s_r4_d365'%g;b=B.get((nm,'points'));body=[]
        if b:md.append(table(['name','build s','label MB','set mean','status'],[[nm,fmt(b['build_s']),fmt(b['label_mb']),fmt(b['set_mean']),b['status']]]))
        for w in [7,30,90,180,365]:
            rs=[find('x2',nm,m,'D10','wk%d'%w) for m in M2]
            body.append([w,next((r['k'] for r in rs if r),'?')]+[tval(r) for r in rs])
        md.append('### '+g+'\n\n'+table(['window','k']+[m+' us' for m in M2],body))
if extra_started('x3'):
    md.append('## X3: rank on BAY (r4 reuses Exp-1)\n\n')
    body=[]
    for rank in [2,4,6,8]:
        nm='BAY_r%d_d90'%rank;b=B.get((nm,'points'))
        body.append([rank,fmt(b['build_s']) if b else 'NA',fmt(b['label_mb']) if b else 'NA']+[tval(find('e1' if rank==4 else 'x3',nm,m)) for m in M2])
    md.append(table(['r','build s','label MB']+[m+' us' for m in M2],body))
if extra_started('x4'):
    md.append('## X4: query pruning on FLA\n\n')
    md.append(table(['set','sighop us','nohub us','countsep us'],[[s]+[tval(find('x4','FLA_r4_d90',m,s)) for m in ['sighop','sighop-nohub','sighop-countsep']] for s in SETS]))
if extra_started('x5') or ('BAY_r4_d90','box') in B or 'build_BAY_r4_d90_box' in status:
    md.append('## X5: box-only index on BAY\n\n');body=[]
    for pr in ['points','box']:
        b=B.get(('BAY_r4_d90',pr));ex='e1' if pr=='points' else 'x5'
        body.append([pr,fmt(b['build_s']) if b else 'NA',fmt(b['label_mb']) if b else 'NA',b['status'] if b else status.get('build_BAY_r4_d90_box','NF')]+[tval(find(ex,'BAY_r4_d90','sighop',s)) for s in SETS])
    md.append(table(['prune','build s','label MB','status']+SETS,body))
if extra_started('x6'):
    md.append('## X6: TTP-adapt parameter, NY\n\n')
    md.append(table(['ttp-k','D1 us','D2 us','D3 us','D4 us'],[[k]+[tval(find('x6','NY_r4_d90_ttpk%d'%k,'ttp','D%d'%i)) for i in range(1,5)] for k in [1,3,5]]))
if extra_started('x7'):
    md.append('## X7: r8 distance curves on NY\n\n')
    md.append(table(['set','sighop us','cchbatch us','cch us'],[[s]+[tval(find('x7','NY_r8_d90',m,s)) for m in ['sighop','cchbatch','cch']] for s in SETS]))

open(os.path.join(R, "RESULTS.md"), "w").write("".join(md))
print("wrote", os.path.join(R, "RESULTS.md"))
