#!/usr/bin/env python3
"""Writes the jobs of one phase of PLAN.md as a TSV for scripts/sched.py.

  python3 scripts/make_jobs.py gen|genq|build|mem|time|threads > jobs_<phase>.tsv

Phases must run in this order; "time" reads results/build.csv and results/mem.csv to size the
memory reservations and to skip CCH-All/CCH-Batch where N copies would not fit in 1 TB.
Environment: BUILD_THREADS (default 32), MEM_BUDGET_MB (default 1e6, i.e. 1 TB).
"""
import csv, os, sys
if '--v4' in sys.argv:
    from v4_jobs import main
    main()
    raise SystemExit(0)
if '--v3' in sys.argv:
    from v3_jobs import main
    main()
    raise SystemExit(0)

GRAPHS = os.environ.get("GRAPHS", "NY BAY COL FLA").split()  # restrict, e.g. GRAPHS="FLA", to work graph by graph
RANKS = [2, 4, 6, 8]
R0, D0 = 4, 90                      # default rank and history length (days)
HIST = [7, 30, 90, 180, 365, 730, 1095]  # Exp-7, NY only
WINDOWS = [7, 14, 30, 60, 90]       # Exp-2: weekdays of the last w days
SLOT = 32                           # 08:00-08:15 in the 96-slot numbering of the CD source
SETS = ",".join("D%d" % i for i in range(1, 11))
SLOTS_DAY = {"03": 12, "08": 32, "12": 48, "17": 68, "21": 84}  # Exp-9: the five windows of the CD source
DOWS = {"wk": "0-4", "we": "5-6", "all": "0-6"}             # Exp-9: weekdays, weekends, all days
NOISE = [0.05, 0.1, 0.2, 0.4]                               # Exp-11: GEN_DAY_NOISE of the generator
N0 = 0.1
THREADS_SWEEP = [1, 8, 16, 32, 64]                          # Exp-10: build threads (NY)
ALPHAS = {1: "0.5", 2: "0.5,0.95", 3: "0.5,0.8,0.95", 5: "0.1,0.25,0.5,0.75,0.9",
          10: "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95"}
REAL = ["CD", "CHI"]
BT = int(os.environ.get("BUILD_THREADS", "32"))
BUDGET_MB = float(os.environ.get("MEM_BUDGET_MB", "1e6"))
# rough build reservations in GB for r=4, d=90; scaled by rank and history below (adjust if needed)
BASE_GB = {"NY": 40, "BAY": 60, "COL": 150, "FLA": 300, "CD": 4, "CHI": 2}

H = "bin/htq"
J = []


def job(jid, kind, mem, thr, deps, cmd):
    J.append("\t".join([jid, kind, "%.0f" % mem, str(thr), ",".join(deps) if deps else "-", cmd]))


def ranks_of(g):  # the rank sweep (Exp-3, Exp-6) runs on NY only, as the parameter studies of PCSP
    return RANKS if g == "NY" else [R0]


def data(g, r, d, sfx=""):
    return "../data/generated/%s_r%d_d%d%s.bin" % (g, r, d, sfx)


def real(g, r):
    return "../data/models/%s_r%d.bin" % (g, r)


def idx(g, r, d, prune="points", sfx=""):
    return "idx/%s_r%d_d%d%s_%s.idx" % (g, r, d, sfx, prune)


def name(g, r, d, sfx=""):
    return "%s_r%d_d%d%s" % (g, r, d, sfx)


def nsfx(nz):  # data variant with another day-level noise
    return "" if nz == N0 else "_n%g" % nz


def phi(tag, lo, hi, slot=SLOT, dow="0-4"):
    return "'%s:dow=%s;slots=%d-%d;days=%d-%d'" % (tag, dow, slot, slot + 1, lo, hi)


def build_gb(g, r, d):
    return min(900, BASE_GB[g] * {2: 0.5, 4: 1, 6: 2, 8: 3}[r] * max(1.0, (d / 90.0) ** 0.5))


phase = sys.argv[1] if len(sys.argv) > 1 else ""
os.makedirs("results", exist_ok=True)
EXTRA = sys.argv[sys.argv.index('--extra')+1].upper() if '--extra' in sys.argv else None
if EXTRA and EXTRA not in ['X%d' % i for i in range(1,8)]:sys.exit('unknown extra experiment')
PY_DATA = os.environ.get('HTQ_DATA_PYTHON',sys.executable)

if '--prepare' in sys.argv or '--prepare-import' in sys.argv or '--real' in sys.argv:
    sys.exit('See ../docs/DATA.md and ../docs/EXPERIMENTS.md for portable data preparation.')

if EXTRA and phase != 'time':
    if phase not in ('gen','genq','build','mem'):sys.exit('unsupported extra phase')
    configs=[]
    if EXTRA=='X2':configs=[(g,R0,365,'points',21600) for g in GRAPHS if g in ('BAY','COL','FLA')]
    elif EXTRA=='X3':configs=[('BAY',r,D0,'points',14400) for r in (2,6,8) if 'BAY' in GRAPHS]
    elif EXTRA=='X5' and 'BAY' in GRAPHS:configs=[('BAY',R0,D0,'box',7200)]
    if '--validate' in sys.argv:
        if EXTRA!='X2' or phase!='gen':sys.exit('prefix validation belongs to X2 gen')
        gs=[g for g,_,_,_,_ in configs]
        if gs:
            job('x2_prefix_audit','gen',16,1,[data(g,R0,d) for g in gs for d in (90,365)],
                PY_DATA+' scripts/validate_history_prefix.py '+' '.join(gs))
    elif phase=='gen' and EXTRA in ('X2','X3'):
        for g,r,d,pr,to in configs:
            gr,co='../data/dimacs/USA-road-t.%s.gr'%g,'../data/dimacs/USA-road-d.%s.co'%g
            job('gen_'+name(g,r,d),'gen',16,1,[gr,co,real('CD',r)],
                'mkdir -p ../data/generated && bin/gen dimacs --gr %s --co %s --r %d --z-from %s --days %d --unit 0.1 --seed 1 --out %s'
                % (gr,co,r,real('CD',r),d,data(g,r,d)))
    elif phase=='build':
        for g,r,d,pr,to in configs:
            job('build_'+name(g,r,d)+('_box' if pr=='box' else ''),'build',build_gb(g,r,d)*(2 if pr=='box' else 1),BT,[data(g,r,d)],
                '%s build --data %s --name %s --prune %s --index %s --timeout %d --threads %d --csv results/build.csv'
                % (H,data(g,r,d),name(g,r,d),pr,idx(g,r,d,pr),to,BT))
    elif phase=='mem' and EXTRA in ('X2','X3'):
        for g,r,d,pr,to in configs:
            job('mem_'+name(g,r,d),'gen',24,1,[data(g,r,d)],
                '%s mem --data %s --name %s --csv results/mem.csv' % (H,data(g,r,d),name(g,r,d)))
    print('\n'.join(J));sys.exit(0)

if phase == "gen":
    for g in GRAPHS:
        gr, co = "../data/dimacs/USA-road-t.%s.gr" % g, "../data/dimacs/USA-road-d.%s.co" % g
        for r in ranks_of(g):
            job("gen_" + name(g, r, D0), "gen", 16, 1, [gr, co, real("CD", r)],
                "mkdir -p ../data/generated && bin/gen dimacs --gr %s --co %s --r %d --z-from %s --days %d --unit 0.1 --seed 1 --out %s"
                % (gr, co, r, real("CD", r), D0, data(g, r, D0)))
        if g == "NY":
            for nz in NOISE:
                if nz != N0:
                    job("gen_" + name(g, R0, D0, nsfx(nz)), "gen", 16, 1, [gr, co, real("CD", R0)],
                        "mkdir -p ../data/generated && GEN_DAY_NOISE=%g bin/gen dimacs --gr %s --co %s --r %d --z-from %s --days %d "
                        "--unit 0.1 --seed 1 --out %s" % (nz, gr, co, R0, real("CD", R0), D0, data(g, R0, D0, nsfx(nz))))
            for d in HIST:
                if d != D0:
                    job("gen_" + name(g, R0, d), "gen", 16, 1, [gr, co, real("CD", R0)],
                        "mkdir -p ../data/generated && bin/gen dimacs --gr %s --co %s --r %d --z-from %s --days %d --unit 0.1 --seed 1 --out %s"
                        % (gr, co, R0, real("CD", R0), d, data(g, R0, d)))

elif phase == "genq":
    for g in GRAPHS:
        job("genq_" + g, "gen", 8, 1, [data(g, R0, D0)],
            "mkdir -p ../data/queries/%s && %s genq --data %s --out ../data/queries/%s --per 1000 --sets 10 --lmin 60 --seed 7"
            % (g, H, data(g, R0, D0), g))

elif phase == "build":
    os.makedirs("idx", exist_ok=True)
    for g in GRAPHS:
        for r in ranks_of(g):
            job("build_" + name(g, r, D0), "build", build_gb(g, r, D0), BT, [data(g, r, D0)],
                "%s build --data %s --name %s --prune points --index %s --timeout 14400 --threads %d --csv results/build.csv"
                % (H, data(g, r, D0), name(g, r, D0), idx(g, r, D0), BT))
    g = "NY"
    if g not in GRAPHS:
        print("\n".join(J)); sys.exit(0)
    job("build_%s_box" % name(g, R0, D0), "build", build_gb(g, R0, D0) * 2, BT, [data(g, R0, D0)],
        "%s build --data %s --name %s --prune box --index %s --timeout 7200 --threads %d --csv results/build.csv"
        % (H, data(g, R0, D0), name(g, R0, D0), idx(g, R0, D0, "box"), BT))
    for nz in NOISE:
        if nz != N0:
            s = nsfx(nz)
            job("build_" + name(g, R0, D0, s), "build", build_gb(g, R0, D0) * 2, BT, [data(g, R0, D0, s)],
                "%s build --data %s --name %s --prune points --index %s --timeout 14400 --threads %d --csv results/build.csv"
                % (H, data(g, R0, D0, s), name(g, R0, D0, s), idx(g, R0, D0, "points", s), BT))
    for d in HIST:
        if d != D0:
            job("build_" + name(g, R0, d), "build", build_gb(g, R0, d), BT, [data(g, R0, d)],
                "%s build --data %s --name %s --prune points --index %s --timeout 21600 --threads %d --csv results/build.csv"
                % (H, data(g, R0, d), name(g, R0, d), idx(g, R0, d), BT))
    for g in (REAL if os.environ.get("REAL", "1") == "1" else []):
        for r in RANKS:
            job("build_%s_r%d" % (g, r), "build", BASE_GB[g] * 2, 8, [real(g, r)],
                "%s build --data %s --name %s_r%d --prune points --timeout 7200 --threads 8 --csv results/build.csv"
                % (H, real(g, r), g, r))

elif phase == "threads":  # Exp-10, run with sched.py --max-jobs 1 so that the builds do not overlap
    if "NY" in GRAPHS:
        for th in THREADS_SWEEP:
            job("threads_%s_t%d" % (name("NY", R0, D0), th), "build", build_gb("NY", R0, D0), th, [data("NY", R0, D0)],
                "%s build --data %s --name %s_t%d --prune points --timeout 21600 --threads %d --csv results/threads.csv"
                % (H, data("NY", R0, D0), name("NY", R0, D0), th, th))

elif phase == "mem":
    for g in GRAPHS:
        for r in ranks_of(g):
            job("mem_" + name(g, r, D0), "gen", 24, 1, [data(g, r, D0)],
                "%s mem --data %s --name %s --csv results/mem.csv" % (H, data(g, r, D0), name(g, r, D0)))
    for d in (HIST if "NY" in GRAPHS else []):
        if d != D0:
            job("mem_" + name("NY", R0, d), "gen", 24, 1, [data("NY", R0, d)],
                "%s mem --data %s --name %s --csv results/mem.csv" % (H, data("NY", R0, d), name("NY", R0, d)))
    for g in (REAL if os.environ.get("REAL", "1") == "1" else []):
        for r in RANKS:
            job("mem_%s_r%d" % (g, r), "gen", 4, 1, [real(g, r)],
                "%s mem --data %s --name %s_r%d --csv results/mem.csv" % (H, real(g, r), g, r))

elif phase == "time":
    builds, mems = {}, {}
    if os.path.exists("results/build.csv"):
        for row in csv.DictReader(open("results/build.csv")):
            if row["status"] == "OK":
                builds[(row["name"], row["prune"])] = row
    if os.path.exists("results/mem.csv"):
        for row in csv.DictReader(open("results/mem.csv")):
            mems[row["name"]] = row

    def kof(g, d, lo):  # weekdays in [lo, d): day 0 is a Monday
        return sum(1 for x in range(lo, d) if x % 7 < 5)

    def cch_fits(nm):
        m = mems.get(nm)
        return m is not None and float(m["cch_all_mb"]) <= BUDGET_MB

    def tjob(exp, g, r, d, method, sets, ph, k, extra="", prune="points", alphas=None, tag="", sfx="", report_name=None):
        nm = name(g, r, d, sfx)
        m = mems.get(name(g, r, d))  # the graph and N are the same for the noise variants
        world = float(m["cch_world_mb"]) / 1e3 if m else 0.5
        ver = k * world  # CCH reference for verification
        deps = [data(g, r, d, sfx), "../data/queries/%s/D10.txt" % g]
        if method.startswith("sighop"):
            b = builds.get((nm, prune))
            if b is None:
                return  # index missing (TO/MO); reported by collect.py
            mem = float(b["label_mb"]) / 1e3 * 1.3 + ver + 6
            extra += " --index " + idx(g, r, d, prune, sfx)
            deps.append(idx(g, r, d, prune, sfx))
        elif method in ("cch", "cchbatch"):
            if not cch_fits(name(g, r, d)):
                return
            mem = ver * (2 if method == "cchbatch" else 1) + 6
        elif method == "ttp":
            mem = k * float(m["n"] if m else 1e6) * 4 / 1e9 + ver + 6
        else:
            mem = ver + 6
        per = "" if "--per" in extra else {"dijkstra": " --per 100", "ttp": " --per 20"}.get(method, " --per 1000")
        if method == "dijkstra":
            extra += per + " --rounds 1 --warmup 1"
        elif method == "ttp":
            extra += per + " --rounds 1 --warmup 0 --timeout 60 --skip-after 5"
        else:
            extra += per + " --rounds 3 --warmup 16"
        if alphas:
            extra += " --alphas " + alphas
        jid = "%s_%s_%s%s%s" % (exp, nm, method, "_" + prune if prune != "points" else "", tag)
        job(jid, "time", mem, 1, deps,
            "%s query --data %s --name %s --method %s --qdir ../data/queries/%s --sets %s --phi %s%s --csv results/%s.csv"
            % (H, data(g, r, d, sfx), report_name or nm, method, g, sets, ph, extra, exp))

    k0 = kof("", D0, 0)
    if EXTRA:
        ex=EXTRA.lower()
        if EXTRA=='X1':
            for g in GRAPHS:
                for method in ['sighop','cchbatch','cch','dijkstra']:
                    tjob(ex,g,R0,D0,method,SETS,phi('wk17',0,D0,68),k0)
        elif EXTRA=='X2':
            for g in [x for x in GRAPHS if x in ('BAY','COL','FLA')]:
                for w in [7,30,90,180,365]:
                    for method in ['sighop','cchbatch','cch','dijkstra']:
                        tjob(ex,g,R0,365,method,'D10',phi('wk%d'%w,365-w,365),kof(g,365,365-w),
                             extra=' --per 20' if method=='dijkstra' else '',tag='_w%d'%w)
        elif EXTRA=='X3' and 'BAY' in GRAPHS:
            for r in (2,6,8):
                for method in ['sighop','cchbatch','cch','dijkstra']:
                    tjob(ex,'BAY',r,D0,method,'D10',phi('wk90',0,D0),k0)
        elif EXTRA=='X4' and 'FLA' in GRAPHS:
            for method in ['sighop','sighop-nohub','sighop-countsep']:
                tjob(ex,'FLA',R0,D0,method,SETS,phi('wk90',0,D0),k0)
        elif EXTRA=='X5' and 'BAY' in GRAPHS:
            tjob(ex,'BAY',R0,D0,'sighop',SETS,phi('wk90',0,D0),k0,prune='box')
        elif EXTRA=='X6' and 'NY' in GRAPHS:
            for parameter in (1,3,5):
                tjob(ex,'NY',R0,D0,'ttp','D1,D2,D3,D4',phi('wk90',0,D0),k0,
                     extra=' --ttp-k %d'%parameter,tag='_ttpk%d'%parameter,report_name=name('NY',R0,D0)+'_ttpk%d'%parameter)
        elif EXTRA=='X7' and 'NY' in GRAPHS:
            for method in ['sighop','cchbatch','cch']:
                tjob(ex,'NY',8,D0,method,SETS,phi('wk90',0,D0),k0)
        print('\n'.join(J));sys.exit(0)
    for g in GRAPHS:  # Exp-1: query distance
        for mth in ["sighop", "cchbatch", "cch", "dijkstra", "ttp"]:
            tjob("e1", g, R0, D0, mth, SETS, phi("wk90", 0, D0), k0)
    for g in GRAPHS:  # Exp-2: number of selected days (windows of the last w days)
        for w in WINDOWS:
            if w == D0:
                continue
            for mth in ["sighop", "cchbatch", "cch", "dijkstra"]:
                tjob("e2", g, R0, D0, mth, "D10", phi("wk%d" % w, D0 - w, D0), kof("", D0, D0 - w), tag="_w%d" % w)
    for g in [x for x in GRAPHS if x == "NY"]:  # Exp-3: rank (NY only)
        for r in RANKS:
            if r == R0:
                continue
            for mth in ["sighop", "cchbatch", "cch", "dijkstra"]:
                tjob("e3", g, r, D0, mth, "D10", phi("wk90", 0, D0), k0)
    if "NY" not in GRAPHS:
        print("\n".join(J)); sys.exit(0)
    for na, al in ALPHAS.items():  # Exp-4: number of quantile levels (NY)
        if na == 3:
            continue
        for mth in ["sighop", "cchbatch"]:
            tjob("e4", "NY", R0, D0, mth, "D10", phi("wk90", 0, D0), k0, alphas=al, tag="_a%d" % na)
    for mth in ["sighop-nohub", "sighop-countsep"]:  # Exp-5: query pruning ablation (NY)
        tjob("e5", "NY", R0, D0, mth, SETS, phi("wk90", 0, D0), k0)
    for d in HIST:  # Exp-7: history length (NY), all weekdays of the history
        if d == D0:
            continue
        for mth in ["sighop", "cchbatch", "cch", "dijkstra"]:
            extra = " --per 20" if mth == "dijkstra" else ""
            tjob("e7", "NY", R0, d, mth, "D10", phi("wk%d" % d, 0, d), kof("", d, 0), extra=extra)
    tjob("e8", "NY", R0, D0, "sighop", SETS, phi("wk90", 0, D0), k0, prune="box")  # Exp-8: box-only index

    def kdow(d, dows):
        lo, hi = map(int, dows.split("-"))
        return sum(1 for x in range(d) if lo <= x % 7 <= hi)

    for hh, sl in SLOTS_DAY.items():  # Exp-9: time description (slot of the day, day type), NY, D10
        for dn, dw in DOWS.items():
            if (dn != "wk" and hh != "08") or (dn == "wk" and hh == "08"):
                continue  # weekdays at every slot; weekends and all days at 08:00; wk/08 is Exp-1
            for mth in ["sighop", "cchbatch", "cch"]:
                tjob("e9", "NY", R0, D0, mth, "D10", phi("%s%s" % (dn, hh), 0, D0, sl, dw), kdow(D0, dw),
                     tag="_%s%s" % (dn, hh))
    for nz in NOISE:  # Exp-11: day-level noise of the semi-synthetic traffic, NY, D10
        if nz == N0:
            continue
        for mth in ["sighop", "cchbatch"]:
            tjob("e11", "NY", R0, D0, mth, "D10", phi("wk90", 0, D0), k0, sfx=nsfx(nz))

else:
    sys.exit(__doc__)

print("\n".join(J))
