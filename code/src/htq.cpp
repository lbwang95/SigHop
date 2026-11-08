// Experiment driver for the revised SigHop (see ../PLAN.md). Subcommands append rows to CSV files.
//   htq genq     --data X --out DIR [--sets 10 --per 1000 --lmin 60 --seed 7]
//   htq build    --data X --name NAME --prune points|box --index OUT.idx --timeout SEC --threads T --csv build.csv
//   htq query    --data X --name NAME --method M --qdir DIR --phi SPEC [--index IDX] --csv query.csv
//   htq mem      --data X --name NAME --csv mem.csv
//   htq selftest --data SMALL.bin
#include "baselines.hpp"

static vector<double> parseAlphas(const string& s) {
    vector<double> v;
    for (const string& x : splitStr(s, ',')) v.push_back(stod(x));
    sort(v.begin(), v.end());
    return v;
}
static vector<int> ranksOf(const vector<double>& alphas, int k) {
    vector<int> r;
    for (double al : alphas) r.push_back(quantileRank(al, k));
    return r;
}
static vector<float> worldMean(const Dataset& d) {
    vector<double> m(d.r, 0);
    for (int w = 0; w < d.N; w++)
        for (int j = 0; j < d.r; j++) m[j] += d.Z(w)[j];
    vector<float> out(d.r);
    for (int j = 0; j < d.r; j++) out[j] = (float)(m[j] / d.N);
    return out;
}
static void selectInto(vector<float>& vals, const vector<int>& ranks, float* out) {
    static thread_local vector<float> sel;
    sel = vals;
    multiSelect(sel.data(), (int)sel.size(), ranks, 0, (int)ranks.size());
    for (size_t i = 0; i < ranks.size(); i++) out[i] = sel[ranks[i]];
}
static bool close(float a, float b) { return fabs(a - b) <= 1e-3f * max(1.f, fabs(b)); }

// Query sets D_1..D_k with increasing travel time at the mean world, as in the PCSP paper.
static int cmdGenq(const Args& a) {
    Dataset d = loadDataset(a.need("data"));
    GraphCsr g;
    g.build(d);
    Dijkstra dj(d, g);
    vector<float> zbar = worldMean(d);
    vector<float> wt(d.m);
    for (int64_t e = 0; e < d.m; e++) wt[e] = dotf(d.U(e), zbar.data(), (int)d.r);
    auto W = [&](int e) { return (double)wt[e]; };
    int sets = a.geti("sets", 10), per = a.geti("per", 1000);
    double lmin = a.getd("lmin", 60);
    mt19937 rng(a.geti("seed", 7));
    uniform_int_distribution<int> pick(0, (int)d.n - 1);
    double lmax = 0;
    for (int i = 0; i < 20; i++) {
        dj.oneToAll(pick(rng), W);
        for (int v = 0; v < d.n; v++)
            if (dj.sf[v] == dj.stamp) lmax = max(lmax, dj.df[v]);
    }
    double x = pow(lmax / lmin, 1.0 / sets);
    vector<vector<pair<int, int>>> q(sets);
    int rounds = 0, maxRounds = a.geti("max-sources", 20000);
    bool done = false;
    while (!done && rounds < maxRounds) {
        rounds++;
        int s = pick(rng);
        dj.oneToAll(s, W);
        vector<vector<int>> cand(sets);
        for (int v = 0; v < d.n; v++) {
            if (dj.sf[v] != dj.stamp || dj.df[v] <= lmin) continue;
            int b = (int)floor(log(dj.df[v] / lmin) / log(x));
            if (b >= 0 && b < sets) cand[b].push_back(v);
        }
        done = true;
        for (int b = 0; b < sets; b++) {
            if ((int)q[b].size() < per && !cand[b].empty()) {
                for (int rep = 0; rep < 3 && (int)q[b].size() < per; rep++)
                    q[b].push_back({s, cand[b][rng() % cand[b].size()]});
            }
            if ((int)q[b].size() < per) done = false;
        }
    }
    string dir = a.need("out");
    for (int b = 0; b < sets; b++) {
        ofstream f(dir + "/D" + to_string(b + 1) + ".txt");
        for (auto& pr : q[b]) f << pr.first << " " << pr.second << "\n";
        printf("D%d: %zu queries in (%.1f, %.1f] s\n", b + 1, q[b].size(), lmin * pow(x, b), lmin * pow(x, b + 1));
    }
    return 0;
}

static string hdrBuild =
    "name,r,days,N,n,m,tw,height_levels,prune,threads,td_s,shortcut_s,label_s,build_s,sig_mb,min_mb,off_mb,"
    "label_mb,shortcut_mb,z_mb,peak_rss_mb,set_mean,set_p50,set_p90,set_p99,set_max,boxkept_mean,boxkept_max,status";

// Index construction of SigHop with its statistics (index-cost table, Exp-6, Exp-7, Exp-8).
static int cmdBuild(const Args& a) {
    if (a.kv.count("stats-csv")) { fprintf(stderr, "V3_STAGE load_dataset\n"); fflush(stderr); }
    Dataset d = loadDataset(a.need("data"));
    string name = a.get("name", a.need("data"));
    string prune = a.get("prune", "points");
    if (prune != "points" && prune != "box" && prune != "pointsonly")
        throw runtime_error("--prune must be points, box, or pointsonly");
    TreeDecomp T;
    if (a.kv.count("stats-csv")) { fprintf(stderr, "V3_STAGE tree_decomposition\n"); fflush(stderr); }
    T.build(d);
    vector<int> worlds(d.N);
    iota(worlds.begin(), worlds.end(), 0);
    SigIndex idx(d, T);
    if (a.kv.count("stats-csv")) { fprintf(stderr, "V3_STAGE index_setup\n"); fflush(stderr); }
    idx.setup(worlds, prune != "box", prune != "pointsonly");
    idx.collectStats = a.kv.count("stats-csv") != 0;
    idx.timeoutSec = a.getd("timeout", 0);
    Timer tb;
    idx.build();
    double bsec = tb.sec();
    string status = idx.timedOut ? "TO" : "OK";
    SizeStats ls;
    double sigMb = 0, minMb = 0, offMb = 0, scMb = idx.shortcutBytes() / 1e6;
    if (!idx.timedOut) {
        ls = idx.labelSizeStats();
        sigMb = idx.sigBytes() / 1e6;
        minMb = idx.minBytes() / 1e6;
        offMb = idx.offBytes() / 1e6;
        if (a.kv.count("index")) {
            if (idx.collectStats) { fprintf(stderr, "V3_STAGE save_index\n"); fflush(stderr); }
            idx.dropShortcuts();
            idx.save(a.get("index"));
        }
    }
    double bkMean = idx.pstat.calls ? (double)idx.pstat.boxKeptSum / idx.pstat.calls : 0;
    Csv csv(a.need("csv"), hdrBuild);
    csv.row(joinCsv(name, d.r, d.days, d.N, d.n, d.m, T.width - 1, T.height + 1, prune, numThreads(), T.buildSec,
                    idx.shortcutSec, idx.labelSec, bsec + T.buildSec, sigMb, minMb, offMb, sigMb + minMb + offMb, scMb,
                    d.N * d.r * 4 / 1e6, peakRssMb(), ls.mean, ls.p50, ls.p90, ls.p99, ls.maxv, bkMean,
                    idx.pstat.boxKeptMax, status));
    if (a.kv.count("stats-csv")) {
        Csv stats(a.need("stats-csv"), "name,prune,threads,N,calls,in_mean,in_max,pre_mean,boxkept_mean,boxkept_max,"
                  "out_mean,out_max,shortcut_s,label_s,build_s,status,phase");
        auto emit = [&](const PruneStats& s, const string& phase, const string& state) {
            double calls = (double)max<uint64_t>(1, s.calls);
            stats.row(joinCsv(name, prune, numThreads(), d.N, s.calls, s.inSum / calls, s.inMax, s.preSum / calls,
                              s.boxKeptSum / calls, s.boxKeptMax, s.outSum / calls, s.outMax,
                              idx.shortcutSec, idx.labelSec, bsec + T.buildSec, state, phase));
        };
        emit(idx.shortcutStats, "shortcut", idx.shortcutCompleted ? "OK" : status);
        emit(idx.labelStats, "label", !idx.labelStarted ? "NOT_RUN" : idx.labelCompleted ? "OK" : status);
        emit(idx.pstat, "all", status);
    }
    printf("%s r=%lld N=%lld prune=%s status=%s build=%.1fs labels=%.1fMB mean_set=%.2f max_set=%zu\n", name.c_str(),
           (long long)d.r, (long long)d.N, prune.c_str(), status.c_str(), bsec + T.buildSec, sigMb + minMb + offMb, ls.mean,
           ls.maxv);
    return idx.timedOut ? 3 : 0;
}

static string hdrQuery =
    "name,r,days,N,k,levels,phi,method,set,queries,completed,timeouts,rounds,mean_us,min_us,max_us,prep_s,custom_ms,"
    "hubs,hubs_kept,cands,cands_kept,verified,mismatches,peak_rss_mb,status";

// Query time of one method on the query sets (Exp-1 to Exp-5). Timing is single-threaded; the
// answers are verified afterwards against CCH (per snapshot), or against Dijkstra for CCH itself.
static int cmdQuery(const Args& a) {
    Dataset d = loadDataset(a.need("data"));
    string name = a.get("name", a.need("data"));
    string method = a.need("method");
    Phi phi = parsePhi(a.need("phi"));
    vector<int> worlds = worldsOf(d, phi);
    if (worlds.empty()) throw runtime_error("the time description selects no snapshot");
    for (int w : worlds)
        if (d.wslot[w] != d.wslot[worlds[0]]) throw runtime_error("--phi must select a single slot");
    const int k = (int)worlds.size(), r = (int)d.r;
    vector<double> alphas = parseAlphas(a.get("alphas", "0.5,0.8,0.95"));
    vector<int> ranks = ranksOf(alphas, k);
    const int A = (int)ranks.size();
    TreeDecomp T;
    T.build(d);
    GraphCsr g;
    g.build(d);
    unique_ptr<SigIndex> idx;
    QueryOpts qo;
    bool isSig = method.rfind("sighop", 0) == 0;
    if (isSig) {
        idx = make_unique<SigIndex>(d, T);
        idx->load(a.need("index"));
        if (method == "sighop-nohub") qo.hubPrune = false;
        else if (method == "sighop-countsep") qo.costSep = false;
        else if (method == "sighop-ffref") qo.ffRef = true;
        else if (method != "sighop") throw runtime_error("unknown method " + method);
    } else if (method != "cch" && method != "cchbatch" && method != "dijkstra" && method != "ttp") {
        throw runtime_error("unknown method " + method);
    }
    CCH cch(d, T);
    CCHBatch cb(d, T);
    vector<vector<float>> Wk;
    double prepSec = 0, customMs = 0;
    auto customizeAll = [&]() {
        Timer tp;
        Wk.assign(k, {});
        for (int i = 0; i < k; i++) cch.customize(d.Z(worlds[i]), Wk[i]);
        return tp.sec();
    };
    if (method == "cch") {
        prepSec = customizeAll();
        customMs = prepSec * 1000 / k;
    } else if (method == "cchbatch") {
        prepSec = customizeAll();
        customMs = prepSec * 1000 / k;
        Timer tp;
        cb.prepare(Wk);
        prepSec += tp.sec();
        vector<vector<float>>().swap(Wk);
    }
    Dijkstra dj(d, g);
    TTPAdapt ttp(d, g);
    ttp.timeoutSec = a.getd("timeout", 60);
    ttp.kTTP = a.geti("ttp-k", 3);
    OmegaCtx ctx;
    vector<float> vals(k);
    QueryStats qst;
    auto run = [&](int s, int t, float* out, bool stats) -> bool {
        if (isSig) {
            ctx.prepare(d, worlds.data(), k);
            idx->answer(s, t, ctx, ranks, qo, out, vals, stats ? &qst : nullptr);
            return true;
        }
        if (method == "cch") {
            for (int i = 0; i < k; i++) vals[i] = cch.query(s, t, Wk[i]);
        } else if (method == "cchbatch") {
            cb.query(s, t, vals.data());
        } else if (method == "dijkstra") {
            for (int i = 0; i < k; i++) {
                const float* z = d.Z(worlds[i]);
                vals[i] = (float)dj.biQuery(s, t, [&](int e) { return dotd(d.U(e), z, r); });
            }
        } else {
            if (!ttp.query(s, t, worlds, vals.data())) return false;
        }
        selectInto(vals, ranks, out);
        return true;
    };
    vector<string> sets = splitStr(a.get("sets", "D1,D2,D3,D4,D5,D6,D7,D8,D9,D10"), ',');
    int per = a.geti("per", 1000), rounds = a.geti("rounds", 3), warmup = a.geti("warmup", 16);
    int skipAfter = a.geti("skip-after", 5), verifyPer = a.geti("verify-per", 1000000);
    bool verify = a.geti("verify", 1) != 0;
    Csv csv(a.need("csv"), hdrQuery);
    for (const string& set : sets) {
        auto qs = loadQueries(a.need("qdir") + "/" + set + ".txt");
        if ((int)qs.size() > per) qs.resize(per);
        int nq = (int)qs.size();
        vector<float> ans((size_t)nq * A, FINF);
        vector<char> to(nq, 0);
        for (int i = 0; i < min(warmup, nq); i++) run(qs[i].first, qs[i].second, &ans[(size_t)i * A], false);
        qst = QueryStats();
        vector<double> roundUs;
        string status = "OK";
        int timeouts = 0;
        for (int rd = 0; rd < rounds; rd++) {
            Timer tm;
            timeouts = 0;
            int done = 0;
            for (int i = 0; i < nq; i++) {
                bool ok = run(qs[i].first, qs[i].second, &ans[(size_t)i * A], rd == 0);
                to[i] = !ok;
                timeouts += !ok;
                done++;
                if (skipAfter > 0 && i + 1 == skipAfter && timeouts == skipAfter) {
                    status = "TO";
                    break;
                }
            }
            roundUs.push_back(tm.us() / max(1, done));
            if (status == "TO") break;
        }
        int completed = 0;
        for (int i = 0; i < nq; i++) completed += !to[i];
        if (status == "TO") completed = 0;
        int verified = 0, mismatches = 0;
        if (verify && status == "OK") {
            int nv = min(nq, verifyPer);
            vector<float> ref(A);
            if (method == "cch") {
                nv = min(nv, a.geti("verify-dijkstra", 10));
                for (int i = 0; i < nv; i++) {
                    for (int w = 0; w < k; w++) {
                        const float* z = d.Z(worlds[w]);
                        vals[w] = (float)dj.biQuery(qs[i].first, qs[i].second, [&](int e) { return dotd(d.U(e), z, r); });
                    }
                    selectInto(vals, ranks, ref.data());
                    verified++;
                    for (int j = 0; j < A; j++)
                        if (!close(ans[(size_t)i * A + j], ref[j])) {
                            mismatches++;
                            break;
                        }
                }
            } else {
                if (Wk.empty()) customizeAll();
                for (int i = 0; i < nv; i++) {
                    if (to[i]) continue;
                    for (int w = 0; w < k; w++) vals[w] = cch.query(qs[i].first, qs[i].second, Wk[w]);
                    selectInto(vals, ranks, ref.data());
                    verified++;
                    for (int j = 0; j < A; j++)
                        if (!close(ans[(size_t)i * A + j], ref[j])) {
                            mismatches++;
                            break;
                        }
                }
                if (method != "cch" && method != "cchbatch" && !isSig) vector<vector<float>>().swap(Wk);
            }
        }
        sort(roundUs.begin(), roundUs.end());
        double med = roundUs.empty() ? 0 : roundUs[roundUs.size() / 2];
        double q = max(1, nq);
        csv.row(joinCsv(name, d.r, d.days, d.N, k, A, phi.name, method, set, nq, completed, timeouts, (int)roundUs.size(), med,
                        roundUs.empty() ? 0 : roundUs.front(), roundUs.empty() ? 0 : roundUs.back(), prepSec, customMs,
                        qst.hubs / q, qst.hubsKept / q, qst.cands / q, qst.candsKept / q, verified, mismatches, peakRssMb(),
                        status));
        printf("%s %s %s k=%d: %.3f us/query (%s) mismatches=%d/%d\n", name.c_str(), method.c_str(), set.c_str(), k, med,
               status.c_str(), mismatches, verified);
        fflush(stdout);
    }
    return 0;
}

static string hdrMem =
    "name,r,days,N,n,m,tw,height_levels,num_shortcuts,cch_world_mb,cch_all_mb,cch_custom_ms,label2hop_world_mb,"
    "label2hop_all_mb,model_mb";

// Space of the per-snapshot indexes: one CCH weight array, and one 2-hop distance labeling with one
// distance per vertex and ancestor in each direction (H2H/CTL-style), times N.
static int cmdMem(const Args& a) {
    Dataset d = loadDataset(a.need("data"));
    string name = a.get("name", a.need("data"));
    TreeDecomp T;
    T.build(d);
    CCH cch(d, T);
    vector<float> w;
    int reps = (int)min<int64_t>(d.N, a.geti("reps", 5));
    Timer tm;
    for (int i = 0; i < reps; i++) cch.customize(d.Z((int)(i * (d.N / max(1, reps)))), w);
    double cms = tm.us() / 1000 / max(1, reps);
    double cchW = 2.0 * T.numShortcuts * 4 / 1e6;
    double lab = 0;
    for (int v = 0; v < T.n; v++) lab += 2.0 * T.depth[v] * 4;
    lab /= 1e6;
    Csv csv(a.need("csv"), hdrMem);
    csv.row(joinCsv(name, d.r, d.days, d.N, d.n, d.m, T.width - 1, T.height + 1, T.numShortcuts, cchW, cchW * d.N, cms,
                    lab, lab * d.N, (d.m + d.N) * d.r * 4 / 1e6));
    printf("%s: CCH %.1f MB/snapshot (all %.1f GB), customization %.1f ms, 2-hop %.1f MB/snapshot (all %.1f GB)\n",
           name.c_str(), cchW, cchW * d.N / 1e3, cms, lab, lab * d.N / 1e3);
    return 0;
}

// Brute-force checks on a small dataset: every method against one Dijkstra search per snapshot.
static int cmdSelftest(const Args& a) {
    Dataset d = loadDataset(a.need("data"));
    const int r = (int)d.r;
    TreeDecomp T;
    T.build(d);
    GraphCsr g;
    g.build(d);
    Dijkstra dj(d, g);
    vector<int> all(d.N);
    iota(all.begin(), all.end(), 0);
    SigIndex ip(d, T), ib(d, T);
    ip.setup(all, true);
    ib.setup(all, false);
    ip.build();
    ib.build();
    string tmp = a.get("tmp", "/tmp/htq_selftest.idx");
    ip.dropShortcuts();
    ip.save(tmp);
    SigIndex il(d, T);
    il.load(tmp);
    unique_ptr<SigIndex> only, onlyLoaded;
    if (a.geti("pointsonly-tests", 0)) {
        only = make_unique<SigIndex>(d, T);
        only->setup(all, true, false);
        only->collectStats = true;
        only->build();
        if (only->pstat.inSum != only->pstat.preSum || only->pstat.preSum != only->pstat.boxKeptSum)
            throw runtime_error("pointsonly unexpectedly filtered candidates before point tests");
        only->dropShortcuts();
        only->save(tmp + ".pointsonly");
        onlyLoaded = make_unique<SigIndex>(d, T);
        onlyLoaded->load(tmp + ".pointsonly");
    }
    CCH cch(d, T);
    CCHBatch cb(d, T);
    TTPAdapt ttp(d, g);
    ttp.timeoutSec = a.getd("timeout", 20);
    mt19937 rng(a.geti("seed", 3));
    int nq = a.geti("queries", 150), fails = 0, checks = 0, ttpChecks = 0;
    set<int> slotSet(d.wslot.begin(), d.wslot.end());
    vector<int> slots(slotSet.begin(), slotSet.end());
    vector<string> specs;
    for (int i = 0; i < min<int>(3, slots.size()); i++) {
        int sl = slots[(i * 7) % slots.size()];
        specs.push_back("wk:dow=0-4;slots=" + to_string(sl) + "-" + to_string(sl + 1));
        specs.push_back("few:days=0-3;slots=" + to_string(sl) + "-" + to_string(sl + 1));
    }
    specs.push_back("one:days=0-1;slots=" + to_string(slots[0]) + "-" + to_string(slots[0] + 1));
    vector<double> alphas = {0.5, 0.8, 0.95};
    for (const string& spec : specs) {
        Phi phi = parsePhi(spec);
        vector<int> ws = worldsOf(d, phi);
        if (ws.empty()) continue;
        int k = (int)ws.size();
        vector<int> ranks = ranksOf(alphas, k);
        vector<vector<float>> Wk(k);
        for (int i = 0; i < k; i++) cch.customize(d.Z(ws[i]), Wk[i]);
        cb.prepare(Wk);
        OmegaCtx ctx;
        ctx.prepare(d, ws.data(), k);
        vector<float> ref(k), vals(k), out(3), bv(k);
        uniform_int_distribution<int> pick(0, (int)d.n - 1);
        for (int q = 0; q < nq; q++) {
            int s = pick(rng), t = q % 17 == 0 ? s : pick(rng);
            for (int i = 0; i < k; i++) {
                const float* z = d.Z(ws[i]);
                ref[i] = (float)dj.biQuery(s, t, [&](int e) { return dotd(d.U(e), z, r); });
            }
            auto cmp = [&](const vector<float>& v, const char* what) {
                checks++;
                for (int i = 0; i < k; i++)
                    if (!close(v[i], ref[i])) {
                        if (fails < 20)
                            fprintf(stderr, "FAIL %s %s s=%d t=%d world=%d got=%g ref=%g\n", what, spec.c_str(), s, t, ws[i],
                                    v[i], ref[i]);
                        fails++;
                        return;
                    }
            };
            vector<const SigIndex*> idxs = {&ip, &ib, &il};
            vector<string> inames = {"points", "box", "loaded"};
            if (only) {
                idxs.push_back(only.get()); inames.push_back("pointsonly");
                idxs.push_back(onlyLoaded.get()); inames.push_back("pointsonly-loaded");
            }
            for (int x = 0; x < (int)idxs.size(); x++) {
                for (int variant = 0; variant < (a.geti("ffref-tests", 0) ? 4 : 3); variant++) {
                    QueryOpts o;
                    if (variant == 1) o.hubPrune = false;
                    if (variant == 2) o.costSep = false;
                    if (variant == 3) o.ffRef = true;
                    idxs[x]->answer(s, t, ctx, ranks, o, out.data(), vals, nullptr);
                    string nm = string("sighop-") + inames[x] + "-v" + to_string(variant);
                    cmp(vals, nm.c_str());
                }
            }
            for (int i = 0; i < k; i++) vals[i] = cch.query(s, t, Wk[i]);
            cmp(vals, "cch");
            cb.query(s, t, bv.data());
            cmp(bv, "cchbatch");
            if (q % 10 == 0 && k <= 8) {
                if (ttp.query(s, t, ws, vals.data())) {
                    cmp(vals, "ttp");
                    ttpChecks++;
                }
            }
        }
    }
    printf("selftest: %d checks (%d TTP), %d failures; mean set size points=%.2f box=%.2f\n", checks, ttpChecks, fails,
           ip.pstat.calls ? (double)ip.pstat.outSum / ip.pstat.calls : 0, ib.pstat.calls ? (double)ib.pstat.outSum / ib.pstat.calls : 0);
    if (fails == 0) printf("ALL TESTS PASSED\n");
    return fails == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) throw runtime_error("usage: htq genq|build|query|mem|selftest [--options]");
        string cmd = argv[1];
        Args a(argc, argv, 2);
#ifdef _OPENMP
        omp_set_num_threads(a.geti("threads", 1));
#endif
        if (cmd == "genq") return cmdGenq(a);
        if (cmd == "build") return cmdBuild(a);
        if (cmd == "query") return cmdQuery(a);
        if (cmd == "mem") return cmdMem(a);
        if (cmd == "selftest") return cmdSelftest(a);
        throw runtime_error("unknown command " + cmd);
    } catch (const bad_alloc&) {
        fprintf(stderr, "error: out of memory\n");
        return 4;
    } catch (const exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
