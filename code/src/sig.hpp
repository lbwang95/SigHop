// SigHop: signature labels on a tree decomposition (algorithm.tex, Sections IV-VI).
//   Index pruning (SigPrune, Algorithm 1): box filter over the box R of all snapshots, then one
//   minimizer per historical snapshot ("points"). With prune=box only the box filter runs.
//   Each set is kept in increasing order of <a, zbar> (zbar = mean of all factor vectors), so the
//   first signature is the representative p, and the componentwise minimum m is stored per set.
//   Query (Algorithm 2): separator by cost, reference c* from the representatives, removal of the
//   separator vertices with L_B(m_sh + m_ht - c*) >= 0, exact evaluation of the rest.
#pragma once
#include "td.hpp"

inline float dotf(const float* a, const float* b, int r) {
    float s = 0;
    for (int j = 0; j < r; j++) s += a[j] * b[j];
    return s;
}

// L_B(x) = min over z in [lo, hi] of <x, z>  (Equation lb).
inline float lowerBox(const float* x, const float* lo, const float* hi, int r) {
    float s = 0;
    for (int j = 0; j < r; j++) s += x[j] > 0 ? lo[j] * x[j] : hi[j] * x[j];
    return s;
}
inline float upperBox(const float* x, const float* lo, const float* hi, int r) {
    float s = 0;
    for (int j = 0; j < r; j++) s += x[j] > 0 ? hi[j] * x[j] : lo[j] * x[j];
    return s;
}
// a dominates b over the box: L_B(b - a) >= 0.
inline bool dominatesBox(const float* a, const float* b, const float* lo, const float* hi, int r) {
    float s = 0;
    for (int j = 0; j < r; j++) {
        float d = b[j] - a[j];
        s += d > 0 ? lo[j] * d : hi[j] * d;
    }
    return s >= 0;
}

struct Box {
    int r = 0;
    vector<float> lo, hi, center;
    void fit(const Dataset& d, const vector<int>& worlds) {
        r = (int)d.r;
        lo.assign(r, FINF);
        hi.assign(r, -FINF);
        vector<double> mean(r, 0);
        for (int w : worlds) {
            const float* z = d.Z(w);
            for (int j = 0; j < r; j++) {
                lo[j] = min(lo[j], z[j]);
                hi[j] = max(hi[j], z[j]);
                mean[j] += z[j];
            }
        }
        center.assign(r, 0);
        for (int j = 0; j < r; j++) center[j] = (float)(mean[j] / max<size_t>(1, worlds.size()));
    }
};

// Factor vectors of the historical snapshots in structure-of-arrays form: soa[j * N + w].
struct PointSet {
    int r = 0;
    int64_t N = 0;
    vector<float> soa;
    void fit(const Dataset& d, const vector<int>& worlds) {
        r = (int)d.r;
        N = (int64_t)worlds.size();
        soa.resize((size_t)N * r);
        for (int64_t i = 0; i < N; i++)
            for (int j = 0; j < r; j++) soa[(size_t)j * N + i] = d.Z(worlds[i])[j];
    }
};

struct PruneStats {
    uint64_t calls = 0, boxKeptSum = 0, boxKeptMax = 0, outSum = 0, outMax = 0;
    uint64_t inSum = 0, inMax = 0, preSum = 0;
    void add(const PruneStats& o) {
        calls += o.calls;
        boxKeptSum += o.boxKeptSum;
        boxKeptMax = max(boxKeptMax, o.boxKeptMax);
        outSum += o.outSum;
        outMax = max(outMax, o.outMax);
        inSum += o.inSum;
        inMax = max(inMax, o.inMax);
        preSum += o.preSum;
    }
};

// Optional, allocation-free progress evidence for TO/MO (not used by query algorithms).
struct PruneProgress {
    atomic<int> phase{0}, completed{0};
    atomic<uint64_t> entered[2], requested[2];
    PruneProgress() {
        for (int i = 0; i < 2; i++) { entered[i].store(0); requested[i].store(0); }
    }
    void log(const char* event) const {
        int p = phase.load(memory_order_relaxed);
        fprintf(stderr, "V3_PROGRESS phase=%s event=%s max_candidates=%llu max_requested=%llu completed_mask=%d\n",
                p == 0 ? "shortcut" : "label", event,
                (unsigned long long)entered[p].load(memory_order_relaxed),
                (unsigned long long)requested[p].load(memory_order_relaxed), completed.load(memory_order_relaxed));
        fflush(stderr);
    }
    void observe(uint64_t count, bool input) {
        int p = phase.load(memory_order_relaxed);
        auto& value = input ? entered[p] : requested[p];
        uint64_t previous = value.load(memory_order_relaxed);
        while (count > previous) {
            if (value.compare_exchange_weak(previous, count, memory_order_relaxed)) {
                log(input ? "input_peak" : "allocation_peak");
                break;
            }
        }
    }
};

// SigPrune (Algorithm 1). The result stays sorted by the value at the center.
struct Pruner {
    int r = 0;
    const float* lo = nullptr;
    const float* hi = nullptr;
    const float* center = nullptr;
    const PointSet* pts = nullptr;  // null: box filter only
    bool useBox = true;
    PruneProgress* progress = nullptr;
    PruneStats st;
    vector<pair<pair<float, float>, int>> ord;
    vector<float> kept, best, row;
    vector<int> arg;
    vector<char> keep;

    void run(vector<float>& pool) {
        int cnt = (int)(pool.size() / r);
        if (cnt == 0) return;
        st.inSum += pool.size() / r;
        st.inMax = max<uint64_t>(st.inMax, pool.size() / r);
        if (progress) progress->observe(pool.size() / r, true);
        // Cheap safe pre-filter: c is dominated if its smallest value over the box exceeds the
        // largest value of another signature over the box.
        float minUb = FINF;
        if (useBox)
            for (int i = 0; i < cnt; i++) minUb = min(minUb, upperBox(&pool[i * r], lo, hi, r));
        ord.clear();
        for (int i = 0; i < cnt; i++) {
            const float* c = &pool[i * r];
            if (!useBox || lowerBox(c, lo, hi, r) <= minUb) ord.push_back({{dotf(c, center, r), c[0]}, i});
        }
        st.preSum += ord.size();
        sort(ord.begin(), ord.end());
        kept.clear();
        int nk = 0;
        for (auto& pr : ord) {
            const float* c = &pool[pr.second * r];
            bool dom = false;
            if (useBox)
                for (int k = 0; k < nk && !dom; k++) dom = dominatesBox(&kept[k * r], c, lo, hi, r);
            if (!dom) {
                kept.insert(kept.end(), c, c + r);
                nk++;
            }
        }
        st.calls++;
        st.boxKeptSum += nk;
        st.boxKeptMax = max<uint64_t>(st.boxKeptMax, nk);
        if (pts && nk > 1) {
            const int64_t N = pts->N;
            best.assign(N, FINF);
            arg.assign(N, 0);
            row.resize(N);
            const float* soa = pts->soa.data();
            for (int k = 0; k < nk; k++) {
                const float* c = &kept[k * r];
                float* t = row.data();
                const float c0 = c[0];
                for (int64_t w = 0; w < N; w++) t[w] = c0 * soa[w];
                for (int j = 1; j < r; j++) {
                    const float cj = c[j];
                    const float* col = soa + (size_t)j * N;
                    for (int64_t w = 0; w < N; w++) t[w] += cj * col[w];
                }
                float* b = best.data();
                int* a = arg.data();
                for (int64_t w = 0; w < N; w++) {
                    if (t[w] < b[w]) {  // strict: ties keep the earlier signature
                        b[w] = t[w];
                        a[w] = k;
                    }
                }
            }
            keep.assign(nk, 0);
            for (int64_t w = 0; w < N; w++) keep[arg[w]] = 1;
            int o = 0;
            for (int k = 0; k < nk; k++) {
                if (!keep[k]) continue;
                if (o != k) copy(&kept[k * r], &kept[k * r] + r, &kept[o * r]);
                o++;
            }
            kept.resize((size_t)o * r);
            nk = o;
        }
        st.outSum += nk;
        st.outMax = max<uint64_t>(st.outMax, nk);
        pool = kept;
    }
};

inline void joinInto(const float* a, size_t na, const float* b, size_t nb, int r, vector<float>& pool,
                     PruneProgress* progress = nullptr) {
    size_t base = pool.size();
    if (progress) progress->observe(base / r + na * nb, false);
    pool.resize(base + na * nb * r);
    float* out = pool.data() + base;
    for (size_t i = 0; i < na; i++) {
        for (size_t j = 0; j < nb; j++) {
            for (int k = 0; k < r; k++) out[k] = a[i * r + k] + b[j * r + k];
            out += r;
        }
    }
}

struct SizeStats {
    double mean = 0, sum = 0;
    size_t count = 0, p50 = 0, p90 = 0, p99 = 0, maxv = 0;
};
inline SizeStats sizeStats(const vector<uint64_t>& hist) {
    SizeStats st;
    for (size_t i = 0; i < hist.size(); i++) {
        st.count += hist[i];
        st.sum += (double)i * hist[i];
        if (hist[i]) st.maxv = i;
    }
    st.mean = st.count ? st.sum / st.count : 0;
    size_t acc = 0;
    bool b50 = false, b90 = false, b99 = false;
    for (size_t i = 0; i < hist.size(); i++) {
        acc += hist[i];
        if (!b50 && acc * 2 >= st.count) { st.p50 = i; b50 = true; }
        if (!b90 && acc * 10 >= st.count * 9) { st.p90 = i; b90 = true; }
        if (!b99 && acc * 100 >= st.count * 99) { st.p99 = i; b99 = true; }
    }
    return st;
}

// The selected snapshots of one query: factor vectors in SoA form, their box and mean.
struct OmegaCtx {
    int k = 0, r = 0;
    vector<float> zsoa;  // r arrays of length k
    vector<float> lo, hi, mean;
    void prepare(const Dataset& d, const int* idx, int cnt) {
        k = cnt;
        r = (int)d.r;
        zsoa.resize((size_t)k * r);
        lo.assign(r, FINF);
        hi.assign(r, -FINF);
        mean.assign(r, 0);
        for (int j = 0; j < r; j++) {
            float* col = &zsoa[(size_t)j * k];
            double sum = 0;
            float mn = FINF, mx = -FINF;
            for (int i = 0; i < k; i++) {
                float v = d.z[(size_t)idx[i] * d.r + j];
                col[i] = v;
                mn = min(mn, v);
                mx = max(mx, v);
                sum += v;
            }
            lo[j] = mn;
            hi[j] = mx;
            mean[j] = k ? (float)(sum / k) : 0;
        }
    }
};

struct QueryOpts {
    bool hubPrune = true;  // Line 8 of Algorithm 2
    bool costSep = true;   // Line 6: separator with fewer signatures (else fewer vertices)
    bool ffRef = false;    // V4 opt-in: choose the reference using free-flow costs.
};

// Per-query averages are taken by the caller; hubs/cands are counted on the chosen separator.
struct QueryStats {
    double hubs = 0, hubsKept = 0, cands = 0, candsKept = 0;
};

struct SigIndex {
    const Dataset& d;
    const TreeDecomp& T;
    int r;
    Box R;
    PointSet P;
    bool usePoints = true;
    bool useBox = true, collectStats = false;
    bool shortcutCompleted = false, labelStarted = false, labelCompleted = false;
    mutable PruneProgress progress;
    vector<vector<float>> S;  // S[2*sc+dir], dir 0: lower->upper, 1: upper->lower
    struct Label {
        vector<uint32_t> off;  // signatures toward the ancestor at depth i: sig[off[i] .. off[i+1])
        vector<float> sig;
        vector<float> mn;      // componentwise minimum of entry i at mn[i*r .. i*r+r)
    };
    vector<Label> L[2];        // L[0][v]: v -> ancestor, L[1][v]: ancestor -> v
    double shortcutSec = 0, labelSec = 0, timeoutSec = 0;
    atomic<bool> timedOut{false};
    PruneStats pstat;
    PruneStats shortcutStats, labelStats;

    SigIndex(const Dataset& data, const TreeDecomp& td) : d(data), T(td), r((int)data.r) {}

    void setup(const vector<int>& worlds, bool points, bool box = true) {
        R.fit(d, worlds);
        usePoints = points;
        useBox = box;
        if (points) P.fit(d, worlds);
    }
    void initPruner(Pruner& p) const {
        p.r = r;
        p.lo = R.lo.data();
        p.hi = R.hi.data();
        p.center = R.center.data();
        p.pts = usePoints ? &P : nullptr;
        p.useBox = useBox;
        p.progress = collectStats ? &progress : nullptr;
    }

    // Bottom-up phase: shortcut signature sets. Vertices are processed by their height in the
    // elimination tree (leaves first). A vertex only reads sets S(v,.) that are final, because every
    // update of S(v,a) comes from a descendant of v. Vertices of one height run in parallel, and the
    // updates of a shared target set are serialized by striped locks.
    void buildShortcuts() {
        Timer tm;
        S.assign(2 * T.numShortcuts, {});
        for (int64_t e = 0; e < d.m; e++) {
            if (T.arcPlace[e] == UINT64_MAX) continue;
            auto& set = S[T.arcPlace[e]];
            set.insert(set.end(), d.U(e), d.U(e) + r);
        }
        {
            Pruner pr;
            initPruner(pr);
            for (auto& set : S) pr.run(set);
            pstat.add(pr.st);
            if (collectStats) shortcutStats.add(pr.st);
        }
        vector<int> hgt(T.n, 0);
        for (int idx = 0; idx < T.n; idx++) {
            int v = T.order[idx];
            if (T.parent[v] >= 0) hgt[T.parent[v]] = max(hgt[T.parent[v]], hgt[v] + 1);
        }
        int maxH = 0;
        for (int v = 0; v < T.n; v++) maxH = max(maxH, hgt[v]);
        vector<vector<int>> byH(maxH + 1);
        for (int idx = 0; idx < T.n; idx++) byH[hgt[T.order[idx]]].push_back(T.order[idx]);
        const int NL = 1 << 16;
        vector<mutex> locks(NL);
        for (int h = 0; h <= maxH; h++) {
            if (timeoutSec > 0 && tm.sec() > timeoutSec) {
                timedOut = true;
                return;
            }
            const auto& vs = byH[h];
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                Pruner pr;
                initPruner(pr);
                vector<float> pool;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif
                for (int ii = 0; ii < (int)vs.size(); ii++) {
                    if (timedOut) continue;
                    if (timeoutSec > 0 && (ii & 63) == 0 && tm.sec() > timeoutSec) {
                        timedOut = true;
                        continue;
                    }
                    int v = vs[ii];
                    int k = (int)T.bag[v].size();
                    uint64_t base = T.scBase[v], tri = T.triBase[v];
                    for (int i = 0; i < k; i++) {
                        const auto& vaDown = S[2 * (base + i) + 1];  // a -> v
                        const auto& vaUp = S[2 * (base + i)];        // v -> a
                        for (int j = i + 1; j < k; j++) {
                            uint64_t sc = T.triSc[tri++];
                            const auto& vbUp = S[2 * (base + j)];        // v -> b
                            const auto& vbDown = S[2 * (base + j) + 1];  // b -> v
                            for (int dir = 0; dir < 2; dir++) {
                                const auto& x = dir == 0 ? vaDown : vbDown;
                                const auto& y = dir == 0 ? vbUp : vaUp;
                                if (x.empty() || y.empty()) continue;
                                pool.clear();
                                joinInto(x.data(), x.size() / r, y.data(), y.size() / r, r, pool,
                                         collectStats ? &progress : nullptr);
                                uint64_t target = 2 * sc + dir;
                                lock_guard<mutex> g(locks[target % NL]);
                                auto& tg = S[target];
                                if (collectStats) progress.observe((pool.size() + tg.size()) / r, false);
                                pool.insert(pool.end(), tg.begin(), tg.end());
                                pr.run(pool);
                                tg = pool;
                            }
                        }
                    }
                }
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    pstat.add(pr.st);
                    if (collectStats) shortcutStats.add(pr.st);
                }
            }
            if (timedOut) return;
        }
        for (auto& set : S) set.shrink_to_fit();
        shortcutSec = tm.sec();
    }

    inline void entry(const Label& l, int depthW, const float*& p, size_t& cnt) const {
        uint32_t a = l.off[depthW], b = l.off[depthW + 1];
        p = l.sig.data() + a;
        cnt = (b - a) / r;
    }
    inline size_t entryCount(const Label& l, int depthW) const { return (l.off[depthW + 1] - l.off[depthW]) / r; }

    // Top-down phase, level by level so that a level runs in parallel (Algorithm 3, Lines 14-22).
    void buildLabels() {
        Timer tm;
        L[0].assign(T.n, {});
        L[1].assign(T.n, {});
        vector<vector<int>> levels(T.height + 1);
        for (int v = 0; v < T.n; v++) levels[T.depth[v]].push_back(v);
        for (int lev = 0; lev <= T.height; lev++) {
            if (timeoutSec > 0 && tm.sec() > timeoutSec) {
                timedOut = true;
                return;
            }
            const auto& vs = levels[lev];
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                Pruner pr;
                initPruner(pr);
                vector<float> pool;
                vector<int> anc;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
                for (int ii = 0; ii < (int)vs.size(); ii++) {
                    if (timedOut) continue;
                    if (timeoutSec > 0 && (ii & 63) == 0 && tm.sec() > timeoutSec) {
                        timedOut = true;
                        continue;
                    }
                    int v = vs[ii], dv = T.depth[v];
                    anc.assign(dv, 0);
                    for (int x = T.parent[v]; x >= 0; x = T.parent[x]) anc[T.depth[x]] = x;
                    const auto& bag = T.bag[v];
                    uint64_t base = T.scBase[v];
                    for (int dir = 0; dir < 2; dir++) {
                        Label& lb = L[dir][v];
                        lb.off.assign(dv + 1, 0);
                        lb.mn.assign((size_t)dv * r, 0.f);
                    }
                    for (int dw = 0; dw < dv; dw++) {
                        int w = anc[dw];
                        for (int dir = 0; dir < 2; dir++) {
                            pool.clear();
                            for (int i = 0; i < (int)bag.size(); i++) {
                                int u = bag[i];
                                const auto& sc = S[2 * (base + i) + dir];
                                if (sc.empty()) continue;
                                if (u == w) {
                                if (collectStats) progress.observe((pool.size() + sc.size()) / r, false);
                                    pool.insert(pool.end(), sc.begin(), sc.end());
                                    continue;
                                }
                                const float* p;
                                size_t cnt;
                                int du = T.depth[u];
                                if (dir == 0) {  // v->u then u->w
                                    if (du > dw) entry(L[0][u], dw, p, cnt);
                                    else entry(L[1][w], du, p, cnt);
                                } else {         // w->u then u->v
                                    if (du > dw) entry(L[1][u], dw, p, cnt);
                                    else entry(L[0][w], du, p, cnt);
                                }
                                if (cnt) joinInto(sc.data(), sc.size() / r, p, cnt, r, pool,
                                                  collectStats ? &progress : nullptr);
                            }
                            pr.run(pool);
                            Label& lb = L[dir][v];
                            lb.sig.insert(lb.sig.end(), pool.begin(), pool.end());
                            lb.off[dw + 1] = (uint32_t)lb.sig.size();
                            float* m = &lb.mn[(size_t)dw * r];
                            size_t c = pool.size() / r;
                            for (int j = 0; j < r; j++) m[j] = c ? FINF : 0.f;
                            for (size_t q = 0; q < c; q++)
                                for (int j = 0; j < r; j++) m[j] = min(m[j], pool[q * r + j]);
                        }
                    }
                    for (int dir = 0; dir < 2; dir++) L[dir][v].sig.shrink_to_fit();
                }
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    pstat.add(pr.st);
                    if (collectStats) labelStats.add(pr.st);
                }
            }
            if (timedOut) return;
        }
        labelSec = tm.sec();
    }

    void build() {
        try {
            if (collectStats) { progress.phase.store(0); progress.log("start"); }
            buildShortcuts();
            shortcutCompleted = !timedOut;
            if (collectStats) {
                if (shortcutCompleted) progress.completed.fetch_or(1);
                progress.log(shortcutCompleted ? "complete" : "TO");
            }
            if (!timedOut) {
                labelStarted = true;
                if (collectStats) { progress.phase.store(1); progress.log("start"); }
                buildLabels();
                labelCompleted = !timedOut;
                if (collectStats) {
                    if (labelCompleted) progress.completed.fetch_or(2);
                    progress.log(labelCompleted ? "complete" : "TO");
                }
            }
        } catch (...) {
            if (collectStats) progress.log("exception");
            throw;
        }
    }
    void dropShortcuts() { vector<vector<float>>().swap(S); }

    // Sizes: signatures, componentwise minima (stored only for sets with two or more signatures,
    // as in Section IV-B), and offsets.
    uint64_t sigBytes() const {
        uint64_t b = 0;
        for (int dir = 0; dir < 2; dir++)
            for (const auto& l : L[dir]) b += l.sig.size() * 4;
        return b;
    }
    uint64_t minBytes() const {
        uint64_t b = 0;
        for (int dir = 0; dir < 2; dir++)
            for (const auto& l : L[dir])
                for (size_t i = 0; i + 1 < l.off.size(); i++)
                    if ((l.off[i + 1] - l.off[i]) / r >= 2) b += (uint64_t)r * 4;
        return b;
    }
    uint64_t offBytes() const {
        uint64_t b = 0;
        for (int dir = 0; dir < 2; dir++)
            for (const auto& l : L[dir]) b += l.off.size() * 4;
        return b;
    }
    uint64_t shortcutBytes() const {
        uint64_t b = 0;
        for (const auto& s : S) b += s.size() * 4;
        return b;
    }
    SizeStats labelSizeStats() const {
        vector<uint64_t> hist(1 << 20, 0);
        for (int dir = 0; dir < 2; dir++)
            for (const auto& l : L[dir])
                for (size_t i = 0; i + 1 < l.off.size(); i++)
                    hist[min<size_t>((l.off[i + 1] - l.off[i]) / r, hist.size() - 1)]++;
        return sizeStats(hist);
    }

    void save(const string& path) const {
        FILE* f = fopen(path.c_str(), "wb");
        if (!f) throw runtime_error("cannot write " + path);
        writeAll(f, "SIGHOP02", 8);
        int64_t hdr[4] = {(int64_t)T.n, (int64_t)r, (int64_t)d.N, usePoints ? 1 : 0};
        writeAll(f, hdr, sizeof(hdr));
        double tms[2] = {shortcutSec, labelSec};
        writeAll(f, tms, sizeof(tms));
        for (int dir = 0; dir < 2; dir++) {
            for (const auto& l : L[dir]) {
                uint64_t a = l.off.size(), b = l.sig.size(), c = l.mn.size();
                writeAll(f, &a, 8);
                writeAll(f, &b, 8);
                writeAll(f, &c, 8);
                writeAll(f, l.off.data(), a * 4);
                writeAll(f, l.sig.data(), b * 4);
                writeAll(f, l.mn.data(), c * 4);
            }
        }
        fclose(f);
    }
    void load(const string& path) {
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) throw runtime_error("cannot read " + path);
        char magic[8];
        readAll(f, magic, 8);
        if (memcmp(magic, "SIGHOP02", 8) != 0) throw runtime_error("bad index file " + path);
        int64_t hdr[4];
        readAll(f, hdr, sizeof(hdr));
        if (hdr[0] != T.n || hdr[1] != r) throw runtime_error("index does not match the dataset");
        usePoints = hdr[3] != 0;
        double tms[2];
        readAll(f, tms, sizeof(tms));
        shortcutSec = tms[0];
        labelSec = tms[1];
        for (int dir = 0; dir < 2; dir++) {
            L[dir].assign(T.n, {});
            for (auto& l : L[dir]) {
                uint64_t a, b, c;
                readAll(f, &a, 8);
                readAll(f, &b, 8);
                readAll(f, &c, 8);
                l.off.resize(a);
                l.sig.resize(b);
                l.mn.resize(c);
                readAll(f, l.off.data(), a * 4);
                readAll(f, l.sig.data(), b * 4);
                readAll(f, l.mn.data(), c * 4);
            }
        }
        fclose(f);
    }

    // min over the signatures of a set, for every selected snapshot (Lines 9-10).
    static void sideMin(const float* set, size_t cnt, int r, const float* zsoa, int k, float* out, float* tmp) {
        fill(out, out + k, FINF);
        for (size_t c = 0; c < cnt; c++) {
            const float* sg = set + c * r;
            const float c0 = sg[0];
            for (int i = 0; i < k; i++) tmp[i] = c0 * zsoa[i];
            for (int j = 1; j < r; j++) {
                const float cj = sg[j];
                const float* col = zsoa + (size_t)j * k;
                for (int i = 0; i < k; i++) tmp[i] += cj * col[i];
            }
            for (int i = 0; i < k; i++) out[i] = tmp[i] < out[i] ? tmp[i] : out[i];
        }
    }

    // A TTQ query (Algorithm 2). vals receives the fastest time of each selected snapshot;
    // out[i] receives the value of ranks[i] (0-based ranks, sorted).
    void answer(int s, int t, const OmegaCtx& ctx, const vector<int>& ranks, const QueryOpts& o, float* out,
                vector<float>& vals, QueryStats* st) const {
        const int k = ctx.k;
        vals.assign(k, FINF);
        static thread_local vector<float> xa, xb, tmp, cstar, x;
        xa.resize(k);
        xb.resize(k);
        tmp.resize(k);
        cstar.resize(r);
        x.resize(r);
        const float* zs = ctx.zsoa.data();
        if (s == t) {
            fill(vals.begin(), vals.end(), 0.f);
        } else {
            const vector<int>* A;
            const vector<int>* B;
            int c = T.separatorPair(s, t, A, B);
            if (!A) {  // Lines 3-4
                const float* p;
                size_t cnt;
                if (c == t) entry(L[0][s], T.depth[t], p, cnt);
                else entry(L[1][t], T.depth[s], p, cnt);
                sideMin(p, cnt, r, zs, k, vals.data(), tmp.data());
                if (st) {
                    st->hubs += 1;
                    st->hubsKept += 1;
                    st->cands += (double)cnt;
                    st->candsKept += (double)cnt;
                }
            } else {
                const Label& ls = L[0][s];
                const Label& lt = L[1][t];
                auto cost = [&](const vector<int>& H) {
                    size_t c2 = 0;
                    for (int h : H) c2 += entryCount(ls, T.depth[h]) + entryCount(lt, T.depth[h]);
                    return c2;
                };
                const vector<int>* H;
                size_t candTotal;
                if (o.costSep) {  // Line 6
                    size_t ca = cost(*A), cb = cost(*B);
                    H = ca <= cb ? A : B;
                    candTotal = min(ca, cb);
                } else {
                    H = A->size() <= B->size() ? A : B;
                    candTotal = cost(*H);
                }
                // Line 7: reference from the representatives (first signature of each set).
                int hstar = -1;
                float bestV = FINF;
                const float* mu = ctx.mean.data();
                if (o.ffRef) {
                    static thread_local vector<float> freeFlow;
                    freeFlow.assign(r, 0.f);
                    freeFlow[0] = 1.f;
                    mu = freeFlow.data();
                }
                for (int idx = 0; idx < (int)H->size(); idx++) {
                    int dh = T.depth[(*H)[idx]];
                    if (!entryCount(ls, dh) || !entryCount(lt, dh)) continue;
                    float v = dotf(ls.sig.data() + ls.off[dh], mu, r) + dotf(lt.sig.data() + lt.off[dh], mu, r);
                    if (v < bestV) {
                        bestV = v;
                        hstar = idx;
                    }
                }
                if (hstar >= 0) {
                    int dh = T.depth[(*H)[hstar]];
                    const float* ps = ls.sig.data() + ls.off[dh];
                    const float* pt = lt.sig.data() + lt.off[dh];
                    for (int j = 0; j < r; j++) cstar[j] = ps[j] + pt[j];
                }
                double kept = 0, keptC = 0;
                for (int idx = 0; idx < (int)H->size(); idx++) {
                    int dh = T.depth[(*H)[idx]];
                    size_t ca = entryCount(ls, dh), cb = entryCount(lt, dh);
                    if (!ca || !cb) continue;
                    if (o.hubPrune && idx != hstar) {  // Line 8
                        // A set with one signature is its own componentwise minimum, so m is only
                        // needed (and only counted in minBytes) for sets with two or more.
                        const float* ms = ca == 1 ? ls.sig.data() + ls.off[dh] : &ls.mn[(size_t)dh * r];
                        const float* mt = cb == 1 ? lt.sig.data() + lt.off[dh] : &lt.mn[(size_t)dh * r];
                        for (int j = 0; j < r; j++) x[j] = ms[j] + mt[j] - cstar[j];
                        if (lowerBox(x.data(), ctx.lo.data(), ctx.hi.data(), r) >= 0) continue;
                    }
                    kept += 1;
                    keptC += (double)(ca + cb);
                    sideMin(ls.sig.data() + ls.off[dh], ca, r, zs, k, xa.data(), tmp.data());
                    sideMin(lt.sig.data() + lt.off[dh], cb, r, zs, k, xb.data(), tmp.data());
                    float* v = vals.data();
                    for (int i = 0; i < k; i++) {
                        float y = xa[i] + xb[i];
                        v[i] = y < v[i] ? y : v[i];
                    }
                }
                if (st) {
                    st->hubs += (double)H->size();
                    st->hubsKept += kept;
                    st->cands += (double)candTotal;
                    st->candsKept += keptC;
                }
            }
        }
        // Lines 11-12: multi-way selection on a copy, so that vals keeps the per-snapshot times.
        static thread_local vector<float> sel;
        sel = vals;
        multiSelect(sel.data(), k, ranks, 0, (int)ranks.size());
        for (size_t i = 0; i < ranks.size(); i++) out[i] = sel[ranks[i]];
    }
};
