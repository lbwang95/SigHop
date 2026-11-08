// Competitors: per-snapshot Dijkstra, CCH-All (one customized weight array per snapshot),
// CCH-Batch (the k-metric query of Dibbelt et al.), and TTP-adapt (Li, Yiu, Mouratidis).
#pragma once
#include "sig.hpp"

struct GraphCsr {
    vector<int64_t> hf, hb;
    vector<int32_t> af, ab;  // arc ids, grouped by tail (forward) and by head (backward)
    void build(const Dataset& d) {
        hf.assign(d.n + 1, 0);
        hb.assign(d.n + 1, 0);
        for (int64_t e = 0; e < d.m; e++) {
            hf[d.from[e] + 1]++;
            hb[d.to[e] + 1]++;
        }
        for (int64_t i = 0; i < d.n; i++) {
            hf[i + 1] += hf[i];
            hb[i + 1] += hb[i];
        }
        af.resize(d.m);
        ab.resize(d.m);
        vector<int64_t> pf(hf.begin(), hf.end() - 1), pb(hb.begin(), hb.end() - 1);
        for (int64_t e = 0; e < d.m; e++) {
            af[pf[d.from[e]]++] = (int32_t)e;
            ab[pb[d.to[e]]++] = (int32_t)e;
        }
    }
};

inline double dotd(const float* a, const float* b, int r) {
    double s = 0;
    for (int j = 0; j < r; j++) s += (double)a[j] * b[j];
    return s;
}

struct Dijkstra {
    const Dataset& d;
    const GraphCsr& g;
    vector<double> df, db;
    vector<int> sf, sb, parentArc;
    int stamp = 0;
    using QE = pair<double, int>;

    Dijkstra(const Dataset& data, const GraphCsr& csr)
        : d(data), g(csr), df(data.n), db(data.n), sf(data.n, 0), sb(data.n, 0), parentArc(data.n, -1) {}

    // Bidirectional search; weight(e) returns the travel time of arc e.
    template <class W>
    double biQuery(int s, int t, W weight) {
        if (s == t) return 0;
        stamp++;
        priority_queue<QE, vector<QE>, greater<QE>> pf, pb;
        df[s] = 0; sf[s] = stamp; pf.push({0, s});
        db[t] = 0; sb[t] = stamp; pb.push({0, t});
        double best = numeric_limits<double>::infinity();
        while (!pf.empty() && !pb.empty()) {
            if (pf.top().first + pb.top().first >= best) break;
            bool forward = pf.top().first <= pb.top().first;
            if (forward) {
                auto [du, x] = pf.top();
                pf.pop();
                if (du > df[x]) continue;
                for (int64_t i = g.hf[x]; i < g.hf[x + 1]; i++) {
                    int e = g.af[i], y = d.to[e];
                    double nd = du + weight(e);
                    if (sf[y] != stamp || nd < df[y]) {
                        df[y] = nd;
                        sf[y] = stamp;
                        pf.push({nd, y});
                    }
                    if (sb[y] == stamp && df[y] + db[y] < best) best = df[y] + db[y];
                }
            } else {
                auto [du, x] = pb.top();
                pb.pop();
                if (du > db[x]) continue;
                for (int64_t i = g.hb[x]; i < g.hb[x + 1]; i++) {
                    int e = g.ab[i], y = d.from[e];
                    double nd = du + weight(e);
                    if (sb[y] != stamp || nd < db[y]) {
                        db[y] = nd;
                        sb[y] = stamp;
                        pb.push({nd, y});
                    }
                    if (sf[y] == stamp && df[y] + db[y] < best) best = df[y] + db[y];
                }
            }
        }
        return best;
    }

    // One-to-all search; the distances stay in df for vertices with sf == stamp.
    template <class W>
    void oneToAll(int s, W weight) {
        stamp++;
        priority_queue<QE, vector<QE>, greater<QE>> pq;
        df[s] = 0; sf[s] = stamp; pq.push({0, s});
        while (!pq.empty()) {
            auto [du, x] = pq.top();
            pq.pop();
            if (du > df[x]) continue;
            for (int64_t i = g.hf[x]; i < g.hf[x + 1]; i++) {
                int e = g.af[i], y = d.to[e];
                double nd = du + weight(e);
                if (sf[y] != stamp || nd < df[y]) {
                    df[y] = nd;
                    sf[y] = stamp;
                    parentArc[y] = e;
                    pq.push({nd, y});
                }
            }
        }
    }

    // Unidirectional search that also returns the arcs of the path.
    template <class W>
    double pathQuery(int s, int t, W weight, vector<int>& arcs) {
        arcs.clear();
        stamp++;
        priority_queue<QE, vector<QE>, greater<QE>> pq;
        df[s] = 0; sf[s] = stamp; parentArc[s] = -1; pq.push({0, s});
        while (!pq.empty()) {
            auto [du, x] = pq.top();
            pq.pop();
            if (du > df[x]) continue;
            if (x == t) break;
            for (int64_t i = g.hf[x]; i < g.hf[x + 1]; i++) {
                int e = g.af[i], y = d.to[e];
                double nd = du + weight(e);
                if (sf[y] != stamp || nd < df[y]) {
                    df[y] = nd;
                    sf[y] = stamp;
                    parentArc[y] = e;
                    pq.push({nd, y});
                }
            }
        }
        if (sf[t] != stamp) return numeric_limits<double>::infinity();
        for (int x = t; x != s; x = d.from[parentArc[x]]) arcs.push_back(parentArc[x]);
        return df[t];
    }

    // Time-dependent search: an arc entered at time x uses the world of day/slot at x.
    // worldAt(day, slot) clamps to the nearest available slot of that day.
    template <class WorldAt>
    double tdQuery(int s, int t, int day, int startSlot, WorldAt worldAt) {
        stamp++;
        priority_queue<QE, vector<QE>, greater<QE>> pq;
        df[s] = 0; sf[s] = stamp; pq.push({0, s});
        int r = (int)d.r;
        while (!pq.empty()) {
            auto [du, x] = pq.top();
            pq.pop();
            if (du > df[x]) continue;
            if (x == t) return du;
            int slot = startSlot + (int)floor(du / d.slotSeconds);
            const float* z = d.Z(worldAt(day, slot));
            for (int64_t i = g.hf[x]; i < g.hf[x + 1]; i++) {
                int e = g.af[i], y = d.to[e];
                double nd = du + dotd(d.U(e), z, r);
                if (sf[y] != stamp || nd < df[y]) {
                    df[y] = nd;
                    sf[y] = stamp;
                    pq.push({nd, y});
                }
            }
        }
        return numeric_limits<double>::infinity();
    }
};

// Customizable contraction hierarchy over the tree decomposition: one weight array per world.
struct CCH {
    const Dataset& d;
    const TreeDecomp& T;
    vector<float> fs, bs;
    CCH(const Dataset& data, const TreeDecomp& td)
        : d(data), T(td), fs(td.height + 1), bs(td.height + 1) {}

    void customize(const float* z, vector<float>& w) const {
        w.assign(2 * T.numShortcuts, FINF);
        int r = (int)d.r;
        for (int64_t e = 0; e < d.m; e++) {
            uint64_t p = T.arcPlace[e];
            if (p == UINT64_MAX) continue;
            float c = dotf(d.U(e), z, r);
            if (c < w[p]) w[p] = c;
        }
        for (int idx = 0; idx < T.n; idx++) {
            int v = T.order[idx];
            int k = (int)T.bag[v].size();
            uint64_t base = T.scBase[v], tri = T.triBase[v];
            for (int i = 0; i < k; i++) {
                float toI = w[2 * (base + i) + 1], fromI = w[2 * (base + i)];
                for (int j = i + 1; j < k; j++) {
                    uint64_t sc = T.triSc[tri++];
                    float toJ = w[2 * (base + j) + 1], fromJ = w[2 * (base + j)];
                    float c1 = toI + fromJ, c2 = toJ + fromI;
                    if (c1 < w[2 * sc]) w[2 * sc] = c1;
                    if (c2 < w[2 * sc + 1]) w[2 * sc + 1] = c2;
                }
            }
        }
    }

    float query(int s, int t, const vector<float>& w) {
        if (s == t) return 0;
        int c = T.lca(s, t), dc = T.depth[c];
        fill(fs.begin(), fs.begin() + T.depth[s] + 1, FINF);
        fill(bs.begin(), bs.begin() + T.depth[t] + 1, FINF);
        fs[T.depth[s]] = 0;
        bs[T.depth[t]] = 0;
        for (int x = s; x >= 0; x = T.parent[x]) {
            float dx = fs[T.depth[x]];
            if (dx == FINF) continue;
            uint64_t base = T.scBase[x];
            const auto& bag = T.bag[x];
            for (int i = 0; i < (int)bag.size(); i++) {
                float nd = dx + w[2 * (base + i)];
                int dy = T.depth[bag[i]];
                if (nd < fs[dy]) fs[dy] = nd;
            }
        }
        for (int x = t; x >= 0; x = T.parent[x]) {
            float dx = bs[T.depth[x]];
            if (dx == FINF) continue;
            uint64_t base = T.scBase[x];
            const auto& bag = T.bag[x];
            for (int i = 0; i < (int)bag.size(); i++) {
                float nd = dx + w[2 * (base + i) + 1];
                int dy = T.depth[bag[i]];
                if (nd < bs[dy]) bs[dy] = nd;
            }
        }
        float best = FINF;
        for (int dd = 0; dd <= dc; dd++) {
            float v = fs[dd] + bs[dd];
            if (v < best) best = v;
        }
        return best;
    }
};

// CCH-Batch: the k-metric variant of CCH (Dibbelt et al., Sec. 7.6). Each shortcut keeps the
// weights of the k selected snapshots contiguously, and one elimination-tree search computes all
// k distances with vector operations.
struct CCHBatch {
    const Dataset& d;
    const TreeDecomp& T;
    int k = 0;
    vector<float> W;  // W[(2*sc+dir)*k + i]
    vector<float> fs, bs;
    vector<char> rf, rb;
    CCHBatch(const Dataset& data, const TreeDecomp& td) : d(data), T(td) {}

    void prepare(const vector<vector<float>>& perWorld) {
        k = (int)perWorld.size();
        uint64_t P = 2 * T.numShortcuts;
        W.assign(P * (uint64_t)k, FINF);
        for (int i = 0; i < k; i++)
            for (uint64_t p = 0; p < P; p++) W[p * k + i] = perWorld[i][p];
        fs.assign((size_t)(T.height + 1) * k, FINF);
        bs.assign((size_t)(T.height + 1) * k, FINF);
        rf.assign(T.height + 1, 0);
        rb.assign(T.height + 1, 0);
    }
    uint64_t bytes() const { return W.size() * 4; }

    void query(int s, int t, float* vals) {
        if (s == t) {
            fill(vals, vals + k, 0.f);
            return;
        }
        int c = T.lca(s, t), dc = T.depth[c];
        auto sweep = [&](int src, vector<float>& D, vector<char>& reached, int dir) {
            int ds = T.depth[src];
            fill(D.begin(), D.begin() + (size_t)(ds + 1) * k, FINF);
            fill(reached.begin(), reached.begin() + ds + 1, 0);
            fill(D.begin() + (size_t)ds * k, D.begin() + (size_t)(ds + 1) * k, 0.f);
            reached[ds] = 1;
            for (int x = src; x >= 0; x = T.parent[x]) {
                int dx = T.depth[x];
                if (!reached[dx]) continue;
                const float* rx = &D[(size_t)dx * k];
                uint64_t base = T.scBase[x];
                const auto& bag = T.bag[x];
                for (int i = 0; i < (int)bag.size(); i++) {
                    int dy = T.depth[bag[i]];
                    float* ry = &D[(size_t)dy * k];
                    const float* w = &W[(2 * (base + i) + dir) * (uint64_t)k];
                    for (int q = 0; q < k; q++) {
                        float nd = rx[q] + w[q];
                        ry[q] = nd < ry[q] ? nd : ry[q];
                    }
                    reached[dy] = 1;
                }
            }
        };
        sweep(s, fs, rf, 0);
        sweep(t, bs, rb, 1);
        fill(vals, vals + k, FINF);
        for (int dd = 0; dd <= dc; dd++) {
            if (!rf[dd] || !rb[dd]) continue;
            const float* a = &fs[(size_t)dd * k];
            const float* b = &bs[(size_t)dd * k];
            for (int q = 0; q < k; q++) {
                float v = a[q] + b[q];
                vals[q] = v < vals[q] ? v : vals[q];
            }
        }
    }
};

// TTP-adapt: Phase I of the k traffic-tolerant paths algorithm of Li, Yiu, and Mouratidis
// (SIGSPATIAL'14, GeoInformatica'17) with the selected snapshots as the m time instants.
//  (1) one backward one-to-all search from t per instant gives tau_j(sp_last_j) for the bounds;
//  (2) the pruning set D holds sp_avg and kTTP paths found with w_s(e) = min_j w_j(e), raising
//      w_s(e) to max_j w_j(e) on each found path;
//  (3) a DFS enumerates simple prefixes from s and discards a prefix whose lower-bound vector is
//      dominated by a path in D (Pruning Rule 2); every completed path is a candidate.
// The fastest time of instant j is the minimum over the candidates and D. Exact for TTQ.
struct TTPAdapt {
    const Dataset& d;
    const GraphCsr& g;
    int kTTP = 3;
    double timeoutSec = 60;
    vector<float> distT;  // k * n
    vector<float> zsel;   // k * r
    long long expanded = 0;
    bool timedOut = false;
    TTPAdapt(const Dataset& data, const GraphCsr& csr) : d(data), g(csr) {}

    inline float wt(int e, int j) const { return dotf(d.U(e), &zsel[(size_t)j * d.r], (int)d.r); }

    // Returns false on timeout. vals receives the k fastest times.
    bool query(int s, int t, const vector<int>& worlds, float* vals) {
        Timer tm;
        const int k = (int)worlds.size(), n = (int)d.n, r = (int)d.r;
        zsel.resize((size_t)k * r);
        for (int j = 0; j < k; j++) copy(d.Z(worlds[j]), d.Z(worlds[j]) + r, &zsel[(size_t)j * r]);
        timedOut = false;
        expanded = 0;
        fill(vals, vals + k, FINF);
        if (s == t) {
            fill(vals, vals + k, 0.f);
            return true;
        }
        // (1) backward searches
        distT.assign((size_t)k * n, FINF);
        using QE = pair<float, int>;
        for (int j = 0; j < k; j++) {
            if (tm.sec() > timeoutSec) return !(timedOut = true);
            float* D = &distT[(size_t)j * n];
            priority_queue<QE, vector<QE>, greater<QE>> pq;
            D[t] = 0;
            pq.push({0.f, t});
            while (!pq.empty()) {
                auto [du, x] = pq.top();
                pq.pop();
                if (du > D[x]) continue;
                for (int64_t i = g.hb[x]; i < g.hb[x + 1]; i++) {
                    int e = g.ab[i], y = d.from[e];
                    float nd = du + wt(e, j);
                    if (nd < D[y]) {
                        D[y] = nd;
                        pq.push({nd, y});
                    }
                }
            }
        }
        // (2) pruning set D
        vector<float> Dvec;  // |D| * k
        auto pathVec = [&](const vector<int>& arcs) {
            size_t base = Dvec.size();
            Dvec.resize(base + k, 0.f);
            for (int e : arcs)
                for (int j = 0; j < k; j++) Dvec[base + j] += wt(e, j);
            for (int j = 0; j < k; j++) vals[j] = min(vals[j], Dvec[base + j]);
        };
        auto spath = [&](auto weight, vector<int>& arcs) {
            vector<float> D(n, FINF);
            vector<int> par(n, -1);
            priority_queue<QE, vector<QE>, greater<QE>> pq;
            D[s] = 0;
            pq.push({0.f, s});
            while (!pq.empty()) {
                auto [du, x] = pq.top();
                pq.pop();
                if (du > D[x]) continue;
                if (x == t) break;
                for (int64_t i = g.hf[x]; i < g.hf[x + 1]; i++) {
                    int e = g.af[i], y = d.to[e];
                    float nd = du + weight(e);
                    if (nd < D[y]) {
                        D[y] = nd;
                        par[y] = e;
                        pq.push({nd, y});
                    }
                }
            }
            arcs.clear();
            for (int x = t; x != s; x = d.from[par[x]]) arcs.push_back(par[x]);
        };
        vector<float> zavg(r, 0.f);
        for (int j = 0; j < k; j++)
            for (int q = 0; q < r; q++) zavg[q] += zsel[(size_t)j * r + q] / k;
        vector<int> arcs;
        spath([&](int e) { return dotf(d.U(e), zavg.data(), r); }, arcs);
        pathVec(arcs);
        vector<float> ws(d.m, -1.f);
        auto wmin = [&](int e) {
            if (ws[e] < 0) {
                float m = FINF;
                for (int j = 0; j < k; j++) m = min(m, wt(e, j));
                ws[e] = m;
            }
            return ws[e];
        };
        for (int it = 0; it < kTTP; it++) {
            if (tm.sec() > timeoutSec) return !(timedOut = true);
            spath(wmin, arcs);
            pathVec(arcs);
            for (int e : arcs) {
                float m = 0;
                for (int j = 0; j < k; j++) m = max(m, wt(e, j));
                ws[e] = m;
            }
        }
        const int nd = (int)(Dvec.size() / k);
        // (3) DFS over simple prefixes with Pruning Rule 2
        vector<char> onPath(n, 0);
        struct Frame {
            int v;
            int64_t it;
        };
        vector<Frame> st;
        vector<float> tau;  // tau[depth * k ..]
        vector<float> lb(k);
        st.push_back({s, g.hf[s]});
        tau.assign(k, 0.f);
        onPath[s] = 1;
        while (!st.empty()) {
            Frame& f = st.back();
            if (f.it == g.hf[f.v + 1]) {
                onPath[f.v] = 0;
                st.pop_back();
                tau.resize(st.size() * k);
                continue;
            }
            int e = g.af[f.it++], y = d.to[e];
            if (onPath[y]) continue;
            if ((++expanded & 1023) == 0 && tm.sec() > timeoutSec) return !(timedOut = true);
            size_t base = (st.size() - 1) * k;
            for (int j = 0; j < k; j++) lb[j] = tau[base + j] + wt(e, j);
            // complete path
            if (y == t) {
                for (int j = 0; j < k; j++) vals[j] = min(vals[j], lb[j]);
                continue;
            }
            bool pruned = false;
            for (int q = 0; q < nd && !pruned; q++) {
                const float* p = &Dvec[(size_t)q * k];
                bool dom = true;
                for (int j = 0; j < k && dom; j++) dom = p[j] <= lb[j] + distT[(size_t)j * n + y];
                pruned = dom;
            }
            if (pruned) continue;
            tau.insert(tau.end(), lb.begin(), lb.end());
            onPath[y] = 1;
            st.push_back({y, g.hf[y]});
        }
        return true;
    }
};
