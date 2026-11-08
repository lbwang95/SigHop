// Competitors: per-world Dijkstra, CCH, H2H, signature CH, typical snapshots, and the
// independent-normal reliable path estimator. All share the tree decomposition of SigHop.
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
            for (int64_t i = g.hf[x]; i < g.hf[x + 1]; i++) {
                int e = g.af[i], y = d.to[e];
                // Integrate progress along the edge under piecewise-constant speeds.
                // The fraction traversed per second is 1 / w_e(slot). Positive
                // speeds make arrival time nondecreasing in departure time (FIFO).
                double remaining = 1, nd = du;
                while (remaining > 0) {
                    int slot = startSlot + (int)floor(nd / d.slotSeconds);
                    int w = worldAt(day, slot);
                    if (w < 0) { nd = numeric_limits<double>::infinity(); break; }
                    double cost = dotd(d.U(e), d.Z(w), r);
                    if (!(cost > 0) || !isfinite(cost)) throw runtime_error("invalid time-dependent edge time");
                    double boundary = (slot - startSlot + 1) * d.slotSeconds;
                    double available = boundary - nd, needed = remaining * cost;
                    if (needed <= available) { nd += needed; break; }
                    remaining -= available / cost;
                    nd = boundary;
                }
                if (!isfinite(nd)) continue;
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

// Scalar 2-hop labels of one world (H2H), built from the customized CCH weights.
struct H2H {
    const TreeDecomp& T;
    vector<vector<float>> up, dn;
    explicit H2H(const TreeDecomp& td) : T(td) {}

    void build(const vector<float>& w) {
        up.assign(T.n, {});
        dn.assign(T.n, {});
        vector<int> anc;
        for (int idx = T.n - 1; idx >= 0; idx--) {
            int v = T.order[idx], dv = T.depth[v];
            anc.assign(dv, 0);
            for (int x = T.parent[v]; x >= 0; x = T.parent[x]) anc[T.depth[x]] = x;
            up[v].assign(dv, FINF);
            dn[v].assign(dv, FINF);
            const auto& bag = T.bag[v];
            uint64_t base = T.scBase[v];
            for (int dw = 0; dw < dv; dw++) {
                int w2 = anc[dw];
                float bu = FINF, bd = FINF;
                for (int i = 0; i < (int)bag.size(); i++) {
                    int u = bag[i], du = T.depth[u];
                    float wu = w[2 * (base + i)], wd = w[2 * (base + i) + 1];
                    float a, b;
                    if (u == w2) {
                        a = b = 0;
                    } else if (du > dw) {
                        a = up[u][dw];
                        b = dn[u][dw];
                    } else {
                        a = dn[w2][du];
                        b = up[w2][du];
                    }
                    bu = min(bu, wu + a);
                    bd = min(bd, b + wd);
                }
                up[v][dw] = bu;
                dn[v][dw] = bd;
            }
        }
    }

    float query(int s, int t) const {
        if (s == t) return 0;
        const vector<int>* sep;
        int c = T.separator(s, t, sep);
        if (!sep) return c == t ? up[s][T.depth[t]] : dn[t][T.depth[s]];
        float best = FINF;
        const float* us = up[s].data();
        const float* dt = dn[t].data();
        for (int h : *sep) {
            float v = us[T.depth[h]] + dt[T.depth[h]];
            if (v < best) best = v;
        }
        return best;
    }

    uint64_t bytes() const {
        uint64_t b = 0;
        for (int v = 0; v < T.n; v++) b += (up[v].size() + dn[v].size()) * 4;
        return b;
    }
};

// Per-world search over the signature shortcuts of SigHop (personalized-CH style).
struct SigCH {
    const SigIndex& I;
    vector<float> fs, bs;
    explicit SigCH(const SigIndex& idx) : I(idx), fs(idx.T.height + 1), bs(idx.T.height + 1) {}

    float query(int s, int t, const float* z) {
        const TreeDecomp& T = I.T;
        int r = I.r;
        if (s == t) return 0;
        int c = T.lca(s, t), dc = T.depth[c];
        fill(fs.begin(), fs.begin() + T.depth[s] + 1, FINF);
        fill(bs.begin(), bs.begin() + T.depth[t] + 1, FINF);
        fs[T.depth[s]] = 0;
        bs[T.depth[t]] = 0;
        for (int dir = 0; dir < 2; dir++) {
            vector<float>& ds = dir == 0 ? fs : bs;
            for (int x = dir == 0 ? s : t; x >= 0; x = T.parent[x]) {
                float dx = ds[T.depth[x]];
                if (dx == FINF) continue;
                uint64_t base = T.scBase[x];
                const auto& bag = T.bag[x];
                for (int i = 0; i < (int)bag.size(); i++) {
                    const auto& set = I.S[2 * (base + i) + dir];
                    float w = FINF;
                    for (size_t q = 0; q < set.size() / r; q++) w = min(w, dotf(&set[q * r], z, r));
                    float nd = dx + w;
                    int dy = T.depth[bag[i]];
                    if (nd < ds[dy]) ds[dy] = nd;
                }
            }
        }
        float best = FINF;
        for (int dd = 0; dd <= dc; dd++) best = min(best, fs[dd] + bs[dd]);
        return best;
    }
};

// Typical snapshots: k-means over the factor vectors, one CCH metric per cluster.
struct TypicalSnapshots {
    int K = 0, r = 0;
    vector<float> centers;
    vector<int> assign;
    vector<vector<float>> metric;

    void fit(const Dataset& d, const CCH& cch, int k, int iters, unsigned seed) {
        K = (int)min<int64_t>(k, d.N);
        r = (int)d.r;
        mt19937 rng(seed);
        centers.resize((size_t)K * r);
        assign.assign(d.N, 0);
        // k-means++ seeding
        vector<double> dist2(d.N, numeric_limits<double>::infinity());
        int first = uniform_int_distribution<int>(0, (int)d.N - 1)(rng);
        copy(d.Z(first), d.Z(first) + r, &centers[0]);
        for (int c = 0; c < K; c++) {
            if (c) {
                double total = accumulate(dist2.begin(), dist2.end(), 0.0);
                int w;
                if (total > 0) {
                    discrete_distribution<int> pick(dist2.begin(), dist2.end());
                    w = pick(rng);
                } else {
                    w = uniform_int_distribution<int>(0, (int)d.N - 1)(rng);
                }
                copy(d.Z(w), d.Z(w) + r, &centers[(size_t)c * r]);
            }
            for (int64_t w = 0; w < d.N; w++) {
                double s = 0;
                for (int j = 0; j < r; j++) {
                    double x = d.Z(w)[j] - centers[(size_t)c * r + j];
                    s += x * x;
                }
                dist2[w] = min(dist2[w], s);
            }
        }
        for (int it = 0; it < iters; it++) {
            vector<double> sum((size_t)K * r, 0);
            vector<int> cnt(K, 0);
            for (int64_t w = 0; w < d.N; w++) {
                double best = numeric_limits<double>::infinity();
                int arg = 0;
                for (int c = 0; c < K; c++) {
                    double s = 0;
                    for (int j = 0; j < r; j++) {
                        double x = d.Z(w)[j] - centers[(size_t)c * r + j];
                        s += x * x;
                    }
                    if (s < best) {
                        best = s;
                        arg = c;
                    }
                }
                assign[w] = arg;
                cnt[arg]++;
                for (int j = 0; j < r; j++) sum[(size_t)arg * r + j] += d.Z(w)[j];
            }
            for (int c = 0; c < K; c++)
                if (cnt[c])
                    for (int j = 0; j < r; j++) centers[(size_t)c * r + j] = (float)(sum[(size_t)c * r + j] / cnt[c]);
        }
        metric.assign(K, {});
        for (int c = 0; c < K; c++) cch.customize(&centers[(size_t)c * r], metric[c]);
    }
    uint64_t bytes() const {
        uint64_t b = centers.size() * 4 + assign.size() * 4;
        for (const auto& w : metric) b += w.size() * 4;
        return b;
    }
};

// Moments of the factors over Omega and the resulting per-arc moments under independence.
struct IndependentModel {
    vector<double> mu, var;  // per arc
    double lamScale = 1;
    void fit(const Dataset& d, const vector<int>& worlds) {
        int r = (int)d.r;
        vector<double> mean(r, 0), cov((size_t)r * r, 0);
        for (int w : worlds)
            for (int j = 0; j < r; j++) mean[j] += d.Z(w)[j];
        for (auto& x : mean) x /= max<size_t>(1, worlds.size());
        for (int w : worlds) {
            const float* z = d.Z(w);
            for (int a = 0; a < r; a++)
                for (int b = 0; b < r; b++) cov[(size_t)a * r + b] += (z[a] - mean[a]) * (z[b] - mean[b]);
        }
        for (auto& x : cov) x /= max<size_t>(1, worlds.size());
        mu.resize(d.m);
        var.resize(d.m);
        double sm = 0, sv = 0;
        for (int64_t e = 0; e < d.m; e++) {
            const float* u = d.U(e);
            double m = 0, v = 0;
            for (int a = 0; a < r; a++) {
                m += u[a] * mean[a];
                for (int b = 0; b < r; b++) v += u[a] * cov[(size_t)a * r + b] * u[b];
            }
            mu[e] = m;
            var[e] = max(0.0, v);
            sm += m;
            sv += var[e];
        }
        lamScale = sv > 0 ? sm / sv : 1;
    }
};
