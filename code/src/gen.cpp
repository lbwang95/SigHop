// Dataset generator: small test grids and DIMACS road networks with semi-synthetic traffic.
//   gen grid   --k 12 --r 4 --days 7 --slots 8 --oneway 0.1 --seed 1 --out t.bin [--raw t.raw]
//   gen dimacs --gr USA-road-t.NY.gr --co USA-road-d.NY.co --r 4 --z-from CD_r4.bin --days 90
//              --seed 1 --out ny_r4_d90.bin      (noise: GEN_DAY_NOISE=0.1 GEN_SLOT_NOISE=0.05)
// With --z-from and --days > 0, the factor vectors of the history are a day-level bootstrap of the
// source model: each synthetic day copies all slots of a random source day of the same class
// (weekday or weekend), the delay factors are scaled so that their 95th percentile over the source is
// 1, and they are perturbed by a log-normal factor per day and per snapshot. The free-flow factor stays 1.
#include "common.hpp"

struct RawGraph {
    int n = 0;
    vector<int> from, to;
    vector<float> base;
    vector<float> x, y;
};

// Keeps the largest strongly connected component (iterative Kosaraju).
static RawGraph largestScc(const RawGraph& g) {
    int n = g.n;
    vector<vector<int>> out(n), in(n);
    for (size_t e = 0; e < g.from.size(); e++) {
        out[g.from[e]].push_back(g.to[e]);
        in[g.to[e]].push_back(g.from[e]);
    }
    vector<int> finish, ptr(n, 0);
    vector<char> seen(n, 0);
    finish.reserve(n);
    for (int s = 0; s < n; s++) {
        if (seen[s]) continue;
        vector<int> st{s};
        seen[s] = 1;
        while (!st.empty()) {
            int x = st.back();
            if (ptr[x] < (int)out[x].size()) {
                int y = out[x][ptr[x]++];
                if (!seen[y]) {
                    seen[y] = 1;
                    st.push_back(y);
                }
            } else {
                finish.push_back(x);
                st.pop_back();
            }
        }
    }
    vector<int> comp(n, -1);
    int nc = 0;
    vector<int> size;
    for (int i = n - 1; i >= 0; i--) {
        int s = finish[i];
        if (comp[s] >= 0) continue;
        vector<int> st{s};
        comp[s] = nc;
        int cnt = 0;
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            cnt++;
            for (int y : in[x]) {
                if (comp[y] < 0) {
                    comp[y] = nc;
                    st.push_back(y);
                }
            }
        }
        size.push_back(cnt);
        nc++;
    }
    int best = (int)(max_element(size.begin(), size.end()) - size.begin());
    vector<int> id(n, -1);
    RawGraph r;
    for (int v = 0; v < n; v++) {
        if (comp[v] == best) {
            id[v] = r.n++;
            if (!g.x.empty()) {
                r.x.push_back(g.x[v]);
                r.y.push_back(g.y[v]);
            }
        }
    }
    for (size_t e = 0; e < g.from.size(); e++) {
        int a = id[g.from[e]], b = id[g.to[e]];
        if (a >= 0 && b >= 0) {
            r.from.push_back(a);
            r.to.push_back(b);
            r.base.push_back(g.base[e]);
        }
    }
    return r;
}

// Delay factors: each is a Gaussian congestion field over space; its intensity follows
// a day-level draw and a morning/evening time profile.
static Dataset withSyntheticTraffic(const RawGraph& g, int r, int days, int slots, int slotLo, const string& zFrom,
                                    unsigned seed) {
    mt19937 rng(seed);
    uniform_real_distribution<float> U01(0, 1);
    normal_distribution<float> N01(0, 1);
    Dataset d;
    d.n = g.n;
    d.m = (int64_t)g.from.size();
    d.r = r;
    d.from.assign(g.from.begin(), g.from.end());
    d.to.assign(g.to.begin(), g.to.end());
    d.x = g.x;
    d.y = g.y;
    float xmin = *min_element(g.x.begin(), g.x.end()), xmax = *max_element(g.x.begin(), g.x.end());
    float ymin = *min_element(g.y.begin(), g.y.end()), ymax = *max_element(g.y.begin(), g.y.end());
    float ext = max(xmax - xmin, ymax - ymin);
    struct Field {
        float cx, cy, s, gamma;
    };
    vector<Field> fields;
    for (int j = 1; j < r; j++) {
        Field f;
        f.cx = xmin + U01(rng) * (xmax - xmin);
        f.cy = ymin + U01(rng) * (ymax - ymin);
        f.s = ext * (0.1f + 0.2f * U01(rng));
        f.gamma = 0.5f + U01(rng);
        fields.push_back(f);
    }
    d.u.resize(d.m * r);
    for (int64_t e = 0; e < d.m; e++) {
        float mx = (g.x[d.from[e]] + g.x[d.to[e]]) / 2, my = (g.y[d.from[e]] + g.y[d.to[e]]) / 2;
        float b = g.base[e];
        d.u[e * r] = b;
        for (int j = 1; j < r; j++) {
            const Field& f = fields[j - 1];
            float dx = mx - f.cx, dy = my - f.cy;
            float val = exp(-(dx * dx + dy * dy) / (2 * f.s * f.s));
            d.u[e * r + j] = b * f.gamma * val * exp(0.3f * N01(rng));
        }
    }
    if (!zFrom.empty() && days > 0) {
        Dataset real = loadDataset(zFrom);
        if (real.r != r) throw runtime_error("--z-from has a different r");
        vector<float> scale(r, 1.f);
        for (int j = 1; j < r; j++) {
            vector<float> col;
            for (int64_t w = 0; w < real.N; w++) col.push_back(real.z[w * r + j]);
            size_t q = (size_t)(0.95 * (col.size() - 1));
            nth_element(col.begin(), col.begin() + q, col.end());
            scale[j] = col[q] > 0 ? 1.f / col[q] : 1.f;
        }
        map<int, vector<int>> byDay;
        for (int64_t w = 0; w < real.N; w++) byDay[real.wday[w]].push_back((int)w);
        size_t full = 0;
        for (auto& kv : byDay) full = max(full, kv.second.size());
        vector<int> wk, we;  // complete source days by class
        for (auto& kv : byDay) {
            if (kv.second.size() != full) continue;
            sort(kv.second.begin(), kv.second.end(), [&](int a, int b) { return real.wslot[a] < real.wslot[b]; });
            (real.wdow[kv.second[0]] < 5 ? wk : we).push_back(kv.first);
        }
        if (wk.empty()) throw runtime_error("--z-from has no complete weekday");
        if (we.empty()) we = wk;
        double dn = getenv("GEN_DAY_NOISE") ? atof(getenv("GEN_DAY_NOISE")) : 0.1;
        double sn = getenv("GEN_SLOT_NOISE") ? atof(getenv("GEN_SLOT_NOISE")) : 0.05;
        mt19937 zr(seed * 7919u + 17u);
        d.days = days;
        d.slotsPerDay = real.slotsPerDay;
        d.slotSeconds = real.slotSeconds;
        d.N = (int64_t)days * (int64_t)full;
        d.z.resize(d.N * r);
        d.wday.resize(d.N);
        d.wslot.resize(d.N);
        d.wdow.resize(d.N);
        int64_t w = 0;
        for (int dd = 0; dd < days; dd++) {
            int dow = dd % 7;
            const auto& pool = dow < 5 ? wk : we;
            int src = pool[zr() % pool.size()];
            vector<float> g(r, 1.f);
            for (int j = 1; j < r; j++) g[j] = exp((float)dn * N01(zr));
            for (int sw : byDay[src]) {
                d.wday[w] = dd;
                d.wslot[w] = real.wslot[sw];
                d.wdow[w] = dow;
                d.z[w * r] = 1.f;
                for (int j = 1; j < r; j++)
                    d.z[w * r + j] = max(0.f, real.z[(int64_t)sw * r + j] * scale[j] * g[j] * exp((float)sn * N01(zr)));
                w++;
            }
        }
        return d;
    }
    if (!zFrom.empty()) {
        Dataset real = loadDataset(zFrom);
        if (real.r != r) throw runtime_error("--z-from has a different r");
        d.z = real.z;
        d.N = real.N;
        d.days = real.days;
        d.slotsPerDay = real.slotsPerDay;
        d.slotSeconds = real.slotSeconds;
        d.wday = real.wday;
        d.wslot = real.wslot;
        d.wdow = real.wdow;
        return d;
    }
    d.days = days;
    d.slotsPerDay = 96;
    d.slotSeconds = 900;
    d.N = (int64_t)days * slots;
    d.z.resize(d.N * r);
    d.wday.resize(d.N);
    d.wslot.resize(d.N);
    d.wdow.resize(d.N);
    vector<float> mix(r), wid(r);
    for (int j = 1; j < r; j++) {
        mix[j] = U01(rng);
        wid[j] = 4 + 6 * U01(rng);
    }
    for (int dd = 0; dd < days; dd++) {
        vector<float> intensity(r);
        int dow = dd % 7;
        for (int j = 1; j < r; j++) intensity[j] = (0.3f + 1.2f * U01(rng)) * (dow < 5 ? 1.0f : 0.4f);
        for (int s = 0; s < slots; s++) {
            int64_t w = (int64_t)dd * slots + s;
            int slot = slotLo + s;
            d.wday[w] = dd;
            d.wslot[w] = slot;
            d.wdow[w] = dow;
            d.z[w * r] = 1;
            for (int j = 1; j < r; j++) {
                float am = (slot - 32.f) / wid[j], pm = (slot - 70.f) / wid[j];
                float g1 = mix[j] * exp(-am * am / 2) + (1 - mix[j]) * exp(-pm * pm / 2);
                d.z[w * r + j] = intensity[j] * g1 + 0.03f * U01(rng);
            }
        }
    }
    return d;
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) throw runtime_error("usage: gen grid|dimacs ...");
        string cmd = argv[1];
        Args a(argc, argv, 2);
        unsigned seed = (unsigned)a.geti("seed", 1);
        mt19937 rng(seed);
        uniform_real_distribution<float> U01(0, 1);
        RawGraph g;
        if (cmd == "grid") {
            int k = a.geti("k", 12);
            double oneway = a.getd("oneway", 0.1);
            g.n = k * k;
            g.x.resize(g.n);
            g.y.resize(g.n);
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    g.x[i * k + j] = i + 0.6f * (U01(rng) - 0.5f);
                    g.y[i * k + j] = j + 0.6f * (U01(rng) - 0.5f);
                }
            }
            auto add = [&](int p, int q) {
                float dx = g.x[p] - g.x[q], dy = g.y[p] - g.y[q];
                float b = sqrt(dx * dx + dy * dy) * (30 + 30 * U01(rng));
                float c = U01(rng);
                if (c >= oneway) {
                    g.from.push_back(p); g.to.push_back(q); g.base.push_back(b);
                    g.from.push_back(q); g.to.push_back(p); g.base.push_back(b * (0.9f + 0.2f * U01(rng)));
                } else if (c < oneway / 2) {
                    g.from.push_back(p); g.to.push_back(q); g.base.push_back(b);
                } else {
                    g.from.push_back(q); g.to.push_back(p); g.base.push_back(b);
                }
            };
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    if (i + 1 < k) add(i * k + j, (i + 1) * k + j);
                    if (j + 1 < k) add(i * k + j, i * k + j + 1);
                    if (i + 1 < k && j + 1 < k && U01(rng) < 0.2f) add(i * k + j, (i + 1) * k + j + 1);
                }
            }
            g = largestScc(g);
        } else if (cmd == "dimacs") {
            ifstream co(a.need("co"));
            if (!co) throw runtime_error("cannot read --co");
            string line;
            while (getline(co, line)) {
                if (line.rfind("p ", 0) == 0) {
                    string p, aux, sp;
                    stringstream ss(line);
                    ss >> p >> aux >> sp >> sp >> g.n;
                    g.x.resize(g.n);
                    g.y.resize(g.n);
                } else if (line.rfind("v ", 0) == 0) {
                    stringstream ss(line);
                    char c;
                    int id;
                    double px, py;
                    ss >> c >> id >> px >> py;
                    g.x[id - 1] = (float)(px / 1e6);
                    g.y[id - 1] = (float)(py / 1e6);
                }
            }
            ifstream gr(a.need("gr"));
            if (!gr) throw runtime_error("cannot read --gr");
            double unit = a.getd("unit", 0.1);  // scales the arc weights of the .gr file to seconds
            while (getline(gr, line)) {
                if (line.rfind("a ", 0) == 0) {
                    stringstream ss(line);
                    char c;
                    int u, v;
                    double w;
                    ss >> c >> u >> v >> w;
                    g.from.push_back(u - 1);
                    g.to.push_back(v - 1);
                    g.base.push_back((float)(w * unit));
                }
            }
            g = largestScc(g);
        } else {
            throw runtime_error("unknown command " + cmd);
        }
        int r = a.geti("r", 4), days = a.geti("days", 7), slots = a.geti("slots", 8);
        Dataset d = withSyntheticTraffic(g, r, days, slots, a.geti("slot-lo", 24), a.get("z-from"), seed);
        saveDataset(d, a.need("out"));
        if (a.kv.count("raw")) {
            // Reference worlds for tests: model times with noise on a fraction of arcs.
            double frac = a.getd("raw-frac", 0.2), noise = a.getd("raw-noise", 0.1);
            RawWorlds rw;
            rw.N = d.N;
            rw.off.assign(d.N + 1, 0);
            normal_distribution<float> N01(0, 1);
            for (int64_t w = 0; w < d.N; w++) {
                for (int64_t e = 0; e < d.m; e++) {
                    if (U01(rng) < frac) {
                        float t = 0;
                        for (int j = 0; j < r; j++) t += d.U(e)[j] * d.Z(w)[j];
                        rw.edge.push_back((int32_t)e);
                        rw.w.push_back(max(0.01f, t * (1 + (float)noise * N01(rng))));
                    }
                }
                rw.off[w + 1] = (int64_t)rw.edge.size();
            }
            rw.nnz = (int64_t)rw.edge.size();
            saveRaw(rw, a.get("raw"));
        }
        printf("n=%lld m=%lld r=%lld N=%lld\n", (long long)d.n, (long long)d.m, (long long)d.r, (long long)d.N);
    } catch (const bad_alloc&) {
        fprintf(stderr, "error: out of memory\n");
        return 4;
    } catch (const exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
