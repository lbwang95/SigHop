// SigHop: signature labels on a tree decomposition, pruned by world-region dominance.
// Sections referred to below are those of the paper (htq-paper/paper/algorithm.tex).
#pragma once
#include "td.hpp"

enum RegionMode { REGION_ORTH, REGION_AABB, REGION_OBB };

inline RegionMode parseRegion(const string& s) {
    if (s == "orth") return REGION_ORTH;
    if (s == "aabb") return REGION_AABB;
    if (s == "obb") return REGION_OBB;
    throw runtime_error("unknown region " + s);
}

// Eigenvectors of a symmetric matrix by cyclic Jacobi. Columns of V are eigenvectors.
inline void jacobiEigen(vector<double> A, int r, vector<double>& V) {
    V.assign(r * r, 0.0);
    for (int i = 0; i < r; i++) V[i * r + i] = 1;
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0;
        for (int p = 0; p < r; p++)
            for (int q = p + 1; q < r; q++) off += A[p * r + q] * A[p * r + q];
        if (off < 1e-24) break;
        for (int p = 0; p < r; p++) {
            for (int q = p + 1; q < r; q++) {
                double apq = A[p * r + q];
                if (fabs(apq) < 1e-300) continue;
                double theta = (A[q * r + q] - A[p * r + p]) / (2 * apq);
                double t = (theta >= 0 ? 1 : -1) / (fabs(theta) + sqrt(theta * theta + 1));
                double c = 1 / sqrt(t * t + 1), s = t * c;
                for (int k = 0; k < r; k++) {
                    double akp = A[k * r + p], akq = A[k * r + q];
                    A[k * r + p] = c * akp - s * akq;
                    A[k * r + q] = s * akp + c * akq;
                }
                for (int k = 0; k < r; k++) {
                    double apk = A[p * r + k], aqk = A[q * r + k];
                    A[p * r + k] = c * apk - s * aqk;
                    A[q * r + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < r; k++) {
                    double vkp = V[k * r + p], vkq = V[k * r + q];
                    V[k * r + p] = c * vkp - s * vkq;
                    V[k * r + q] = s * vkp + c * vkq;
                }
            }
        }
    }
}

// Orthogonal change of basis on all loadings and factors: inner products, hence travel
// times, are unchanged, and an axis-aligned box in the new basis is an oriented box.
inline void rotateDataset(Dataset& d, const vector<int>& train) {
    int r = (int)d.r;
    vector<double> mean(r, 0), cov(r * r, 0);
    for (int w : train)
        for (int j = 0; j < r; j++) mean[j] += d.Z(w)[j];
    for (auto& x : mean) x /= max<size_t>(1, train.size());
    for (int w : train) {
        const float* z = d.Z(w);
        for (int a = 0; a < r; a++)
            for (int b = 0; b < r; b++) cov[a * r + b] += (z[a] - mean[a]) * (z[b] - mean[b]);
    }
    for (auto& x : cov) x /= max<size_t>(1, train.size());
    vector<double> V;
    jacobiEigen(cov, r, V);
    auto apply = [&](float* row) {
        vector<double> out(r, 0);
        for (int a = 0; a < r; a++)
            for (int b = 0; b < r; b++) out[a] += V[b * r + a] * row[b];
        for (int a = 0; a < r; a++) row[a] = (float)out[a];
    };
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int64_t e = 0; e < d.m; e++) apply(&d.u[e * r]);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int64_t w = 0; w < d.N; w++) apply(&d.z[w * r]);
}

// Axis-aligned box in the (possibly rotated) basis; REGION_ORTH means the whole orthant.
struct Region {
    int r = 0;
    RegionMode mode = REGION_AABB;
    vector<float> lo, hi, center;

    void fit(const Dataset& d, const vector<int>& train, double margin, RegionMode m) {
        r = (int)d.r;
        mode = m;
        lo.assign(r, FINF);
        hi.assign(r, -FINF);
        center.assign(r, 0);
        vector<double> mean(r, 0);
        for (int w : train) {
            const float* z = d.Z(w);
            for (int j = 0; j < r; j++) {
                lo[j] = min(lo[j], z[j]);
                hi[j] = max(hi[j], z[j]);
                mean[j] += z[j];
            }
        }
        for (int j = 0; j < r; j++) {
            float w = hi[j] - lo[j];
            lo[j] -= (float)margin * w;
            hi[j] += (float)margin * w;
            center[j] = (float)(mean[j] / max<size_t>(1, train.size()));
        }
    }
    bool contains(const float* z) const {
        for (int j = 0; j < r; j++)
            if (z[j] < lo[j] - 1e-6f || z[j] > hi[j] + 1e-6f) return false;
        return true;
    }
};

// a dominates b over the box [lo, hi]: the support function of a-b is non-positive.
inline bool dominatesBox(const float* a, const float* b, const float* lo, const float* hi, int r) {
    float s = 0;
    for (int j = 0; j < r; j++) {
        float d = a[j] - b[j];
        s += d > 0 ? hi[j] * d : lo[j] * d;
    }
    return s <= 0;
}
inline bool dominatesOrth(const float* a, const float* b, int r) {
    for (int j = 0; j < r; j++)
        if (a[j] > b[j]) return false;
    return true;
}
inline float dotf(const float* a, const float* b, int r) {
    float s = 0;
    for (int j = 0; j < r; j++) s += a[j] * b[j];
    return s;
}

// SigPrune of the paper: keeps a set of non-dominated signatures in place.
struct Pruner {
    int r = 0;
    RegionMode mode = REGION_AABB;
    const float* lo = nullptr;
    const float* hi = nullptr;
    const float* center = nullptr;
    const Dataset* worldFilter = nullptr;  // when set, also keep only per-world minimizers
    vector<pair<float, size_t>> ord;
    vector<float> kept;

    void run(vector<float>& pool) {
        if (r <= 0 || pool.size() % (size_t)r) throw runtime_error("invalid signature array shape");
        size_t cnt = pool.size() / (size_t)r;
        if (cnt <= 1) return;
#ifdef HTQ_TRACE_LARGE_POOLS
        bool largePool = pool.size() > (size_t)numeric_limits<int>::max();
        if (largePool) fprintf(stderr, "large candidate pool: %zu signatures, %zu floats\n", cnt, pool.size());
#endif
        ord.clear();
        if (mode == REGION_ORTH) {
            for (size_t i = 0; i < cnt; i++) {
                float s = 0;
                for (int j = 0; j < r; j++) s += pool[i * r + j];
                ord.push_back({s, i});
            }
        } else {
            // Interval filter: lb(c) > min ub means some other signature is always faster.
            float minUb = FINF;
            for (size_t i = 0; i < cnt; i++) {
                const float* c = &pool[i * r];
                float ub = 0;
                for (int j = 0; j < r; j++) ub += c[j] > 0 ? hi[j] * c[j] : lo[j] * c[j];
                minUb = min(minUb, ub);
            }
#ifdef HTQ_TRACE_LARGE_POOLS
            if (largePool) fprintf(stderr, "large-pool interval scan completed\n");
#endif
            for (size_t i = 0; i < cnt; i++) {
                const float* c = &pool[i * r];
                float lb = 0;
                for (int j = 0; j < r; j++) lb += c[j] > 0 ? lo[j] * c[j] : hi[j] * c[j];
                if (lb <= minUb) ord.push_back({dotf(c, center, r), i});
            }
        }
        sort(ord.begin(), ord.end());
        kept.clear();
        size_t nk = 0;
        for (auto& pr : ord) {
            const float* c = &pool[pr.second * r];
            bool dom = false;
            for (size_t k = 0; k < nk && !dom; k++) {
                const float* a = &kept[k * r];
                dom = mode == REGION_ORTH ? dominatesOrth(a, c, r) : dominatesBox(a, c, lo, hi, r);
            }
            if (!dom) {
                kept.insert(kept.end(), c, c + r);
                nk++;
            }
        }
        if (worldFilter && nk > 1) {
            const Dataset& d = *worldFilter;
            vector<char> keep(nk, 0);
            for (int64_t w = 0; w < d.N; w++) {
                const float* z = d.Z(w);
                float best = FINF;
                size_t arg = 0;
                for (size_t k = 0; k < nk; k++) {
                    float v = dotf(&kept[k * r], z, r);
                    if (v < best) {
                        best = v;
                        arg = k;
                    }
                }
                keep[arg] = 1;
            }
            size_t o = 0;
            for (size_t k = 0; k < nk; k++) {
                if (keep[k]) {
                    if (o != k) copy(&kept[k * r], &kept[k * r] + r, &kept[o * r]);
                    o++;
                }
            }
            kept.resize(o * r);
        }
        pool = kept;
    }
};

inline void joinInto(const float* a, size_t na, const float* b, size_t nb, int r, vector<float>& pool) {
    size_t base = pool.size();
    if (r <= 0 || (na && nb > (pool.max_size() - base) / na / (size_t)r))
        throw length_error("index size limit: signature join output exceeds address space");
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
inline SizeStats sizeStats(const vector<uint64_t>& hist, const map<size_t, uint64_t>& overflow = {}) {
    SizeStats st;
    for (size_t i = 0; i < hist.size(); i++) {
        st.count += hist[i];
        st.sum += (double)i * hist[i];
        if (hist[i]) st.maxv = i;
    }
    for (const auto& item : overflow) {
        st.count += item.second; st.sum += (double)item.first * item.second;
        if (item.second) st.maxv = item.first;
    }
    st.mean = st.count ? st.sum / st.count : 0;
    size_t acc = 0;
    bool b50 = false, b90 = false, b99 = false;
    auto accumulateBin = [&](size_t i, uint64_t count) {
        acc += count;
        if (!b50 && acc * 2 >= st.count) { st.p50 = i; b50 = true; }
        if (!b90 && acc * 10 >= st.count * 9) { st.p90 = i; b90 = true; }
        if (!b99 && acc * 100 >= st.count * 99) { st.p99 = i; b99 = true; }
    };
    for (size_t i = 0; i < hist.size(); i++) accumulateBin(i, hist[i]);
    for (const auto& item : overflow) accumulateBin(item.first, item.second);
    return st;
}

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
                float v = d.z[idx[i] * d.r + j];
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

struct HubCand {
    vector<float> A, B;  // surviving signatures; an empty side with zeroSide set contributes 0
    bool zeroA = false, zeroB = false;
};

struct QueryStats {
    long long hubs = 0, hubsKept = 0, cands = 0, candsKept = 0;
};

struct SigIndex {
    const Dataset& d;
    const TreeDecomp& T;
    Region R;
    int r;
    vector<vector<float>> S;  // S[2*sc+dir], dir 0: lower->upper, 1: upper->lower
    struct Label {
        vector<uint32_t> off;  // off[depth(w)] .. off[depth(w)+1] indexes sig, in floats
        vector<float> sig;
    };
    vector<Label> L[2];        // L[0][v]: v->ancestor, L[1][v]: ancestor->v
    bool useWorldFilter = false;
    double shortcutSec = 0, labelSec = 0;
    double timeoutSec = 0;
    atomic<bool> timedOut{false};

    SigIndex(const Dataset& data, const TreeDecomp& td) : d(data), T(td), r((int)data.r) {}

    void initPruner(Pruner& p) const {
        p.r = r;
        p.mode = R.mode;
        p.lo = R.lo.data();
        p.hi = R.hi.data();
        p.center = R.center.data();
        p.worldFilter = useWorldFilter ? &d : nullptr;
    }

    // The elimination part of the index construction: shortcut signature sets.
    void buildShortcuts() {
        Timer tm;
        S.assign(2 * T.numShortcuts, {});
        for (int64_t e = 0; e < d.m; e++) {
            if (T.arcPlace[e] == UINT64_MAX) continue;
            auto& set = S[T.arcPlace[e]];
            set.insert(set.end(), d.U(e), d.U(e) + r);
        }
        Pruner pr;
        initPruner(pr);
        for (auto& set : S) pr.run(set);
        vector<float> pool;
        for (int idx = 0; idx < T.n; idx++) {
            if ((idx & 1023) == 0 && timeoutSec > 0 && tm.sec() > timeoutSec) {
                timedOut = true;
                return;
            }
            int v = T.order[idx];
            int k = (int)T.bag[v].size();
            uint64_t base = T.scBase[v], tri = T.triBase[v];
            for (int i = 0; i < k; i++) {
                const auto& vaDown = S[2 * (base + i) + 1];  // a -> v
                const auto& vaUp = S[2 * (base + i)];        // v -> a
                for (int j = i + 1; j < k; j++) {
                    uint64_t sc = T.triSc[tri++];
                    const auto& vbUp = S[2 * (base + j)];        // v -> b
                    const auto& vbDown = S[2 * (base + j) + 1];  // b -> v
                    if (!vaDown.empty() && !vbUp.empty()) {
                        auto& target = S[2 * sc];  // a -> b
                        pool = target;
                        joinInto(vaDown.data(), vaDown.size() / r, vbUp.data(), vbUp.size() / r, r, pool);
                        pr.run(pool);
                        target = pool;
                    }
                    if (!vbDown.empty() && !vaUp.empty()) {
                        auto& target = S[2 * sc + 1];  // b -> a
                        pool = target;
                        joinInto(vbDown.data(), vbDown.size() / r, vaUp.data(), vaUp.size() / r, r, pool);
                        pr.run(pool);
                        target = pool;
                    }
                }
            }
        }
        for (auto& set : S) set.shrink_to_fit();
        shortcutSec = tm.sec();
    }

    inline void entry(const Label& l, int depthW, const float*& p, size_t& cnt) const {
        uint32_t a = l.off[depthW], b = l.off[depthW + 1];
        p = l.sig.data() + a;
        cnt = (b - a) / r;
    }

    // The top-down part of the index construction, level by level so that a level can run in parallel.
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
#pragma omp for schedule(dynamic, 16)
#endif
                for (int ii = 0; ii < (int)vs.size(); ii++) {
                    if (timedOut) continue;
                    if (timeoutSec > 0 && (ii & 255) == 0 && tm.sec() > timeoutSec) {
                        timedOut = true;
                        continue;
                    }
                    int v = vs[ii], dv = T.depth[v];
                    anc.assign(dv, 0);
                    for (int x = T.parent[v]; x >= 0; x = T.parent[x]) anc[T.depth[x]] = x;
                    Label& up = L[0][v];
                    Label& dn = L[1][v];
                    up.off.assign(dv + 1, 0);
                    dn.off.assign(dv + 1, 0);
                    const auto& bag = T.bag[v];
                    uint64_t base = T.scBase[v];
                    for (int dw = 0; dw < dv; dw++) {
                        int w = anc[dw];
                        for (int dir = 0; dir < 2; dir++) {
                            pool.clear();
                            for (int i = 0; i < (int)bag.size(); i++) {
                                int u = bag[i];
                                const auto& sc = S[2 * (base + i) + dir];
                                if (sc.empty()) continue;
                                if (u == w) {
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
                                if (cnt) joinInto(sc.data(), sc.size() / r, p, cnt, r, pool);
                            }
                            pr.run(pool);
                            Label& lb = dir == 0 ? up : dn;
                            size_t offsetLimit = numeric_limits<uint32_t>::max();
                            if (lb.sig.size() > offsetLimit || pool.size() > offsetLimit - lb.sig.size())
                                throw length_error("index size limit: 32-bit label offsets");
                            lb.sig.insert(lb.sig.end(), pool.begin(), pool.end());
                            lb.off[dw + 1] = (uint32_t)lb.sig.size();
                        }
                    }
                    up.sig.shrink_to_fit();
                    dn.sig.shrink_to_fit();
                }
            }
            if (timedOut) return;
        }
        labelSec = tm.sec();
    }

    void build() {
        buildShortcuts();
        if (!timedOut) buildLabels();
    }

    uint64_t labelBytes() const {
        uint64_t b = 0;
        for (int dir = 0; dir < 2; dir++)
            for (const auto& l : L[dir]) b += l.off.size() * 4 + l.sig.size() * 4;
        return b;
    }
    uint64_t shortcutBytes() const {
        uint64_t b = 0;
        for (const auto& s : S) b += s.size() * 4;
        return b;
    }
    SizeStats labelSizeStats() const {
        vector<uint64_t> hist(1 << 16, 0);
        map<size_t, uint64_t> overflow;
        for (int dir = 0; dir < 2; dir++) {
            for (const auto& l : L[dir]) {
                for (size_t i = 0; i + 1 < l.off.size(); i++) {
                    size_t c = (l.off[i + 1] - l.off[i]) / r;
                    if (c < hist.size()) hist[c]++; else overflow[c]++;
                }
            }
        }
        return sizeStats(hist, overflow);
    }
    SizeStats shortcutSizeStats() const {
        vector<uint64_t> hist(1 << 16, 0);
        map<size_t, uint64_t> overflow;
        for (const auto& s : S) {
            size_t count = s.size() / r;
            if (count < hist.size()) hist[count]++; else overflow[count]++;
        }
        return sizeStats(hist, overflow);
    }

    // Lines 2-11 of the query algorithm. Fills hubs with the pairs that must be evaluated on Omega.
    // refSig receives c*, the signature of the fastest route at the mean world of Omega.
    void collect(int s, int t, const OmegaCtx& ctx, bool hubPrune, bool candPrune, vector<HubCand>& hubs,
                 vector<float>* refSig, QueryStats* st) const {
        hubs.clear();
        const vector<int>* sep;
        int c = T.separator(s, t, sep);
        const float* lo = ctx.lo.data();
        const float* hi = ctx.hi.data();
        const float* mean = ctx.mean.data();
        if (!sep) {  // ancestor-descendant case
            hubs.emplace_back();
            HubCand& h = hubs.back();
            const float* p;
            size_t cnt;
            if (c == t) {
                entry(L[0][s], T.depth[t], p, cnt);
                h.A.assign(p, p + cnt * r);
                h.zeroB = true;
            } else {
                entry(L[1][t], T.depth[s], p, cnt);
                h.B.assign(p, p + cnt * r);
                h.zeroA = true;
            }
            if (st) {
                st->hubs += 1;
                st->hubsKept += 1;
                st->cands += (long long)cnt;
            }
        }
        // Per-separator candidates and the reference pair (Line 8).
        struct Item {
            const float* a;
            size_t na;
            const float* b;
            size_t nb;
        };
        static thread_local vector<Item> items;
        items.clear();
        float bestVal = FINF;
        int bestItem = -1;
        const float* bestA = nullptr;
        const float* bestB = nullptr;
        if (sep) {
            for (int h : *sep) {
                int dh = T.depth[h];
                Item it;
                entry(L[0][s], dh, it.a, it.na);
                entry(L[1][t], dh, it.b, it.nb);
                if (!it.na || !it.nb) continue;
                float va = FINF, vb = FINF;
                const float* ba = nullptr;
                const float* bb = nullptr;
                for (size_t i = 0; i < it.na; i++) {
                    float x = dotf(it.a + i * r, mean, r);
                    if (x < va) {
                        va = x;
                        ba = it.a + i * r;
                    }
                }
                for (size_t i = 0; i < it.nb; i++) {
                    float x = dotf(it.b + i * r, mean, r);
                    if (x < vb) {
                        vb = x;
                        bb = it.b + i * r;
                    }
                }
                if (va + vb < bestVal) {
                    bestVal = va + vb;
                    bestItem = (int)items.size();
                    bestA = ba;
                    bestB = bb;
                }
                items.push_back(it);
            }
            if (st) {
                st->hubs += (long long)sep->size();
                for (const auto& it : items) st->cands += (long long)(it.na + it.nb);
            }
            static thread_local vector<float> ref;
            ref.resize(r);
            if (bestItem >= 0)
                for (int j = 0; j < r; j++) ref[j] = bestA[j] + bestB[j];
            if (refSig && bestItem >= 0) *refSig = ref;
            for (int idx = 0; idx < (int)items.size(); idx++) {
                const Item& it = items[idx];
                bool keep = !hubPrune || idx == bestItem;
                if (!keep) {  // Line 10: drop the separator if c* dominates all its pairs over B
                    for (size_t i = 0; i < it.na && !keep; i++) {
                        const float* a = it.a + i * r;
                        for (size_t q = 0; q < it.nb && !keep; q++) {
                            const float* b = it.b + q * r;
                            float sv = 0;
                            for (int j = 0; j < r; j++) {
                                float dd = ref[j] - a[j] - b[j];
                                sv += dd > 0 ? hi[j] * dd : lo[j] * dd;
                            }
                            if (sv > 0) keep = true;
                        }
                    }
                }
                if (!keep) continue;
                hubs.emplace_back();
                HubCand& hc = hubs.back();
                hc.A.assign(it.a, it.a + it.na * r);
                hc.B.assign(it.b, it.b + it.nb * r);
            }
        } else if (refSig) {
            // The single pair of an ancestor-descendant query.
            const HubCand& h = hubs[0];
            const vector<float>& side = h.zeroA ? h.B : h.A;
            size_t cnt = side.size() / (size_t)r;
            float best = FINF;
            refSig->assign(r, 0);
            for (size_t i = 0; i < cnt; i++) {
                float x = dotf(&side[i * r], mean, r);
                if (x < best) {
                    best = x;
                    copy(&side[i * r], &side[i * r] + r, refSig->begin());
                }
            }
        }
        if (candPrune) {  // Line 11: SigPrune with respect to B
            Pruner p;
            p.r = r;
            p.mode = REGION_AABB;
            p.lo = lo;
            p.hi = hi;
            p.center = mean;
            for (auto& hc : hubs) {
                if (!hc.A.empty()) p.run(hc.A);
                if (!hc.B.empty()) p.run(hc.B);
            }
        }
        if (st) {
            if (sep) st->hubsKept += (long long)hubs.size();
            for (const auto& hc : hubs) st->candsKept += (long long)((hc.A.size() + hc.B.size()) / r);
        }
    }

    // Lines 12-13 of the query algorithm: d(w) for every world of the context, written to vals.
    void evaluate(const vector<HubCand>& hubs, const float* zsoa, int k, float* vals) const {
        static thread_local vector<float> xa, xb, tmp;
        xa.resize(k);
        xb.resize(k);
        tmp.resize(k);
        fill(vals, vals + k, FINF);
        auto sideMin = [&](const vector<float>& set, float* out) {
            size_t cnt = set.size() / (size_t)r;
            fill(out, out + k, FINF);
            float* t = tmp.data();
            for (size_t c = 0; c < cnt; c++) {
                const float* sig = &set[c * r];
                fill(t, t + k, 0.f);
                for (int j = 0; j < r; j++) {
                    const float cj = sig[j];
                    const float* col = zsoa + (size_t)j * k;
                    for (int i = 0; i < k; i++) t[i] += cj * col[i];
                }
                for (int i = 0; i < k; i++) out[i] = t[i] < out[i] ? t[i] : out[i];
            }
        };
        for (const auto& hc : hubs) {
            if (hc.zeroA) fill(xa.begin(), xa.end(), 0.f);
            else sideMin(hc.A, xa.data());
            if (hc.zeroB) fill(xb.begin(), xb.end(), 0.f);
            else sideMin(hc.B, xb.data());
            for (int i = 0; i < k; i++) {
                float v = xa[i] + xb[i];
                vals[i] = v < vals[i] ? v : vals[i];
            }
        }
    }

    // A complete TTQ query. ranks must be sorted; out[i] receives the value of ranks[i].
    void answer(int s, int t, const int* idx, int cnt, const vector<int>& ranks, bool hubPrune, bool candPrune,
                vector<float>& out, vector<float>& vals, QueryStats* st, OmegaCtx& ctx, vector<HubCand>& hubs) const {
        out.resize(ranks.size());
        vals.resize(cnt);
        if (s == t) {
            fill(vals.begin(), vals.end(), 0.f);
            fill(out.begin(), out.end(), 0.f);
            return;
        }
        ctx.prepare(d, idx, cnt);
        collect(s, t, ctx, hubPrune, candPrune, hubs, nullptr, st);
        evaluate(hubs, ctx.zsoa.data(), cnt, vals.data());
        multiSelect(vals.data(), cnt, ranks, 0, (int)ranks.size());
        for (size_t i = 0; i < ranks.size(); i++) out[i] = vals[ranks[i]];
    }
};
