// Fits the low-rank traffic model of the paper on observed segment speeds.
// Input (from scripts/prepare_mets10.py) in --in DIR:
//   meta.txt  : n m days slots_per_day slot_seconds slot_lo slot_hi nnz classes
//   arcs.bin  : int32 from[m], int32 to[m], float len_m[m], float ff_kph[m], int32 cls[m]
//   coords.bin: float x[n], float y[n]
//   days.bin  : int32 dow[days]
//   obs.bin   : records of (int32 edge, int32 world, float speed_kph)
// Output: --out DIR/dataset_r<r>.bin for each r in --ranks, DIR/raw.bin, DIR/nmf_report.csv
#include "common.hpp"

struct Obs {
    int32_t e, w;
    float y;  // delay in s/m
    float pace;  // observed pace, before clipping the delay at zero
};

int main(int argc, char** argv) {
    try {
        Args a(argc, argv, 1);
        string in = a.need("in"), outDir = a.get("out", in);
        vector<int> ranks;
        for (const string& s : splitStr(a.get("ranks", "2,4,8,16,32"), ',')) ranks.push_back(stoi(s));
        int iters = a.geti("iters", 60);
        double minSpeed = a.getd("min-speed", 3.0), holdout = a.getd("holdout", 0.05);
        int minObs = a.geti("min-obs", 20);
        string name = a.get("name", in);

        int64_t n, m, days, slots, slotLo, slotHi, nnz, classes;
        double slotSeconds;
        {
            ifstream f(in + "/meta.txt");
            if (!f) throw runtime_error("cannot read meta.txt");
            f >> n >> m >> days >> slots >> slotSeconds >> slotLo >> slotHi >> nnz >> classes;
        }
        vector<int32_t> from(m), to(m), cls(m);
        vector<float> lenM(m), ff(m);
        {
            FILE* f = fopen((in + "/arcs.bin").c_str(), "rb");
            if (!f) throw runtime_error("cannot read arcs.bin");
            readAll(f, from.data(), m * 4);
            readAll(f, to.data(), m * 4);
            readAll(f, lenM.data(), m * 4);
            readAll(f, ff.data(), m * 4);
            readAll(f, cls.data(), m * 4);
            fclose(f);
        }
        vector<float> x(n), y(n);
        {
            FILE* f = fopen((in + "/coords.bin").c_str(), "rb");
            if (!f) throw runtime_error("cannot read coords.bin");
            readAll(f, x.data(), n * 4);
            readAll(f, y.data(), n * 4);
            fclose(f);
        }
        vector<int32_t> dow(days);
        {
            FILE* f = fopen((in + "/days.bin").c_str(), "rb");
            if (!f) throw runtime_error("cannot read days.bin");
            readAll(f, dow.data(), days * 4);
            fclose(f);
        }
        int64_t slotCnt = slotHi - slotLo, N = days * slotCnt;
        vector<int32_t> slotTags(slotCnt);
        iota(slotTags.begin(), slotTags.end(), (int32_t)slotLo);
        if (FILE* tags = fopen((in + "/slots.bin").c_str(), "rb")) {
            readAll(tags, slotTags.data(), slotCnt * sizeof(int32_t));
            if (fgetc(tags) != EOF) throw runtime_error("slot tags have extra records");
            fclose(tags);
            for (int64_t i = 0; i < slotCnt; i++)
                if (slotTags[i] < 0 || slotTags[i] >= slots || (i && slotTags[i] <= slotTags[i-1]))
                    throw runtime_error("invalid actual slot tags");
        }
        vector<float> pace0(m);  // free-flow pace, s/m
        for (int64_t e = 0; e < m; e++) pace0[e] = 3.6f / max(ff[e], (float)minSpeed);

        struct Rec {
            int32_t e, w;
            float kph;
        };
        vector<Obs> fit, test;
        RawWorlds raw;
        raw.N = N;
        raw.off.assign(N + 1, 0);
        {
            FILE* f = fopen((in + "/obs.bin").c_str(), "rb");
            if (!f) throw runtime_error("cannot read obs.bin");
            vector<Rec> all(nnz);
            readAll(f, all.data(), nnz * sizeof(Rec));
            fclose(f);
            vector<uint8_t> holdoutMask;
            if (a.kv.count("holdout-mask")) {
                holdoutMask.resize(nnz);
                FILE* mask = fopen(a.need("holdout-mask").c_str(), "rb");
                if (!mask) throw runtime_error("cannot read holdout mask");
                readAll(mask, holdoutMask.data(), nnz);
                if (fgetc(mask) != EOF) throw runtime_error("holdout mask has extra records");
                fclose(mask);
                for (uint8_t flag : holdoutMask)
                    if (flag > 1) throw runtime_error("holdout mask must contain only 0 or 1");
            }
            mt19937_64 rng(a.geti("seed", 1));
            uniform_real_distribution<double> U(0, 1);
            fit.reserve(nnz);
            for (size_t qi = 0; qi < all.size(); qi++) {
                const Rec& q = all[qi];
                if (q.e < 0 || q.e >= m || q.w < 0 || q.w >= N || !(q.kph > 0)) continue;
                float pace = 3.6f / max(q.kph, (float)minSpeed);
                raw.off[q.w + 1]++;
                Obs o{q.e, q.w, max(0.f, pace - pace0[q.e]), pace};
                bool isTest = holdoutMask.empty() ? U(rng) < holdout : holdoutMask[qi] != 0;
                if (isTest) test.push_back(o);
                else fit.push_back(o);
            }
            for (int64_t w = 0; w < N; w++) raw.off[w + 1] += raw.off[w];
            raw.nnz = raw.off[N];
            raw.edge.resize(raw.nnz);
            raw.w.resize(raw.nnz);
            vector<int64_t> pos(raw.off.begin(), raw.off.end() - 1);
            for (const Rec& q : all) {
                if (q.e < 0 || q.e >= m || q.w < 0 || q.w >= N || !(q.kph > 0)) continue;
                float pace = 3.6f / max(q.kph, (float)minSpeed);
                int64_t p = pos[q.w]++;
                raw.edge[p] = q.e;
                raw.w[p] = pace * lenM[q.e];
            }
        }
        saveRaw(raw, outDir + "/raw.bin");

        // Group the fitting observations by edge and by world.
        vector<int64_t> he(m + 1, 0), hw(N + 1, 0);
        for (const Obs& o : fit) {
            he[o.e + 1]++;
            hw[o.w + 1]++;
        }
        for (int64_t i = 0; i < m; i++) he[i + 1] += he[i];
        for (int64_t i = 0; i < N; i++) hw[i + 1] += hw[i];
        vector<int32_t> ew(fit.size()), we(fit.size());
        vector<float> ey(fit.size()), wy(fit.size());
        {
            vector<int64_t> pe(he.begin(), he.end() - 1), pw(hw.begin(), hw.end() - 1);
            for (const Obs& o : fit) {
                ew[pe[o.e]] = o.w;
                ey[pe[o.e]++] = o.y;
                we[pw[o.w]] = o.e;
                wy[pw[o.w]++] = o.y;
            }
        }

        Csv report(outDir + "/nmf_report.csv", "dataset,r,fit_obs,test_obs,fit_rmse_spm,test_mape_pace,base_mape_pace,seconds");
        double baseMape = 0;
        for (const Obs& o : test) {
            double p0 = pace0[o.e], p = o.pace;
            baseMape += fabs(p - p0) / p;
        }
        baseMape /= max<size_t>(1, test.size());

        for (int r : ranks) {
            Timer tm;
            int rd = r - 1;  // delay factors; factor 0 is the free-flow time with z = 1
            vector<float> U((size_t)m * max(rd, 1), 0), Z((size_t)N * max(rd, 1), 0);
            double rmse = 0, mape = baseMape;
            if (rd > 0) {
                mt19937 rng(a.geti("seed", 1) + r);
                uniform_real_distribution<float> U01(0.01f, 1.f);
                double meanY = 0;
                for (const Obs& o : fit) meanY += o.y;
                meanY /= max<size_t>(1, fit.size());
                float initScale = (float)sqrt(max(1e-9, meanY) / rd);
                for (auto& v : U) v = U01(rng) * initScale;
                for (auto& v : Z) v = U01(rng) * initScale;
                // Masked multiplicative updates (Lee and Seung).
                for (int it = 0; it < iters; it++) {
#ifdef _OPENMP
#pragma omp parallel
#endif
                    {
                        vector<double> num(rd), den(rd);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
                        for (int64_t e = 0; e < m; e++) {
                            if (he[e] == he[e + 1]) continue;
                            fill(num.begin(), num.end(), 0.0);
                            fill(den.begin(), den.end(), 0.0);
                            float* ue = &U[(size_t)e * rd];
                            for (int64_t i = he[e]; i < he[e + 1]; i++) {
                                const float* zw = &Z[(size_t)ew[i] * rd];
                                double pred = 0;
                                for (int j = 0; j < rd; j++) pred += ue[j] * zw[j];
                                for (int j = 0; j < rd; j++) {
                                    num[j] += ey[i] * zw[j];
                                    den[j] += pred * zw[j];
                                }
                            }
                            for (int j = 0; j < rd; j++) ue[j] = (float)(ue[j] * num[j] / (den[j] + 1e-12));
                        }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
                        for (int64_t w = 0; w < N; w++) {
                            if (hw[w] == hw[w + 1]) continue;
                            fill(num.begin(), num.end(), 0.0);
                            fill(den.begin(), den.end(), 0.0);
                            float* zw = &Z[(size_t)w * rd];
                            for (int64_t i = hw[w]; i < hw[w + 1]; i++) {
                                const float* ue = &U[(size_t)we[i] * rd];
                                double pred = 0;
                                for (int j = 0; j < rd; j++) pred += ue[j] * zw[j];
                                for (int j = 0; j < rd; j++) {
                                    num[j] += wy[i] * ue[j];
                                    den[j] += pred * ue[j];
                                }
                            }
                            for (int j = 0; j < rd; j++) zw[j] = (float)(zw[j] * num[j] / (den[j] + 1e-12));
                        }
                    }
                }
                // Worlds without observations take the mean factor vector of their slot.
                {
                    vector<double> slotSum((size_t)slotCnt * rd, 0);
                    vector<int> slotCount(slotCnt, 0);
                    for (int64_t w = 0; w < N; w++) {
                        if (hw[w] == hw[w + 1]) continue;
                        int s = (int)(w % slotCnt);
                        slotCount[s]++;
                        for (int j = 0; j < rd; j++) slotSum[(size_t)s * rd + j] += Z[(size_t)w * rd + j];
                    }
                    for (int64_t w = 0; w < N; w++) {
                        if (hw[w] != hw[w + 1]) continue;
                        int s = (int)(w % slotCnt);
                        for (int j = 0; j < rd; j++)
                            Z[(size_t)w * rd + j] = slotCount[s] ? (float)(slotSum[(size_t)s * rd + j] / slotCount[s]) : 0;
                    }
                }
                // Edges with too few observations take the mean loading of their road class.
                {
                    vector<double> sum((size_t)classes * rd, 0);
                    vector<int64_t> cnt(classes, 0);
                    for (int64_t e = 0; e < m; e++) {
                        if (he[e + 1] - he[e] < minObs) continue;
                        cnt[cls[e]]++;
                        for (int j = 0; j < rd; j++) sum[(size_t)cls[e] * rd + j] += U[(size_t)e * rd + j];
                    }
                    for (int64_t e = 0; e < m; e++) {
                        if (he[e + 1] - he[e] >= minObs) continue;
                        for (int j = 0; j < rd; j++)
                            U[(size_t)e * rd + j] = cnt[cls[e]] ? (float)(sum[(size_t)cls[e] * rd + j] / cnt[cls[e]]) : 0;
                    }
                }
                double se = 0;
                for (const Obs& o : fit) {
                    double p = 0;
                    for (int j = 0; j < rd; j++) p += U[(size_t)o.e * rd + j] * Z[(size_t)o.w * rd + j];
                    se += (p - o.y) * (p - o.y);
                }
                rmse = sqrt(se / max<size_t>(1, fit.size()));
                mape = 0;
                for (const Obs& o : test) {
                    double p = 0;
                    for (int j = 0; j < rd; j++) p += U[(size_t)o.e * rd + j] * Z[(size_t)o.w * rd + j];
                    double truth = o.pace;
                    mape += fabs(pace0[o.e] + p - truth) / truth;
                }
                mape /= max<size_t>(1, test.size());
            }
            Dataset d;
            d.n = n; d.m = m; d.r = r; d.N = N; d.days = days; d.slotsPerDay = slots;
            d.slotSeconds = slotSeconds;
            d.from = from; d.to = to; d.x = x; d.y = y;
            d.u.resize((size_t)m * r);
            d.z.resize((size_t)N * r);
            d.wday.resize(N); d.wslot.resize(N); d.wdow.resize(N);
            for (int64_t e = 0; e < m; e++) {
                float len = max(lenM[e], 0.1f);
                d.u[(size_t)e * r] = len * pace0[e];
                for (int j = 0; j < rd; j++) d.u[(size_t)e * r + j + 1] = len * U[(size_t)e * rd + j];
            }
            for (int64_t w = 0; w < N; w++) {
                d.z[(size_t)w * r] = 1;
                for (int j = 0; j < rd; j++) d.z[(size_t)w * r + j + 1] = Z[(size_t)w * rd + j];
                d.wday[w] = (int)(w / slotCnt);
                d.wslot[w] = slotTags[w % slotCnt];
                d.wdow[w] = dow[d.wday[w]];
            }
            saveDataset(d, outDir + "/dataset_r" + to_string(r) + ".bin");
            report.row(joinCsv(name, r, fit.size(), test.size(), rmse, mape, baseMape, tm.sec()));
            printf("r=%d fit_rmse=%.5g test_mape=%.4f base_mape=%.4f (%.1fs)\n", r, rmse, mape, baseMape, tm.sec());
            fflush(stdout);
        }
    } catch (const exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
