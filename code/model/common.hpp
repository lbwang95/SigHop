// Shared data structures: dataset I/O, time descriptions, CLI arguments, timers.
#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

static const float FINF = numeric_limits<float>::infinity();

struct Timer {
    chrono::steady_clock::time_point t0 = chrono::steady_clock::now();
    void reset() { t0 = chrono::steady_clock::now(); }
    double us() const { return chrono::duration<double, micro>(chrono::steady_clock::now() - t0).count(); }
    double sec() const { return us() / 1e6; }
};

// Directed road network plus low-rank traffic model:
//   travel time of arc e in world w = <u[e*r .. e*r+r), z[w*r .. w*r+r)>.
struct Dataset {
    int64_t n = 0, m = 0, r = 0, N = 0, days = 0, slotsPerDay = 0;
    double slotSeconds = 900;
    vector<int32_t> from, to;
    vector<float> u, z;
    vector<int32_t> wday, wslot, wdow;  // tags of each world; wdow: 0 = Monday
    vector<float> x, y;                 // optional vertex coordinates
    inline const float* U(int64_t e) const { return &u[e * r]; }
    inline const float* Z(int64_t w) const { return &z[w * r]; }
};

inline void writeAll(FILE* f, const void* p, size_t bytes) {
    if (bytes && fwrite(p, 1, bytes, f) != bytes) throw runtime_error("short write");
}
inline void readAll(FILE* f, void* p, size_t bytes) {
    if (bytes && fread(p, 1, bytes, f) != bytes) throw runtime_error("short read");
}

inline void saveDataset(const Dataset& d, const string& path) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) throw runtime_error("cannot write " + path);
    writeAll(f, "HTQDATA1", 8);
    int64_t hdr[7] = {d.n, d.m, d.r, d.N, d.days, d.slotsPerDay, d.x.empty() ? 0 : 1};
    writeAll(f, hdr, sizeof(hdr));
    writeAll(f, &d.slotSeconds, sizeof(double));
    writeAll(f, d.from.data(), d.m * 4);
    writeAll(f, d.to.data(), d.m * 4);
    writeAll(f, d.u.data(), d.m * d.r * 4);
    writeAll(f, d.z.data(), d.N * d.r * 4);
    writeAll(f, d.wday.data(), d.N * 4);
    writeAll(f, d.wslot.data(), d.N * 4);
    writeAll(f, d.wdow.data(), d.N * 4);
    if (!d.x.empty()) {
        writeAll(f, d.x.data(), d.n * 4);
        writeAll(f, d.y.data(), d.n * 4);
    }
    fclose(f);
}

inline Dataset loadDataset(const string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) throw runtime_error("cannot read " + path);
    char magic[8];
    readAll(f, magic, 8);
    if (memcmp(magic, "HTQDATA1", 8) != 0) throw runtime_error("bad magic in " + path);
    int64_t hdr[7];
    readAll(f, hdr, sizeof(hdr));
    Dataset d;
    d.n = hdr[0]; d.m = hdr[1]; d.r = hdr[2]; d.N = hdr[3]; d.days = hdr[4]; d.slotsPerDay = hdr[5];
    readAll(f, &d.slotSeconds, sizeof(double));
    d.from.resize(d.m); d.to.resize(d.m); d.u.resize(d.m * d.r); d.z.resize(d.N * d.r);
    d.wday.resize(d.N); d.wslot.resize(d.N); d.wdow.resize(d.N);
    readAll(f, d.from.data(), d.m * 4);
    readAll(f, d.to.data(), d.m * 4);
    readAll(f, d.u.data(), d.m * d.r * 4);
    readAll(f, d.z.data(), d.N * d.r * 4);
    readAll(f, d.wday.data(), d.N * 4);
    readAll(f, d.wslot.data(), d.N * 4);
    readAll(f, d.wdow.data(), d.N * 4);
    if (hdr[6]) {
        d.x.resize(d.n); d.y.resize(d.n);
        readAll(f, d.x.data(), d.n * 4);
        readAll(f, d.y.data(), d.n * 4);
    }
    fclose(f);
    return d;
}

// Observed travel times per world, used as the reference "raw worlds".
struct RawWorlds {
    int64_t N = 0, nnz = 0;
    vector<int64_t> off;
    vector<int32_t> edge;
    vector<float> w;
};
inline void saveRaw(const RawWorlds& rw, const string& path) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) throw runtime_error("cannot write " + path);
    writeAll(f, "HTQRAW01", 8);
    int64_t hdr[2] = {rw.N, rw.nnz};
    writeAll(f, hdr, sizeof(hdr));
    writeAll(f, rw.off.data(), (rw.N + 1) * 8);
    writeAll(f, rw.edge.data(), rw.nnz * 4);
    writeAll(f, rw.w.data(), rw.nnz * 4);
    fclose(f);
}
inline RawWorlds loadRaw(const string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) throw runtime_error("cannot read " + path);
    char magic[8];
    readAll(f, magic, 8);
    if (memcmp(magic, "HTQRAW01", 8) != 0) throw runtime_error("bad magic in " + path);
    int64_t hdr[2];
    readAll(f, hdr, sizeof(hdr));
    RawWorlds rw;
    rw.N = hdr[0]; rw.nnz = hdr[1];
    rw.off.resize(rw.N + 1); rw.edge.resize(rw.nnz); rw.w.resize(rw.nnz);
    readAll(f, rw.off.data(), (rw.N + 1) * 8);
    readAll(f, rw.edge.data(), rw.nnz * 4);
    readAll(f, rw.w.data(), rw.nnz * 4);
    fclose(f);
    return rw;
}

// Time description: conjunction of a day-of-week set, a slot range and a day range.
struct Phi {
    string name;
    int dowMask = 127;
    int slotLo = 0, slotHi = numeric_limits<int>::max();  // [slotLo, slotHi)
    int dayLo = 0, dayHi = numeric_limits<int>::max();    // [dayLo, dayHi)
    bool match(int day, int slot, int dow) const {
        return ((dowMask >> dow) & 1) && slot >= slotLo && slot < slotHi && day >= dayLo && day < dayHi;
    }
};

inline vector<string> splitStr(const string& s, char sep) {
    vector<string> out;
    string cur;
    stringstream ss(s);
    while (getline(ss, cur, sep))
        if (!cur.empty()) out.push_back(cur);
    return out;
}

// Grammar: "name:dow=0-4;slots=32-36;days=0-100". Absent fields match everything.
inline Phi parsePhi(const string& spec) {
    Phi p;
    string body = spec;
    size_t c = spec.find(':');
    if (c != string::npos) {
        p.name = spec.substr(0, c);
        body = spec.substr(c + 1);
    } else {
        p.name = spec;
    }
    for (const string& kv : splitStr(body, ';')) {
        size_t e = kv.find('=');
        if (e == string::npos) continue;
        string k = kv.substr(0, e), v = kv.substr(e + 1);
        auto range = [&](int& lo, int& hi) {
            size_t d = v.find('-');
            lo = stoi(v.substr(0, d));
            hi = stoi(v.substr(d + 1));
        };
        if (k == "dow") {
            p.dowMask = 0;
            for (const string& part : splitStr(v, ',')) {
                size_t d = part.find('-');
                int a = stoi(part.substr(0, d)), b = d == string::npos ? a : stoi(part.substr(d + 1));
                for (int i = a; i <= b; i++) p.dowMask |= 1 << i;
            }
        } else if (k == "slots") {
            range(p.slotLo, p.slotHi);
        } else if (k == "days") {
            range(p.dayLo, p.dayHi);
        }
    }
    return p;
}

inline vector<int> worldsOf(const Dataset& d, const Phi& p) {
    vector<int> idx;
    for (int w = 0; w < d.N; w++)
        if (p.match(d.wday[w], d.wslot[w], d.wdow[w])) idx.push_back(w);
    return idx;
}

struct Args {
    map<string, string> kv;
    Args(int argc, char** argv, int start) {
        for (int i = start; i < argc; i++) {
            string a = argv[i];
            if (a.rfind("--", 0) != 0) throw runtime_error("unexpected argument " + a);
            if (i + 1 < argc && string(argv[i + 1]).rfind("--", 0) != 0) {
                kv[a.substr(2)] = argv[++i];
            } else {
                kv[a.substr(2)] = "1";
            }
        }
    }
    string get(const string& k, const string& def = "") const {
        auto it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
    string need(const string& k) const {
        auto it = kv.find(k);
        if (it == kv.end()) throw runtime_error("missing --" + k);
        return it->second;
    }
    int geti(const string& k, int def) const { return kv.count(k) ? stoi(kv.at(k)) : def; }
    double getd(const string& k, double def) const { return kv.count(k) ? stod(kv.at(k)) : def; }
};

inline vector<pair<int, int>> loadQueries(const string& path) {
    vector<pair<int, int>> q;
    ifstream in(path);
    if (!in) throw runtime_error("cannot read " + path);
    int s, t;
    while (in >> s >> t) q.push_back({s, t});
    return q;
}

// CSV writer with a fixed header; appends when the file already has content.
struct Csv {
    FILE* f = nullptr;
    Csv(const string& path, const string& header) {
        bool exists = false;
        {
            ifstream in(path);
            exists = in.good() && in.peek() != ifstream::traits_type::eof();
        }
        f = fopen(path.c_str(), "a");
        if (!f) throw runtime_error("cannot write " + path);
        if (!exists) fprintf(f, "%s\n", header.c_str());
        fflush(f);
    }
    ~Csv() {
        if (f) fclose(f);
    }
    void row(const string& s) {
        fprintf(f, "%s\n", s.c_str());
        fflush(f);
    }
};

template <class... T>
inline string joinCsv(const T&... v) {
    stringstream ss;
    ss.precision(10);
    int i = 0;
    ((ss << (i++ ? "," : "") << v), ...);
    return ss.str();
}

// Multi-selection: after the call, a[r] holds the r-th smallest for every rank r in ranks[lo, hi).
inline void multiSelect(float* a, int n, const vector<int>& ranks, int lo, int hi, int base = 0) {
    // Equal ranks can occur when several levels fall on one value; those are already placed.
    while (lo < hi && ranks[lo] - base < 0) lo++;
    while (lo < hi && ranks[hi - 1] - base >= n) hi--;
    if (lo >= hi || n <= 0) return;
    int mid = (lo + hi) / 2;
    int k = ranks[mid] - base;
    nth_element(a, a + k, a + n);
    multiSelect(a, k, ranks, lo, mid, base);
    multiSelect(a + k + 1, n - k - 1, ranks, mid + 1, hi, base + k + 1);
}

// 0-based rank of Q_alpha among k values: ceil(alpha * k) - 1, with a guard for rounding.
inline int quantileRank(double alpha, int k) {
    int c = (int)ceil(alpha * k - 1e-9);
    return max(0, min(k - 1, c - 1));
}

// Inverse standard normal CDF (Acklam's rational approximation).
inline double invNorm(double p) {
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                               3.754408661907416e+00};
    const double pl = 0.02425;
    if (p <= 0) return -1e30;
    if (p >= 1) return 1e30;
    if (p < pl) {
        double q = sqrt(-2 * log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    if (p <= 1 - pl) {
        double q = p - 0.5, t = q * q;
        return (((((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4]) * t + a[5]) * q /
               (((((b[0] * t + b[1]) * t + b[2]) * t + b[3]) * t + b[4]) * t + 1);
    }
    double q = sqrt(-2 * log(1 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
}
