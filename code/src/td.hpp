// Tree decomposition by minimum degree elimination (MDE), shared by all indexes.
#pragma once
#include "common.hpp"

struct TreeDecomp {
    int n = 0;
    vector<int> order, rank_, parent, depth, jump;
    vector<vector<int>> bag;        // X(v) \ {v}, sorted by elimination rank (parent first)
    vector<uint64_t> scBase;        // shortcut id of (v, bag[v][i]) is scBase[v] + i
    uint64_t numShortcuts = 0;
    vector<uint64_t> triBase;       // per vertex, offset into triSc
    vector<uint64_t> triSc;         // for v and i<j: shortcut id of (bag[v][i], bag[v][j])
    vector<uint64_t> arcPlace;      // per arc: 2*sc + dir, dir 0 = lower->upper, 1 = upper->lower
    int width = 0, height = 0;      // treewidth (max |X(v)|), treeheight (max depth)
    double buildSec = 0;

    void build(const Dataset& d) {
        Timer tm;
        n = (int)d.n;
        vector<vector<int>> adj(n);
        for (int64_t e = 0; e < d.m; e++) {
            int a = d.from[e], b = d.to[e];
            if (a == b) continue;
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
        for (auto& v : adj) {
            sort(v.begin(), v.end());
            v.erase(unique(v.begin(), v.end()), v.end());
        }
        using P = pair<int, int>;
        priority_queue<P, vector<P>, greater<P>> pq;
        for (int v = 0; v < n; v++) pq.push({(int)adj[v].size(), v});
        order.clear();
        order.reserve(n);
        rank_.assign(n, -1);
        bag.assign(n, {});
        vector<int> merged;
        while (!pq.empty()) {
            auto [dg, v] = pq.top();
            pq.pop();
            if (rank_[v] >= 0 || dg != (int)adj[v].size()) continue;
            rank_[v] = (int)order.size();
            order.push_back(v);
            bag[v] = adj[v];
            for (int a : bag[v]) {
                merged.clear();
                const auto& x = adj[a];
                const auto& y = bag[v];
                size_t i = 0, j = 0;
                while (i < x.size() || j < y.size()) {
                    int w;
                    if (j == y.size() || (i < x.size() && x[i] <= y[j])) {
                        w = x[i];
                        if (j < y.size() && y[j] == w) j++;
                        i++;
                    } else {
                        w = y[j++];
                    }
                    if (w != a && w != v) merged.push_back(w);
                }
                adj[a] = merged;
                pq.push({(int)adj[a].size(), a});
            }
            vector<int>().swap(adj[v]);
        }
        parent.assign(n, -1);
        depth.assign(n, 0);
        jump.assign(n, 0);
        scBase.assign(n, 0);
        numShortcuts = 0;
        width = 0;
        height = 0;
        for (int v = 0; v < n; v++) {
            sort(bag[v].begin(), bag[v].end(), [&](int a, int b) { return rank_[a] < rank_[b]; });
            if (!bag[v].empty()) parent[v] = bag[v][0];
            scBase[v] = numShortcuts;
            numShortcuts += bag[v].size();
            width = max(width, (int)bag[v].size() + 1);
        }
        int roots = 0;
        for (int i = n - 1; i >= 0; i--) {
            int v = order[i], p = parent[v];
            if (p < 0) {
                roots++;
                depth[v] = 0;
                jump[v] = v;
            } else {
                depth[v] = depth[p] + 1;
                int jp = jump[p];
                jump[v] = (depth[p] - depth[jp] == depth[jp] - depth[jump[jp]]) ? jump[jp] : p;
            }
            height = max(height, depth[v]);
        }
        if (roots != 1) throw runtime_error("graph must be connected, found " + to_string(roots) + " roots");
        triBase.assign(n + 1, 0);
        uint64_t tot = 0;
        for (int v = 0; v < n; v++) {
            triBase[v] = tot;
            uint64_t k = bag[v].size();
            tot += k * (k - 1) / 2;
        }
        triBase[n] = tot;
        triSc.resize(tot);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
        for (int v = 0; v < n; v++) {
            uint64_t t = triBase[v];
            int k = (int)bag[v].size();
            for (int i = 0; i < k; i++) {
                for (int j = i + 1; j < k; j++) {
                    triSc[t++] = scOf(bag[v][i], bag[v][j]);
                }
            }
        }
        arcPlace.resize(d.m);
        for (int64_t e = 0; e < d.m; e++) {
            int a = d.from[e], b = d.to[e];
            if (a == b) {
                arcPlace[e] = UINT64_MAX;
                continue;
            }
            if (rank_[a] < rank_[b]) arcPlace[e] = 2 * scOf(a, b);
            else arcPlace[e] = 2 * scOf(b, a) + 1;
        }
        buildSec = tm.sec();
    }

    // Shortcut id of (lower, upper), where upper is in bag[lower].
    uint64_t scOf(int lower, int upper) const {
        const auto& b = bag[lower];
        int lo = 0, hi = (int)b.size() - 1, ru = rank_[upper];
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (rank_[b[mid]] < ru) lo = mid + 1;
            else hi = mid;
        }
        if (b.empty() || b[lo] != upper) throw runtime_error("shortcut lookup failed");
        return scBase[lower] + lo;
    }

    int levelAncestor(int v, int d) const {
        while (depth[v] > d) {
            if (depth[jump[v]] >= d) v = jump[v];
            else v = parent[v];
        }
        return v;
    }

    int lca(int a, int b) const {
        if (depth[a] < depth[b]) swap(a, b);
        a = levelAncestor(a, depth[b]);
        while (a != b) {
            if (jump[a] != jump[b]) {
                a = jump[a];
                b = jump[b];
            } else {
                a = parent[a];
                b = parent[b];
            }
        }
        return a;
    }

    // Both candidate separators X(c_s)\{c_s} and X(c_t)\{c_t} of Lemma 1. Returns the LCA.
    // When s and t are in an ancestor relation, a and b stay null.
    int separatorPair(int s, int t, const vector<int>*& a, const vector<int>*& b) const {
        int c = lca(s, t);
        a = b = nullptr;
        if (c == s || c == t) return c;
        a = &bag[levelAncestor(s, depth[c] + 1)];
        b = &bag[levelAncestor(t, depth[c] + 1)];
        return c;
    }

    // Separator of Lemma 3: the smaller of X(c_s)\{c_s}, X(c_t)\{c_t}. Returns the LCA.
    // When s and t are in an ancestor relation, sep stays null.
    int separator(int s, int t, const vector<int>*& sep) const {
        int c = lca(s, t);
        sep = nullptr;
        if (c == s || c == t) return c;
        int cs = levelAncestor(s, depth[c] + 1), ct = levelAncestor(t, depth[c] + 1);
        sep = bag[cs].size() <= bag[ct].size() ? &bag[cs] : &bag[ct];
        return c;
    }
};
