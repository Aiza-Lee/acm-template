#pragma once
#include "aizalib.h"

/*
 * CycleCounting·三四元环计数
 *
 * Overview:
 *     基于按顶点度数定向的无向简单图三元环与四元环统计模板。
 *     定义全序 u < v 当且仅当 deg[u] < deg[v] 或 (deg[u] == deg[v] 且 u < v):
 *     - 三元环: 将无向边从序小点指向序大点构建 DAG，每个点出度不超过 sqrt(2m)，
 *       从而枚举所有定向路径 u -> v -> w 并检查边 (u, w) 在 O(m sqrt(m))
 *       内不重不漏访问每个三元环一次。
 *     - 四元环: 固定最大点 u，经所有序小邻点 v (< u) 再走向序小邻点 w (< u) 统计
 *       2-路径数 cnt[w]，总计耗费 O(m sqrt(m))；亦可在同复杂度下统计每个点、
 *       每条边所在的四元环数量或直接定位一个四元环。
 *
 * API:
 *     struct Cycle4Stats:
 *         total      — 图中四元环总数
 *         vertex_cnt — vertex_cnt[u] 为顶点 u 所在的四元环个数 (1-based)
 *         edge_cnt   — edge_cnt[eid] 为边 eid 所在的四元环个数 (1-based)
 *     CycleCounter(n)         — 构造 n 个顶点的图 (1-based)
 *     CycleCounter(n, edges)  — 从边表批量构造图
 *     add_edge(u, v)          — 添加无向边 (u, v)，返回 1-based 边编号
 *     for_each_3_cycle(cb)    — 枚举每个三元环一次，支持 cb(u, v, w) 或 cb(u, v, w,
 *                                e_uv, e_vw, e_uw)，复杂度 O(m sqrt(m))
 *     count_3_cycles()        — 统计无向简单图三元环总数，复杂度 O(m sqrt(m))
 *     count_4_cycles()        — 统计无向简单图四元环总数，复杂度 O(m sqrt(m))
 *     count_4_cycles_detail() — 统计四元环总数及各点/各边贡献，复杂度 O(m sqrt(m))
 *     find_4_cycle()          — 查找任意一个四元环 {u, v1, w, v2}，无则空，复杂度
 *                                O(min(n^2, m sqrt(m)))
 *
 * Notes:
 *     1. Time: 预处理 O(n + m)，三元环与四元环查询均摊/最坏 O(m sqrt(m))。
 *     2. Space: O(n + m)。
 *     3. 1-based indexing。构造时自动忽略自环并按 O(n + m)
 *        桶去重重边（仅保留首条边编号）。
 *     4. 若用于有向图三元环计数，可在外部记录每条边编号对应的原方向，通过
 *        for_each_3_cycle(u, v, w, e_uv, e_vw, e_uw) 核对三条边是否同向成环。
 *
 * Related:
 *     - BronKerbosch·最大团.cpp::MaxClique — 一般图极大团枚举
 */
struct CycleCounter {
    struct Edge {
        int to, id;
    };
    struct Cycle4Stats {
        i64 total = 0;
        std::vector<i64> vertex_cnt;
        std::vector<i64> edge_cnt;
    };

    int n, raw_m;
    bool built;
    std::vector<std::pair<int, int>> raw_edges;
    std::vector<int> deg;
    std::vector<std::vector<Edge>> adj, dag;

    CycleCounter(int n = 0) : n(n), raw_m(0), built(false) {}

    CycleCounter(int n, const std::vector<std::pair<int, int>>& edges)
        : n(n), raw_m(0), built(false) {
        raw_edges.reserve(edges.size());
        for (auto [u, v] : edges) add_edge(u, v);
    }

    int add_edge(int u, int v) {
        AST(1 <= u && u <= n && 1 <= v && v <= n);
        raw_edges.emplace_back(u, v);
        built = false;
        return ++raw_m;
    }

    bool _cmp(int u, int v) const {
        return deg[u] != deg[v] ? deg[u] < deg[v] : u < v;
    }

    void _build() {
        if (built) return;
        built = true;
        deg.assign(n + 1, 0);
        adj.assign(n + 1, {});
        dag.assign(n + 1, {});

        std::vector<std::vector<Edge>> bucket(n + 1);
        rep(id, 1, raw_m) {
            auto [u, v] = raw_edges[id - 1];
            if (u == v) continue;
            if (u > v) std::swap(u, v);
            bucket[u].push_back({v, id});
        }

        std::vector<int> seen(n + 1, 0);
        std::vector<std::tuple<int, int, int>> simple_edges;
        simple_edges.reserve(raw_m);
        rep(u, 1, n) {
            for (auto [v, id] : bucket[u]) {
                if (seen[v] == u) continue;
                seen[v] = u;
                simple_edges.emplace_back(u, v, id);
                ++deg[u];
                ++deg[v];
            }
        }

        for (auto [u, v, id] : simple_edges) {
            adj[u].push_back({v, id});
            adj[v].push_back({u, id});
            if (_cmp(u, v)) dag[u].push_back({v, id});
            else dag[v].push_back({u, id});
        }
    }

    template<typename F>
    bool _for_2_paths(int u, F&& fn) const {
        for (auto [v, e1] : adj[u]) {
            if (!_cmp(v, u)) continue;
            for (auto [w, e2] : adj[v]) {
                if (!_cmp(w, u)) continue;
                using R = std::invoke_result_t<F, int, int, int, int>;
                if constexpr (std::is_same_v<R, bool>) {
                    if (fn(v, w, e1, e2)) return true;
                } else {
                    fn(v, w, e1, e2);
                }
            }
        }
        return false;
    }

    Cycle4Stats _solve_4_cycles(bool need_detail) {
        _build();
        Cycle4Stats res;
        if (need_detail) {
            res.vertex_cnt.assign(n + 1, 0);
            res.edge_cnt.assign(raw_m + 1, 0);
        }
        std::vector<int> cnt(n + 1, 0), touched;
        rep(u, 1, n) {
            _for_2_paths(u, [&](int, int w, int, int) {
                if (cnt[w] == 0) touched.push_back(w);
                res.total += cnt[w];
                ++cnt[w];
            });
            if (need_detail && !touched.empty()) {
                i64 u_add = 0;
                for (int w : touched) {
                    i64 pairs = 1LL * cnt[w] * (cnt[w] - 1) / 2;
                    res.vertex_cnt[w] += pairs;
                    u_add += pairs;
                }
                res.vertex_cnt[u] += u_add;
                _for_2_paths(u, [&](int v, int w, int e1, int e2) {
                    int ways = cnt[w] - 1;
                    res.vertex_cnt[v] += ways;
                    res.edge_cnt[e1] += ways;
                    res.edge_cnt[e2] += ways;
                });
            }
            for (int w : touched) cnt[w] = 0;
            touched.clear();
        }
        return res;
    }

    template<typename F>
    void for_each_3_cycle(F&& cb) {
        _build();
        std::vector<int> vis(n + 1, 0), vis_e(n + 1, 0);
        rep(u, 1, n) {
            for (auto [v, e] : dag[u]) {
                vis[v] = u;
                vis_e[v] = e;
            }
            for (auto [v, e_uv] : dag[u]) {
                for (auto [w, e_vw] : dag[v]) {
                    if (vis[w] != u) continue;
                    if constexpr (std::invocable<F, int, int, int, int, int, int>) {
                        cb(u, v, w, e_uv, e_vw, vis_e[w]);
                    } else {
                        cb(u, v, w);
                    }
                }
            }
        }
    }

    i64 count_3_cycles() {
        i64 ans = 0;
        for_each_3_cycle([&](int, int, int) { ++ans; });
        return ans;
    }

    i64 count_4_cycles() { return _solve_4_cycles(false).total; }

    Cycle4Stats count_4_cycles_detail() { return _solve_4_cycles(true); }

    std::vector<int> find_4_cycle() {
        _build();
        std::vector<int> first_v(n + 1, 0), touched, cyc;
        rep(u, 1, n) {
            bool found = _for_2_paths(u, [&](int v, int w, int, int) -> bool {
                if (first_v[w]) {
                    cyc = {u, first_v[w], w, v};
                    return true;
                }
                first_v[w] = v;
                touched.push_back(w);
                return false;
            });
            if (found) return cyc;
            for (int w : touched) first_v[w] = 0;
            touched.clear();
        }
        return {};
    }
};
