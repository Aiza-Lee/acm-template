#include "aizalib.h"
#include "Dinic·分层图最大流.cpp"
/*
 * Gomory-Hu Tree
 *
 * Overview:
 *     无向带权图的 Gomory-Hu 最小割树（Gusfield 简易算法实现）。
 *     构造一棵包含原图全部顶点的带权树，满足原图中任意两点 u, v
 *     的最小割容量恰好等于树上 u 到 v 路径上的最小边权。
 *     仅需在原图上调用 n - 1 次 Dinic 最大流即可完成构建；树上预处理倍增 RMQ，
 *     支持全源点对最小割单次 O(log n) 查询。
 *     常用于全点对网络可靠性分析、最小割聚类与多点对瓶颈连通度判定。
 *
 * API:
 *     GomoryHu(n)         — 初始化 1..n 的图，O(n)。
 *     add_edge(u, v, w)   — 无向容量，允许重边/零容量，忽略自环；均摊 O(1)。
 *     build()             — 重建割树与倍增表；O(n*T_flow + n^2 + n log(n+1))，空间
 *                            O(n + m + n log(n+1))。
 *     parent / cut / tree — 根 1，parent[1]=0；cut[v] 为父边容量，tree 为邻接表。
 *     min_cut(u, v)       — build 后 O(log(n+1)) 查询；u=v 返回 i64
 *                            最大值（空路径）。
 *
 * Notes:
 *     1. 每轮新建残量图；添加边后必须重新 build，重复 build 不改变原容量。
 *     2. 原图全部非自环容量之和须 <= LLONG_MAX；Dinic 一般图 T_flow=O(n^2*m)。
 *     3. 空图/单点允许构造与 build；不连通点对最小割为 0；不适用于有向图。
 * Related:
 *     Dinic·分层图最大流.cpp — 最大流核心。
 *     StoerWagner·全局最小割.cpp — 只求单个无向全局最小割。
 */
struct GomoryHu {
    std::vector<int> parent;
    std::vector<i64> cut;
    std::vector<std::vector<std::pair<int, i64>>> tree;

    explicit GomoryHu(int n) : n(n), parent(n + 1), cut(n + 1), tree(n + 1) {}
    void add_edge(int u, int v, i64 w) {
        AST(1 <= u && u <= n && 1 <= v && v <= n && w >= 0);
        built = false;
        if (u != v) edges.push_back({u, v, w});
    }
    void build() {
        std::fill(parent.begin(), parent.end(), 1);
        std::fill(cut.begin(), cut.end(), 0);
        parent[0] = 0; if (n) parent[1] = 0;
        for (auto& a : tree) a.clear();
        rep(s, 2, n) {
            int t = parent[s];
            Dinic<i64> flow(n);
            for (auto [u, v, w] : edges) { flow.add_edge(u, v, w); flow.add_edge(v, u, w); }
            i64 f = flow.solve(s, t);
            std::vector<char> side(n + 1);
            std::deque<int> q{s}; side[s] = 1;
            while (!q.empty()) {
                int u = q.front(); q.pop_front();
                for (int i = flow.g.head[u]; i; i = flow.g.e[i].nxt) {
                    auto e = flow.g.e[i];
                    if (e.w > 0 && !side[e.v]) { side[e.v] = 1; q.push_back(e.v); }
                }
            }
            rep(v, s + 1, n) if (parent[v] == t && side[v]) parent[v] = s;
            if (t != 1 && side[parent[t]]) {
                parent[s] = parent[t]; parent[t] = s;
                cut[s] = cut[t]; cut[t] = f;
            } else cut[s] = f;
        }
        rep(v, 2, n) {
            tree[v].push_back({parent[v], cut[v]});
            tree[parent[v]].push_back({v, cut[v]});
        }
        int lg = std::bit_width((unsigned)std::max(n, 1));
        up.assign(lg, std::vector<int>(n + 1));
        mn.assign(lg, std::vector<i64>(n + 1, std::numeric_limits<i64>::max()));
        dep.assign(n + 1, 0);
        std::deque<int> q; if (n) q.push_back(1);
        while (!q.empty()) {
            int u = q.front(); q.pop_front();
            for (auto [v, w] : tree[u]) if (v != up[0][u]) {
                dep[v] = dep[u] + 1; up[0][v] = u; mn[0][v] = w; q.push_back(v);
            }
        }
        rep(k, 1, lg - 1) rep(v, 1, n) {
            up[k][v] = up[k - 1][up[k - 1][v]];
            mn[k][v] = std::min(mn[k - 1][v], mn[k - 1][up[k - 1][v]]);
        }
        built = true;
    }
    i64 min_cut(int u, int v) const {
        AST(built && 1 <= u && u <= n && 1 <= v && v <= n);
        i64 ans = std::numeric_limits<i64>::max();
        if (dep[u] < dep[v]) std::swap(u, v);
        int d = dep[u] - dep[v];
        rep(k, 0, (int)up.size() - 1) if ((d >> k) & 1) {
            ans = std::min(ans, mn[k][u]); u = up[k][u];
        }
        if (u == v) return ans;
        per(k, (int)up.size() - 1, 0) if (up[k][u] != up[k][v]) {
            ans = std::min({ans, mn[k][u], mn[k][v]}); u = up[k][u]; v = up[k][v];
        }
        return std::min({ans, mn[0][u], mn[0][v]});
    }
private:
    struct Edge { int u, v; i64 w; };
    int n;
    bool built = false;
    std::vector<Edge> edges;
    std::vector<int> dep;
    std::vector<std::vector<int>> up;
    std::vector<std::vector<i64>> mn;
};
