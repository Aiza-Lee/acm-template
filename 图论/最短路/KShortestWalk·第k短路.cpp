#include "aizalib.h"
/*
 * K-th Shortest Walk (A*)
 *
 * Overview:
 *     求解有向非负权图从起点 s 到终点 t 的第 k 短游走（允许重复经过点与边）。
 *     首先在反图上以 t 为源点跑 Dijkstra，以各点到 t 的最短路 h[u] 作为 A*
 *     启发函数的精确下界估计（Admissible & Consistent）。
 *     优先队列按估计总代价 f = d + h[u] 顺序搜索，终点第 k 次出队即为第 k 短长度。
 *     对每个节点设置扩展次数上限 k+1，杜绝零权环或自环导致的死循环。
 *     常用于多备用路径规划、受限约束路径搜索及前 k 优状态生成。
 *
 * API:
 *     KShortestWalk(n)              — 建立 1..n 的有向图，O(n)。
 *     add_edge(u, v, w)             — 非负边权，重边独立计数；均摊 O(1)。
 *     solve(s, t, k, include_empty) — k>=1，返回第 k 小长度 optional<i128>；不足 k
 *                                      条返回 nullopt；s=t 的空游走由 include_empty
 *                                      控制。时间 O((n+m)log(n+1)+K*m*log(K*m+2))，
 *                                      空间 O(n+m+K*m)，K<=k+1。
 *
 * Notes:
 *     1. 不是简单路径；同长游走按重数计数，零权环可产生无穷多个等长游走。
 *     2. 每点至多扩展 K 次，避免零环无限扩展；到达 t 后仍可离开并再次到达。
 *     3. 内部距离 i128，避免长游走溢出 i64；solve 不改图，可重复调用。
 */
struct KShortestWalk {
    explicit KShortestWalk(int n) : n(n), adj(n + 1), rev(n + 1) {}
    void add_edge(int u, int v, i64 w) {
        AST(1 <= u && u <= n && 1 <= v && v <= n && w >= 0);
        adj[u].push_back({v, w}); rev[v].push_back({u, w});
    }
    std::optional<i128> solve(int s, int t, int k, bool include_empty = true) const {
        AST(1 <= s && s <= n && 1 <= t && t <= n && k >= 1);
        constexpr i128 INF = (i128)1 << 126;
        using Item = std::pair<i128, int>;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
        std::vector<i128> h(n + 1, INF);
        h[t] = 0; pq.push({0, t});
        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d != h[u]) continue;
            for (auto [v, w] : rev[u]) if (d + w < h[v]) {
                h[v] = d + w; pq.push({h[v], v});
            }
        }
        if (h[s] == INF) return std::nullopt;
        i64 need = (i64)k + (s == t && !include_empty);
        using State = std::tuple<i128, i128, int>;
        std::priority_queue<State, std::vector<State>, std::greater<State>> q;
        std::vector<i64> cnt(n + 1);
        q.push({h[s], 0, s});
        while (!q.empty()) {
            auto [f, d, u] = q.top(); q.pop();
            if (cnt[u] == need) continue;
            ++cnt[u];
            if (u == t && cnt[u] == need) return d;
            for (auto [v, w] : adj[u]) if (h[v] != INF && cnt[v] < need)
                q.push({d + w + h[v], d + w, v});
        }
        return std::nullopt;
    }
private:
    struct Edge { int v; i64 w; };
    int n;
    std::vector<std::vector<Edge>> adj, rev;
};
