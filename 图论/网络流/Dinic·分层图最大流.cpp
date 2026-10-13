#include "aizalib.h"

/*
 * Dinic·分层图最大流
 *
 * Overview:
 *     Dinic's Algorithm 基于分层图的最大流算法，每轮先 BFS 建立分层图，再在分层图上
 *     DFS 多路增广。
 *
 * API:
 *     Dinic(int n, int m = 0)                  — 初始化 n 个点、预估 m 条原图边的网络
 *     int add_edge(int u, int v, Cap w)        — 添加一条容量为 w 的有向边，返回正向边编号 id (从 2 开始)
 *     Cap solve(int s, int t, Cap limit = INF) — 返回至多增广 limit 流量后的最大流
 *     Cap edge_flow(int id) const              — 读取正向边 id 的实际流量 (等于反向边当前残量)
 *
 * Notes:
 *     模板参数: Cap (容量类型)
 *     1. Time: 一般图 O(V^2E)；单位容量网络 O(min(V^(2/3), E^(1/2))E)；二分图匹配
 *        O(E\sqrt{V})
 *     2. Space: O(V + E)
 *     3. 1-based indexing (点编号 1~n).
 *     4. 方案构造与最小割:
 *        4.1 流量方案: add_edge 返回正向边编号 id (偶数，从 2 开始)；调用 solve 后
 *            可通过 edge_flow(id) 读取该边实际推流值 (底层即 g.e[id ^ 1].w)。
 *        4.2 最小割方案: solve 结束后残量网络中 dep[u] > 0 的点属于源割集 S；
 *            若正向边满足 dep[u] > 0 且 dep[v] == 0，则该边为最小割的饱和割边。
 *     5. 用法/技巧: 二分图最大匹配可直接建模求解；部分流量限制可传入 solve(s, t, limit)。
 */

template<typename Cap>
struct Graph {
    struct Edge {
        int v, nxt;
        Cap w;
    };

    int n;                 // 点数
    std::vector<int> head; // 链式前向星表头
    std::vector<Edge> e;   // 残量网络边集
    int ec;                // 当前边计数，边下标从 2 开始，便于 i ^ 1 找反边

    Graph(int n, int m = 0)
        : n(n), head(n + 1, 0), e(std::max(2 * m + 2, 2)), ec(1) {}

    int add_edge(int u, int v, Cap w) {
        AST(1 <= u && u <= n);
        AST(1 <= v && v <= n);
        if (ec + 2 >= (int)e.size()) e.resize(std::max((int)e.size() * 2, ec + 3));
        e[++ec] = {v, head[u], w};
        head[u] = ec;
        e[++ec] = {u, head[v], 0};
        head[v] = ec;
        return ec - 1;
    }
};

template<typename Cap = i64>
struct Dinic {
    static constexpr Cap INF = std::numeric_limits<Cap>::max();

    Graph<Cap> g;
    std::vector<int> dep; // dep[u]: 分层图中点 u 的层数，0 表示当前不可达
    std::vector<int> cur; // cur[u]: 当前弧优化指针
    int n;                // 点数

    Dinic(int n, int m = 0) : g(n, m), dep(n + 1), cur(n + 1), n(n) {}

    int add_edge(int u, int v, Cap w) {
        return g.add_edge(u, v, w);
    }

    Cap edge_flow(int id) const {
        return g.e[id ^ 1].w;
    }

    bool bfs(int s, int t) {
        std::fill(dep.begin(), dep.end(), 0);
        std::deque<int> q;
        dep[s] = 1;
        q.push_back(s);

        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            for (int i = g.head[u]; i; i = g.e[i].nxt) {
                auto& e = g.e[i];
                if (e.w == 0 || dep[e.v]) continue;
                dep[e.v] = dep[u] + 1;
                if (e.v == t) return true;
                q.push_back(e.v);
            }
        }
        return dep[t] != 0;
    }

    Cap dfs(int u, int t, Cap flow) {
        if (u == t || flow == 0) return flow;

        Cap used = 0;
        for (int& i = cur[u]; i; i = g.e[i].nxt) {
            auto& e = g.e[i];
            if (e.w == 0 || dep[e.v] != dep[u] + 1) continue;
            Cap pushed = dfs(e.v, t, std::min(flow - used, e.w));
            if (pushed == 0) continue;
            e.w -= pushed;
            g.e[i ^ 1].w += pushed;
            used += pushed;
            if (used == flow) return used;
        }
        if (used == 0) dep[u] = 0; // Pruning
        return used;
    }

    Cap solve(int s, int t, Cap limit = INF) {
        AST(1 <= s && s <= n);
        AST(1 <= t && t <= n);
        if (s == t || limit == 0) return 0;

        Cap max_flow = 0;
        while (max_flow < limit && bfs(s, t)) {
            std::copy(g.head.begin(), g.head.end(), cur.begin());
            max_flow += dfs(s, t, limit - max_flow);
        }
        return max_flow;
    }
};
