#include "aizalib.h"
/*
 * Doubling LCA (倍增法 LCA)
 *
 * Overview:
 *     基于二进制倍增的树上 LCA 在线算法。支持查询任意两点 LCA、树上距离及 k 级祖先。
 *
 * API:
 *     struct Graph(n)          — 树的邻接表表示，1-based
 *     Graph::add_edge(u, v)    — 添加无向边 (u, v)
 *     DoublingLCA(G, root = 1) — 预处理倍增数组，时间复杂度 O(N log N)
 *     lca(u, v)                — 查询 u 和 v 的最近公共祖先，复杂度 O(log N)
 *     dist(u, v)               — 查询 u 和 v 的树上距离（边数），复杂度 O(log N)
 *     kth_ancestor(u, k)       — 查询 u 的第 k 级祖先，越界返回 0，复杂度 O(log N)
 *
 * Notes:
 *     1. 1-based indexing，root 深度为 0。
 *     2. Time: 预处理 O(N log N)，单次查询 O(log N)；Space: O(N log N)。
 */

struct Graph {
    int n;
    std::vector<std::vector<int>> adj;
    Graph(int n) : n(n), adj(n + 1) {}
    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

struct DoublingLCA {
    const Graph& G;                   // 图引用
    std::vector<std::vector<int>> up; // up[u][i]: u的第2^i个祖先
    std::vector<int> dep;             // 节点深度
    int LOG;                          // 最大倍增层数

    DoublingLCA(const Graph& G, int root = 1) : G(G), dep(G.n + 1) {
        LOG = std::bit_width((unsigned)G.n);
        up.assign(G.n + 1, std::vector<int>(LOG));
        dfs(root, 0, 0);
    }

    void dfs(int u, int p, int d) {
        dep[u] = d;
        up[u][0] = p;
        rep(i, 1, LOG - 1) up[u][i] = up[up[u][i-1]][i-1];
        for (int v : G.adj[u]) {
            if (v != p) dfs(v, u, d + 1);
        }
    }

    int lca(int u, int v) {
        if (dep[u] < dep[v]) std::swap(u, v);
        per(i, LOG - 1, 0) {
            if (dep[u] - (1 << i) >= dep[v]) u = up[u][i];
        }
        if (u == v) return u;
        per(i, LOG - 1, 0) {
            if (up[u][i] != up[v][i]) {
                u = up[u][i];
                v = up[v][i];
            }
        }
        return up[u][0];
    }

    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }

    int kth_ancestor(int u, int k) {
        if (k < 0 || k > dep[u]) return 0;
        rep(i, 0, LOG - 1) {
            if ((k >> i) & 1) u = up[u][i];
        }
        return u;
    }
};
