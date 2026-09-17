#include "aizalib.h"
/*
 * RMQ 转化 LCA (RMQ LCA)
 *
 * Overview:
 *     利用欧拉遍历回溯序列将树上 LCA 问题等价转化为静态数组上的区间深度极小值查询
 *     （RMQ），通过 ST 表实现严格 O(1) 的在线 LCA 查询。
 *     - 欧拉序与区间等价性：深度优先遍历树，每次进入节点及从子节点回溯时均将当前节点
 *       追加至序列，得到长度为 2N - 1 的欧拉序列。
 *     - LCA 与深度最小点对应：对于任意两点 u 和 v，
 *       设其在欧拉序列中首次出现的位置分别为 first[u] 与 first[v]（假设 first[u] <=
 *       first[v]）。则两点在树上的简单路径在欧拉序列中完全覆盖在 [first[u],
 *       first[v]] 区间内，且该区间内深度最小（dep
 *       最小）的节点严格等于两点的最近公共祖先。
 *     - 工具：Graph 结构、RMQLCA 求解器、lca、dist。
 *
 * API:
 *     struct Graph(n)       — 树的邻接表表示（1-based）。
 *     Graph::add_edge(u, v) — 添加无向树边 (u, v)。
 *     RMQLCA(G, root = 1)   — 预处理欧拉序与 ST 表，时间复杂度 O(N log N)。
 *     lca(u, v)             — 在线查询节点 u 和 v 的最近公共祖先，时间复杂度严格
 *                              O(1)。
 *     dist(u, v)            — 查询两点树上边数距离，时间复杂度 O(1)。
 *
 * Notes:
 *     1. 下标统一为 1-based，根节点 root 的深度定义为 0。
 *     2. Time: 预处理 O(N log N)，单次查询 O(1)；Space: O(N log N)。
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

struct RMQLCA {
    const Graph& G;                   // 图引用
    std::vector<int> first;           // 节点在欧拉序中首次出现的位置
    std::vector<int> euler;           // 欧拉序序列
    std::vector<int> dep;             // 节点深度
    std::vector<std::vector<int>> st; // ST 表
    int LOG;                          // ST 表层数

    RMQLCA(const Graph& G, int root = 1)
        : G(G), first(G.n + 1), dep(G.n + 1) {
        // euler tour size is 2*N - 1
        euler.reserve(2 * G.n);
        dfs(root, 0, 0);

        int m = euler.size();
        if (m == 0) return;
        LOG = std::bit_width((unsigned)m);
        st.assign(m, std::vector<int>(LOG));

        rep(i, 0, m - 1) st[i][0] = i;

        rep(j, 1, LOG - 1) {
            rep(i, 0, m - (1 << j)) {
                int l = st[i][j - 1];
                int r = st[i + (1 << (j - 1))][j - 1];
                st[i][j] = dep[euler[l]] < dep[euler[r]] ? l : r;
            }
        }
    }

    void dfs(int u, int p, int d) {
        first[u] = euler.size();
        euler.push_back(u);
        dep[u] = d;

        for (int v : G.adj[u]) {
            if (v == p) continue;
            dfs(v, u, d + 1);
            euler.push_back(u);
        }
    }

    int lca(int u, int v) {
        int l = first[u], r = first[v];
        if (l > r) std::swap(l, r);
        int len = r - l + 1;
        int k = std::bit_width((unsigned)len) - 1;

        int x = st[l][k];
        int y = st[r - (1 << k) + 1][k];
        return dep[euler[x]] < dep[euler[y]] ? euler[x] : euler[y];
    }

    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }
};
