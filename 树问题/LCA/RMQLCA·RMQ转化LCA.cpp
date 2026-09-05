#include "aizalib.h"
/*
 * RMQ LCA (RMQ 转化 LCA)
 *
 * Overview:
 *     将 LCA 问题转化为欧拉序上的区间深度 RMQ 问题，通过 ST 表实现 O(1) 在线查询。
 *
 * API:
 *     struct Graph(n)       — 树的邻接表表示，1-based
 *     Graph::add_edge(u, v) — 添加无向边 (u, v)
 *     RMQLCA(G, root = 1)   — 预处理欧拉序与 ST 表，时间复杂度 O(N log N)
 *     lca(u, v)             — 查询 u 和 v 的最近公共祖先，复杂度 O(1)
 *     dist(u, v)            — 查询 u 和 v 的树上距离（边数），复杂度 O(1)
 *
 * Notes:
 *     1. 1-based indexing，root 深度为 0。
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
    std::vector<std::vector<int>> st; // ST表
    int LOG;                          // ST表层数

    RMQLCA(const Graph& G, int root = 1) : G(G), first(G.n + 1), dep(G.n + 1) {
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
                int l = st[i][j-1];
                int r = st[i + (1 << (j-1))][j-1];
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
