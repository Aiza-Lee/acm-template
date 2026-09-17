#include "aizalib.h"
/*
 * 重链剖分 LCA (HLD LCA)
 *
 * Overview:
 *     基于轻重链剖分（Heavy-Light Decomposition）的树上最近公共祖先在线查询算法。
 *     - 剖分与链顶结构：每个非叶节点选择子树规模最大的子节点作为重儿子，
 *       由重边相连形成若干极大重链。每个节点属于且仅属于一条重链，
 *       其所在链最高点记为 top[u]。
 *     - 链顶跳转与跳步上限：
 *       1. 链顶跳转：比较两点所在重链链顶深度 dep[top[u]] 与 dep[top[v]]，
 *          将链顶较深的点一步跃升至链顶的父节点 fa[top[u]]。
 *       2. 终止判断：当两点落入同一条重链（top[u] == top[v]）时，深度较浅者即为
 *          LCA。
 *       3. 复杂度与常数优势：由轻边性质知，任意节点到根路径上至多经过 log2(N)
 *          条轻边，因此链跳转至多 O(log N) 次，且实际运行中常数极小，
 *          预处理严格线性。
 *     - 工具：Graph 结构、HLDLCA 求解器、lca、dist。
 *
 * API:
 *     struct Graph(n)       — 树的邻接表表示（1-based）。
 *     Graph::add_edge(u, v) — 添加无向树边 (u, v)。
 *     HLDLCA(G, root = 1)   — 两次 DFS 线性预处理重链信息，时间复杂度 O(N)。
 *     lca(u, v)             — 查询节点 u 和 v 的最近公共祖先，最坏时间复杂度 O(log
 *                              N)。
 *     dist(u, v)            — 查询两点树上边数距离，最坏时间复杂度 O(log N)。
 *
 * Notes:
 *     1. 下标统一为 1-based，根节点 root 的深度定义为 0。
 *     2. Time: 预处理 O(N)，单次查询最坏 O(log N) 且常数显著优于倍增；Space: O(N)。
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

struct HLDLCA {
    const Graph& G;       // 图引用
    std::vector<int> fa;  // 父节点
    std::vector<int> dep; // 深度
    std::vector<int> siz; // 子树大小
    std::vector<int> son; // 重儿子
    std::vector<int> top; // 所在重链顶端节点

    HLDLCA(const Graph& G, int root = 1)
        : G(G), fa(G.n + 1), dep(G.n + 1), siz(G.n + 1), son(G.n + 1),
          top(G.n + 1) {
        dfs1(root, 0, 0);
        dfs2(root, root);
    }

    void dfs1(int u, int p, int d) {
        dep[u] = d; fa[u] = p; siz[u] = 1; son[u] = 0;
        int max_siz = -1;
        for (int v : G.adj[u]) {
            if (v == p) continue;
            dfs1(v, u, d + 1);
            siz[u] += siz[v];
            if (siz[v] > max_siz) {
                max_siz = siz[v];
                son[u] = v;
            }
        }
    }

    void dfs2(int u, int t) {
        top[u] = t;
        if (!son[u]) return;
        dfs2(son[u], t);
        for (int v : G.adj[u]) {
            if (v != fa[u] && v != son[u]) {
                dfs2(v, v);
            }
        }
    }

    int lca(int u, int v) {
        while (top[u] != top[v]) {
            if (dep[top[u]] < dep[top[v]]) std::swap(u, v);
            u = fa[top[u]];
        }
        return dep[u] < dep[v] ? u : v;
    }

    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }
};
