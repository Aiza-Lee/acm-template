#include "aizalib.h"
/*
 * 点分治 (Centroid Decomposition)
 *
 * Overview:
 *     基于树重心的分治算法框架，用于解决与无向树上路径或距离相关的离线统计问题。
 *     - 树重心与深度保证：对于大小为 S 的无向连通块，其重心 c 满足删除 c
 *       后产生的所有子连通块大小均不超过 S / 2。每次递归选取重心可保证分治树最大深度
 *       严格不超过 log2(N)。
 *     - 路径划分与容斥结构：任意树上路径要么完全落在某个子连通块内部（子问题），
 *       要么跨越当前层重心 c（以 c 为端点或经过 c）。框架通过外部传入的回调函数
 *       calc(c) 统计跨越 c 的贡献，然后将 c 标记为已删除（done[c] =
 *       true）并递归分治所有相邻子块。
 *     - 工具：CentroidDecomposition 结构、solve 分治驱动、get_size 与
 *       get_centroid。
 *
 * API:
 *     CentroidDecomposition(n)   — 初始化包含 n 个节点的点分治结构体（1-based）。
 *     add_edge(u, v)             — 添加树上的无向边 (u, v)。
 *     solve(root, calc)          — 从 root 出发执行分治，对每个分治重心调用
 *                                   calc(c)。
 *     get_size(u, fa = 0)        — 计算包含 u 的未删除连通块大小。
 *     get_centroid(u, fa, total) — 寻找包含 u 且大小为 total 的连通块重心。
 *
 * Notes:
 *     1. 下标统一为 1-based。
 *     2. Time: 分治树深度 O(log N)，框架总时间开销为 O(N log N)；Space: O(N)。
 *     3. 回调函数 calc(int c) 中遍历 c 的邻接表时需检查 !done[v]，
 *        避免跨越已分治节点。
 */

struct CentroidDecomposition {
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<int> siz;
    std::vector<bool> done;

    CentroidDecomposition(int n)
        : n(n), adj(n + 1), siz(n + 1, 0), done(n + 1, false) {}

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    int get_size(int u, int fa = 0) {
        siz[u] = 1;
        for (int v : adj[u]) {
            if (v == fa || done[v]) continue;
            siz[u] += get_size(v, u);
        }
        return siz[u];
    }

    int get_centroid(int u, int fa, int total) {
        for (int v : adj[u]) {
            if (v == fa || done[v]) continue;
            if (siz[v] > total / 2) return get_centroid(v, u, total);
        }
        return u;
    }

    template<class Callback>
    void solve(int root, Callback&& calc) {
        auto decomp = [&](auto&& self, int u) -> void {
            int total = get_size(u, 0);
            int c = get_centroid(u, 0, total);
            calc(c);
            done[c] = true;
            for (int v : adj[c]) {
                if (!done[v]) {
                    self(self, v);
                }
            }
        };
        decomp(decomp, root);
    }
};
