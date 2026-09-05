#include "aizalib.h"
/*
 * Euler Tour for RMQ-LCA (RMQ-LCA 用 Euler 序)
 *
 * Overview:
 *     DFS 每遍历到一个点就记录一次，从子树回溯时再记录一次父亲，序列总长度为 2n-1。
 *     用于将 LCA 查询转化为序列区间的深度 RMQ 查询。
 *
 * API:
 *     RMQLCAEulerDFS(g, root = 1) — 构造欧拉序，以 root 为根预处理
 *     dfs(u, p)                   — 内部 DFS 遍历函数
 *
 * Notes:
 *     1. 1-based indexing；Euler 序有效位置为 1..2n-1。
 *     2. Time: O(N)；Space: O(N)。
 *     3. 性质: 设 first[u] <= first[v]，则 lca(u, v) 为 euler[first[u]..first[v]] 中深度最小者。
 *     4. 用法: 对 Euler 序按 dep 建立 ST 表，即可实现 O(N log N) 预处理、O(1) 查询 LCA。
 */
struct RMQLCAEulerDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, dep, first, euler;

    RMQLCAEulerDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), dep(n + 1), first(n + 1), euler(2 * n) {
        dfs(root, 0);
    }

    void dfs(int u, int p) {
        fa[u] = p;
        dep[u] = dep[p] + 1;
        euler[++timer] = u;
        if (!first[u]) first[u] = timer;
        for (int v : g[u]) if (v != p) {
            dfs(v, u);
            euler[++timer] = u;
        }
    }
};
