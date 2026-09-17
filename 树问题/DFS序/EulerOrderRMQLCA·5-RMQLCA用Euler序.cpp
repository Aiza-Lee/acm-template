#include "aizalib.h"
/*
 * RMQ-LCA 用 Euler 序 (Euler Tour for RMQ-LCA)
 *
 * Overview:
 *     在树的 DFS 遍历过程中，每次初次进入节点及从子节点回溯时均记录当前节点编号，
 *     生成总长度为 2n - 1 的欧拉序列，将 LCA 查询转化为区间深度 RMQ。
 *     - 欧拉回路与区间映射：序列完整记录树边被向下访问与向上回溯的完整轨迹。设节点
 *       u 初次出现的位置为 first[u]。
 *     - RMQ 等价原理：对于任意两点 u, v（假设 first[u] <= first[v]），两点之间的
 *       LCA 必处于遍历序列子区间 [first[u], first[v]] 中，
 *       且恰为该区间内深度最小（dep 最小）的节点。
 *     - 工具：RMQLCAEulerDFS 预处理器，包含 euler、dep、first、fa 数组。
 *
 * API:
 *     RMQLCAEulerDFS(g, root = 1) — 以 root 为根预处理欧拉序与深度数组。
 *     dfs(u, p)                   — 内部遍历推进函数。
 *
 * Notes:
 *     1. 下标统一为 1-based；欧拉序有效下标范围为 [1, 2n - 1]。
 *     2. Time: O(N)；Space: O(N)。
 *     3. 后续搭配 ST 表可实现 O(N log N) 预处理、O(1) 无回溯查询 LCA。
 */

struct RMQLCAEulerDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, dep, first, euler;

    RMQLCAEulerDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), dep(n + 1), first(n + 1),
          euler(2 * n) {
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
