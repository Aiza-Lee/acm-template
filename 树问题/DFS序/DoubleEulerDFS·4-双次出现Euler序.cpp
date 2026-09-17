#include "aizalib.h"
/*
 * 双次出现 Euler 序 (Double Occurrence Euler Tour)
 *
 * Overview:
 *     在进入节点 u 时和离开节点 u 时分别记录一次编号，生成长度为 2n
 *     的括号化欧拉序列，将树上路径对称差映射为一维区间的奇偶翻转（Toggle）。
 *     - 序列与括号结构：每个节点 u 出现且仅出现两次，对应进入位置 st[u] 与离开位置
 *       ed[u]，形成树结构的完全括号匹配。子树 subtree(u) 对应闭区间 [st[u], ed[u]]。
 *     - 路径转区间（树上莫队核心）：设 p = lca(u, v) 且 st[u] <= st[v]：
 *       1. 若 p = u（即 u 为 v 的祖先）：路径 u -> v 对应区间 [st[u], st[v]]，
 *          区间内出现奇数次（即恰好出现 1 次）的节点集恰好构成路径点集。
 *       2. 若 p != u：路径 u -> v 对应区间 [ed[u], st[v]]，
 *          区间内出现奇数次的节点集加上公共祖先 p 恰好构成路径点集。
 *     - 工具：DoubleEulerDFS 结构，包含 st, ed, euler, fa 序列。
 *
 * API:
 *     DoubleEulerDFS(g, root = 1) — 以 root 为根预处理双次出现欧拉序。
 *     dfs(u, p)                   — 内部遍历推进函数。
 *
 * Notes:
 *     1. 下标统一为 1-based；欧拉序有效下标范围为 [1, 2n]。
 *     2. Time: O(N)；Space: O(N)。
 */

struct DoubleEulerDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, st, ed, euler;

    DoubleEulerDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), st(n + 1), ed(n + 1),
          euler(2 * n + 1) {
        dfs(root, 0);
    }

    void dfs(int u, int p) {
        fa[u] = p;
        st[u] = ++timer;
        euler[timer] = u;
        for (int v : g[u]) if (v != p) dfs(v, u);
        ed[u] = ++timer;
        euler[timer] = u;
    }
};
