#include "aizalib.h"
/*
 * Double Occurrence Euler Tour (双次出现 Euler 序)
 *
 * Overview:
 *     进入 u 和离开 u 时各记录一次节点编号，序列总长度为 2n。
 *     通过按区间出现次数奇偶性翻转（Toggle），支持树上莫队将树上路径查询转化为序列区间查询。
 *
 * API:
 *     DoubleEulerDFS(g, root = 1) — 构造双次出现 Euler 序，以 root 为根预处理
 *     dfs(u, p)                   — 内部 DFS 遍历函数
 *
 * Notes:
 *     1. 1-based indexing；Euler 序有效位置为 1..2n。
 *     2. Time: O(N)；Space: O(N)。
 *     3. 性质: subtree(u) 对应区间 [st[u], ed[u]]。
 *     4. 路径转区间 (设 p = lca(u, v) 且 st[u] <= st[v]):
 *        - 若 p = u，路径 (u, v) 对应 [st[u], st[v]]，出现奇数次的节点恰为路径点集。
 *        - 若 p != u，路径 (u, v) 对应 [ed[u], st[v]]，出现奇数次的节点集加上 p 恰为路径点集。
 */
struct DoubleEulerDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, st, ed, euler;

    DoubleEulerDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), st(n + 1), ed(n + 1), euler(2 * n + 1) {
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
