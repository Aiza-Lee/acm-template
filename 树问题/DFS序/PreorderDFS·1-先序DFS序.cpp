#include "aizalib.h"
/*
 * Preorder DFS Order (先序 DFS 序)
 *
 * Overview:
 *     首次进入点 u 时记录时间戳 pre[u]，并维护子树大小 sz[u]，将子树点集映射为连续区间。
 *
 * API:
 *     PreorderDFS(g, root = 1) — 构造先序 DFS 序，以 root 为根预处理
 *     dfs(u, p)                — 内部 DFS 遍历函数
 *
 * Notes:
 *     1. 1-based indexing；输入 g 为树的邻接表，默认根为 1。
 *     2. Time: O(N)；Space: O(N)。
 *     3. 性质: subtree(u) 对应区间 [pre[u], pre[u] + sz[u] - 1]。
 *     4. 性质: v 在 u 子树内 <=> pre[u] <= pre[v] < pre[u] + sz[u]。
 *     5. 用法: id[pre[u]] = u，可将子树上的点权直接映射到树状数组或线段树等连续结构上。
 */
struct PreorderDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, pre, sz, id;

    PreorderDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), pre(n + 1), sz(n + 1), id(n + 1) {
        dfs(root, 0);
    }

    void dfs(int u, int p) {
        fa[u] = p;
        pre[u] = ++timer;
        id[timer] = u;
        sz[u] = 1;
        for (int v : g[u]) if (v != p) {
            dfs(v, u);
            sz[u] += sz[v];
        }
    }
};
