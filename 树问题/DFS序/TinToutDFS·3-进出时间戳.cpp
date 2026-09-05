#include "aizalib.h"
/*
 * In-Out Timestamp DFS (进出时间戳)
 *
 * Overview:
 *     进入 u 时记录进入时间 tin[u]，遍历完子树后记录子树内最大进入时间 tout[u]。
 *     提供 O(1) 的祖先关系判定与子树区间表达。
 *
 * API:
 *     TinToutDFS(g, root = 1) — 构造进出时间戳，以 root 为根预处理
 *     dfs(u, p)               — 内部 DFS 遍历函数
 *
 * Notes:
 *     1. 1-based indexing；输入 g 为树的邻接表，默认根为 1。
 *     2. Time: O(N)；Space: O(N)。
 *     3. 性质: u 是 v 的祖先 <=> tin[u] <= tin[v] 且 tout[v] <= tout[u]（等价于 tin[v] <= tout[u]）。
 *     4. 性质: subtree(u) 对应时间戳闭区间 [tin[u], tout[u]]，区间长度恰好等于 siz[u]。
 */
struct TinToutDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, tin, tout;

    TinToutDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), tin(n + 1), tout(n + 1) {
        dfs(root, 0);
    }

    void dfs(int u, int p) {
        fa[u] = p;
        tin[u] = ++timer;
        for (int v : g[u]) if (v != p) dfs(v, u);
        tout[u] = timer;
    }
};
