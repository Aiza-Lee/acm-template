#include "aizalib.h"
/*
 * 进出时间戳 (In-Out Timestamp DFS)
 *
 * Overview:
 *     在进入节点 u 时记录递增访问时间戳 tin[u]，在完成其整棵子树遍历后记录最大时间戳
 *     tout[u]，以区间的严格包含关系刻画树上的祖先-后代偏序。
 *     - 祖先关系判定：节点 u 是节点 v 的祖先，当且仅当 tin[u] <= tin[v] 且 tout[v]
 *       <= tout[u]（在有根树中由于进入时间有序，等价于 tin[v] <= tout[u]）。
 *     - 子树区间映射：子树 subtree(u) 对应的节点进入时间恰好落在连续闭区间 [tin[u],
 *       tout[u]]，区间长度 tout[u] - tin[u] + 1 恰好等于子树点数。
 *     - 工具：TinToutDFS 结构，包含 tin、tout、fa 属性。
 *
 * API:
 *     TinToutDFS(g, root = 1) — 以 root 为根预处理进出时间戳。
 *     dfs(u, p)               — 内部遍历推进函数。
 *
 * Notes:
 *     1. 下标统一为 1-based；输入 g 为树的邻接表，默认以 1 为根。
 *     2. Time: O(N)；Space: O(N)。
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
