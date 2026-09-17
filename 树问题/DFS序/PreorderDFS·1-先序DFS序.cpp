#include "aizalib.h"
/*
 * 先序 DFS 序 (Preorder DFS Order)
 *
 * Overview:
 *     在树的先序遍历过程中，首次访问节点 u 时记录递增时间戳 pre[u]，并统计子树大小
 *     sz[u]，将整棵树的子树拓扑映射为连续的一维序列区间。
 *     - 映射与区间性质：
 *       1. 子树连续性：以 u 为根的整棵子树在 DFS 序中恰好对应连续闭区间 [pre[u],
 *          pre[u] + sz[u] - 1]。
 *       2. 祖先后代判定：节点 v 属于 u 的子树当且仅当 pre[u] <= pre[v] < pre[u] +
 *          sz[u]。
 *       3. 双向索引：id[pre[u]] = u，支持将树上子树操作直接映射到线段树或树状数组等
 *          连续数据结构。
 *     - 工具：PreorderDFS 求解器，包含 pre、sz、id、fa 等树结构属性。
 *
 * API:
 *     PreorderDFS(g, root = 1) — 以 root 为根预处理树的先序遍历 DFS 序。
 *     dfs(u, p)                — 内部先序遍历推进函数。
 *
 * Notes:
 *     1. 下标统一为 1-based；输入 g 为树的邻接表，默认以 1 为根。
 *     2. Time: O(N)；Space: O(N)。
 */

struct PreorderDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, pre, sz, id;

    PreorderDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), pre(n + 1), sz(n + 1),
          id(n + 1) {
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
