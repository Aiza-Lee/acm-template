#include "aizalib.h"
/*
 * Postorder DFS Order (后序 DFS 序)
 *
 * Overview:
 *     处理完全部子节点后记录时间戳 post[u]，保证任意后代节点的离开时间戳均早于其祖先。
 *
 * API:
 *     PostorderDFS(g, root = 1) — 构造后序 DFS 序，以 root 为根预处理
 *     dfs(u, p)                 — 内部 DFS 遍历函数
 *
 * Notes:
 *     1. 1-based indexing；输入 g 为树的邻接表，默认根为 1。
 *     2. Time: O(N)；Space: O(N)。
 *     3. 性质: 若 v 是 u 的真后代，则 post[v] < post[u]。
 *     4. 用法: 适合按子树后序遍历进行自底向上的树形 DP 或拓扑推导。
 */
struct PostorderDFS {
    const std::vector<std::vector<int>>& g;
    int n, timer = 0;
    std::vector<int> fa, post;

    PostorderDFS(const std::vector<std::vector<int>>& g, int root = 1)
        : g(g), n((int)g.size() - 1), fa(n + 1), post(n + 1) {
        dfs(root, 0);
    }

    void dfs(int u, int p) {
        fa[u] = p;
        for (int v : g[u]) if (v != p) dfs(v, u);
        post[u] = ++timer;
    }
};
