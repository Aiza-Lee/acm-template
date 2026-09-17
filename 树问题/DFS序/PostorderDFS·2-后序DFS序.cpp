#include "aizalib.h"
/*
 * 后序 DFS 序 (Postorder DFS Order)
 *
 * Overview:
 *     在遍历完节点 u 的所有子树后记录其离开时间戳 post[u]，
 *     建立满足自底向上依赖关系的拓扑时序。
 *     - 拓扑依赖性质：对于树上的任意节点 u 及其真后代 v，必有 post[v] < post[u]。
 *       这保证了在按 post[u] 升序遍历时，子树的所有信息均已被前置计算完毕。
 *     - 树形 DP 与自底向上规约：适用于无需递归的循环式树形 DP、
 *       拓扑推导与消除深层递归栈溢出的后序遍历重构。
 *     - 工具：PostorderDFS 求解器，包含 post、fa 等树结构属性。
 *
 * API:
 *     PostorderDFS(g, root = 1) — 以 root 为根预处理树的后序遍历 DFS 序。
 *     dfs(u, p)                 — 内部后序遍历推进函数。
 *
 * Notes:
 *     1. 下标统一为 1-based；输入 g 为树的邻接表，默认以 1 为根。
 *     2. Time: O(N)；Space: O(N)。
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
