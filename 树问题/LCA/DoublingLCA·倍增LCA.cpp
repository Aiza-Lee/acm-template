#include "aizalib.h"
/*
 * 倍增法 LCA (Doubling LCA)
 *
 * Overview:
 *     基于二进制倍增思想的树上最近公共祖先（LCA）在线查询算法，
 *     同时支持两点距离计算与 k 级祖先查询。
 *     - 状态转移与跳跃结构：定义 up[u][i] 为节点 u 向上跳 2^i 步到达的祖先。
 *       满足倍增转移方程 up[u][i] = up[up[u][i - 1]][i - 1]。
 *     - LCA 判定过程：
 *       1. 深度对齐：若两点深度不同，利用二进制拆分将较深节点倍增提升至相同深度。
 *       2. 同步上跳：若重合则已是 LCA；否则两点从最大幂次向下枚举，
 *          同步上跳至不同祖先处，循环结束后其直接父节点 up[u][0] 即为 LCA。
 *     - 工具：Graph 结构、DoublingLCA 求解器、lca、dist、kth_ancestor。
 *
 * API:
 *     struct Graph(n)          — 树的邻接表表示（1-based）。
 *     Graph::add_edge(u, v)    — 添加无向树边 (u, v)。
 *     DoublingLCA(G, root = 1) — 预处理树上倍增表，时间复杂度 O(N log N)。
 *     lca(u, v)                — 查询节点 u 和 v 的最近公共祖先，时间复杂度 O(log
 *                                 N)。
 *     dist(u, v)               — 查询两点树上简单路径边数距离，时间复杂度 O(log
 *                                 N)。
 *     kth_ancestor(u, k)       — 查询节点 u 的第 k 级祖先，越界返回 0，复杂度
 *                                 O(log N)。
 *
 * Notes:
 *     1. 下标统一为 1-based，根节点 root 的深度定义为 0。
 *     2. Time: 预处理 O(N log N)，单次查询 O(log N)；Space: O(N log N)。
 *     3. 倍增表第一维为 log 层（up[LOG][N + 1]），避免 N 次小 vector 堆内存分配与碎片化。
 */

struct Graph {
    int n;
    std::vector<std::vector<int>> adj;
    Graph(int n) : n(n), adj(n + 1) {}
    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

struct DoublingLCA {
    const Graph& G;                   // 图引用
    std::vector<std::vector<int>> up; // up[i][u]: u 的第 2^i 个祖先
    std::vector<int> dep;             // 节点深度
    int LOG;                          // 最大倍增层数

    DoublingLCA(const Graph& G, int root = 1) : G(G), dep(G.n + 1) {
        LOG = std::bit_width((unsigned)G.n);
        up.assign(LOG, std::vector<int>(G.n + 1));
        dfs(root, 0, 0);
    }

    void dfs(int u, int p, int d) {
        dep[u] = d;
        up[0][u] = p;
        rep(i, 1, LOG - 1) up[i][u] = up[i - 1][up[i - 1][u]];
        for (int v : G.adj[u]) {
            if (v != p) dfs(v, u, d + 1);
        }
    }

    int lca(int u, int v) {
        if (dep[u] < dep[v]) std::swap(u, v);
        per(i, LOG - 1, 0) {
            if (dep[u] - (1 << i) >= dep[v]) u = up[i][u];
        }
        if (u == v) return u;
        per(i, LOG - 1, 0) {
            if (up[i][u] != up[i][v]) {
                u = up[i][u];
                v = up[i][v];
            }
        }
        return up[0][u];
    }

    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }

    int kth_ancestor(int u, int k) {
        if (k < 0 || k > dep[u]) return 0;
        rep(i, 0, LOG - 1) {
            if ((k >> i) & 1) u = up[i][u];
        }
        return u;
    }
};
