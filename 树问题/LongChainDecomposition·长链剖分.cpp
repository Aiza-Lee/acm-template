#include "aizalib.h"
/*
 * 长链剖分 (Long Chain Decomposition)
 *
 * Overview:
 *     按子树最大深度（树高）选取长儿子，将树剖分为若干互不相交的长链。
 *     结合倍增与链头上下预处理数组，支持 O(1) 在线查询任意点 k 级祖先。
 *     - 剖分机制与链长性质：
 *       1. 长儿子选取：子树最大深度最大的子节点选为长儿子（heavy），
 *          长边相连构成极大长链。
 *       2. 链长定理：若长链链头 h 的长链长度为 len[h]，则以 h 为根的子树最大深度为
 *          len[h]，且 h 向上至少存在 len[h] 级祖先。
 *     - Level Ancestor O(1) 原理：
 *       1. 预处理：对每个链头 h，预处理其向上 len[h]
 *          级祖先（up_nodes）与向下长链上的 len[h] 个节点（down_nodes），
 *          总空间严格为 O(N)。
 *       2. O(1) 定位：查询 u 的第 k 级祖先时，设 2^t <= k < 2^{t + 1}。
 *          利用倍增一步跳 2^t 步到达 x，剩余步数 rem = k - 2^t < 2^t。由于 x
 *          所在长链长度 >= 其子树深度 >= 2^t > rem，因此目标祖先必定落在 x
 *          所在长链链头的向上或向下预处理数组内，实现 O(1) 索引。
 *     - 工具：Graph 结构、LongChainDecomposition 求解器、kth_ancestor、lca、
 *       dist、jump。
 *
 * API:
 *     struct Graph(n)                     — 树的邻接表表示（1-based）。
 *     Graph::add_edge(u, v)               — 添加无向树边 (u, v)。
 *     LongChainDecomposition(G, root = 1) — 预处理长链信息，时间复杂度 O(N log N)。
 *     kth_ancestor(u, k)                  — 在线查询节点 u 的第 k 级祖先，
 *                                            时间复杂度严格 O(1)。
 *     lca(u, v)                           — 查询节点 u 和 v 的最近公共祖先，
 *                                            时间复杂度 O(log N)。
 *     dist(u, v)                          — 查询两点树上边数距离，时间复杂度 O(log
 *                                            N)。
 *     jump(u, v, k)                       — 返回路径 u -> v 上从 u 出发第 k
 *                                            步到达的节点，越界返回 0。
 *
 * Notes:
 *     1. 全部下标均为 1-based，根节点深度定义为 0。
 *     2. Time: 预处理 O(N log N)，kth_ancestor 为 O(1)，lca / dist / jump 为 O(log
 *        N)；Space: O(N log N)。
 *     3. 长链剖分亦常用于树上与深度相关的 DP 空间优化（利用指针复用达到 O(N)
 *        空间）。
 */

struct Graph {
    int n;
    std::vector<std::vector<int>> adj;

    Graph(int n) : n(n), adj(n + 1) {}

    void add_edge(int u, int v) {
        adj[u].emplace_back(v);
        adj[v].emplace_back(u);
    }
};

struct LongChainDecomposition {
    const Graph& G;
    int n, root, LOG;
    std::vector<std::vector<int>> up;
    std::vector<std::vector<int>> up_nodes;   // 链头 h 的第 k 级祖先
    std::vector<std::vector<int>> down_nodes; // 链头 h 向下长链深度 +k 的点
    std::vector<int> fa;
    std::vector<int> dep;
    std::vector<int> len;   // len[u]: 以 u 为起点的最长向下链长（按点数）
    std::vector<int> heavy; // 长儿子
    std::vector<int> top;   // 所在长链链头
    std::vector<int> pos;   // 在所在长链中的位置

    LongChainDecomposition(const Graph& G, int root = 1)
        : G(G), n(G.n), root(root),
          LOG(std::bit_width((unsigned)std::max(1, G.n))),
          up(n + 1, std::vector<int>(LOG)), up_nodes(n + 1),
          down_nodes(n + 1), fa(n + 1), dep(n + 1), len(n + 1),
          heavy(n + 1), top(n + 1), pos(n + 1) {
        _dfs(root, 0);
        _decompose(root, root);
        _build_chain_vectors();
    }

    void _dfs(int u, int p) {
        fa[u] = p;
        up[u][0] = p;
        rep(i, 1, LOG - 1) up[u][i] = up[up[u][i - 1]][i - 1];
        len[u] = 1;
        heavy[u] = 0;
        for (int v : G.adj[u]) {
            if (v == p) continue;
            dep[v] = dep[u] + 1;
            _dfs(v, u);
            if (len[v] + 1 > len[u]) {
                len[u] = len[v] + 1;
                heavy[u] = v;
            }
        }
    }

    void _decompose(int u, int h) {
        top[u] = h;
        pos[u] = down_nodes[h].size();
        down_nodes[h].emplace_back(u);
        if (heavy[u]) _decompose(heavy[u], h);
        for (int v : G.adj[u]) {
            if (v == fa[u] || v == heavy[u]) continue;
            _decompose(v, v);
        }
    }

    void _build_chain_vectors() {
        rep(h, 1, n) if (top[h] == h) {
            up_nodes[h].resize(len[h]);
            int u = h;
            rep(i, 0, len[h] - 1) {
                up_nodes[h][i] = u;
                u = fa[u];
            }
        }
    }

    int _lift(int u, int k) const {
        rep(i, 0, LOG - 1) if ((k >> i) & 1) u = up[u][i];
        return u;
    }

    int kth_ancestor(int u, int k) const {
        if (k < 0 || k > dep[u]) return 0;
        if (k == 0) return u;
        int t = std::bit_width((unsigned)k) - 1;
        int x = up[u][t];
        int rem = k - (1 << t);
        int h = top[x];
        int d = pos[x];
        if (rem <= d) return down_nodes[h][d - rem];
        rem -= d;
        return rem < (int)up_nodes[h].size() ? up_nodes[h][rem] : 0;
    }

    int lca(int u, int v) const {
        if (dep[u] < dep[v]) std::swap(u, v);
        u = _lift(u, dep[u] - dep[v]);
        if (u == v) return u;
        per(i, LOG - 1, 0) if (up[u][i] != up[v][i]) {
            u = up[u][i];
            v = up[v][i];
        }
        return fa[u];
    }

    int dist(int u, int v) const {
        int w = lca(u, v);
        return dep[u] + dep[v] - 2 * dep[w];
    }

    int jump(int u, int v, int k) const {
        int w = lca(u, v);
        int up_len = dep[u] - dep[w];
        int down_len = dep[v] - dep[w];
        int len = up_len + down_len;
        if (k < 0 || k > len) return 0;
        if (k <= up_len) return kth_ancestor(u, k);
        return kth_ancestor(v, len - k);
    }
};
