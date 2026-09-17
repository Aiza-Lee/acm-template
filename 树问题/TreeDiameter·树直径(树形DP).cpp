#include "aizalib.h"
/*
 * 树直径（树形 DP） (Tree Diameter - Tree DP)
 *
 * Overview:
 *     利用树形动态规划自底向上合并最长向下链求解树的直径，
 *     天然支持任意包含正负边权的带权树。
 *     - 状态设计与链拼合结构：
 *       1. 向下最长链：以 root 定根，对于节点 u，维护其子树内从 u 出发向下的最长链
 *          mx1 与次长链 mx2（必须属于不同的子树分支）。
 *       2. 全局拼合更新：在每个节点 u 处，通过 mx1 + mx2 拼接成以 u
 *          为折返点的树上简单路径，并松弛全局最优直径。
 *       3. 拓扑与迭代优化：基于 BFS 序倒序自底向上转移，无需递归栈，
 *          具备极高的运行速度且防止深度递归爆栈。
 *     - 工具：Graph 结构、TreeDiameterDP 求解器、solve。
 *
 * API:
 *     struct Graph<T>(n)           — 树的带权邻接表表示（1-based）。
 *     Graph::add_edge(u, v, w = 1) — 添加无向边 (u, v, w)。
 *     TreeDiameterDP<T>(G)         — 构造树形 DP 直径求解器。
 *     solve(root = 1)              — 以 root 为定根求解树直径，返回 Result {len, u,
 *                                     v, edges}。
 *
 * Notes:
 *     1. 下标统一为 1-based。
 *     2. Time: O(N)；Space: O(N)。
 */

template<typename T>
concept TreeDiameterWeight =
    std::default_initializable<T> && std::totally_ordered<T> &&
    requires(T a, T b) {
        { a + b } -> std::convertible_to<T>;
    };

template<typename T = i64>
    requires TreeDiameterWeight<T>
struct Graph {
    struct Edge { int v; T w; };
    int n;
    std::vector<std::vector<Edge>> adj;

    Graph(int n) : n(n), adj(n + 1) { AST(n >= 1); }

    void add_edge(int u, int v, T w = T{1}) {
        AST(1 <= u && u <= n && 1 <= v && v <= n);
        AST(!(w < T{}));
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w);
    }
};

template<typename T = i64>
    requires TreeDiameterWeight<T>
struct TreeDiameterDP {
    struct Result {
        T len{};
        int u = 1, v = 1, edges = 0;
    };

    const Graph<T>& G;

    TreeDiameterDP(const Graph<T>& G) : G(G) {}

    Result solve(int root = 1) const {
        AST(1 <= root && root <= G.n);
        std::vector<int> fa(G.n + 1), ord;
        ord.reserve(G.n);
        ord.emplace_back(root);
        rep(i, 0, (int)ord.size() - 1) {
            int u = ord[i];
            for (auto [v, w] : G.adj[u]) if (v != fa[u]) {
                fa[v] = u;
                ord.emplace_back(v);
            }
        }

        Result ans{T{}, root, root, 0};
        std::vector<Down> down(G.n + 1);
        per(i, (int)ord.size() - 1, 0) {
            int u = ord[i];
            Down mx1{T{}, u, 0}, mx2{T{}, u, 0};
            for (auto [v, w] : G.adj[u]) if (fa[v] == u) {
                Down cur{down[v].len + w, down[v].end, down[v].cnt + 1};
                if (_better_down(cur, mx1)) {
                    mx2 = mx1;
                    mx1 = cur;
                } else if (_better_down(cur, mx2)) {
                    mx2 = cur;
                }
            }
            down[u] = mx1;
            _relax(ans, mx1.len + mx2.len, mx1.end, mx2.end, mx1.cnt + mx2.cnt);
        }
        return ans;
    }

private:
    struct Down {
        T len{};
        int end = 1, cnt = 0;
    };

    static bool _better_down(const Down& a, const Down& b) {
        if (b.len < a.len) return true;
        if (a.len < b.len) return false;
        if (a.cnt != b.cnt) return a.cnt > b.cnt;
        return a.end < b.end;
    }

    static void _relax(Result& ans, T len, int u, int v, int edges) {
        if (v < u) std::swap(u, v);
        if (ans.len < len) {
            ans = {len, u, v, edges};
            return;
        }
        if (len < ans.len) return;
        if (edges != ans.edges) {
            if (edges > ans.edges) ans = {len, u, v, edges};
            return;
        }
        if (std::pair{u, v} < std::pair{ans.u, ans.v}) {
            ans = {len, u, v, edges};
        }
    }
};
