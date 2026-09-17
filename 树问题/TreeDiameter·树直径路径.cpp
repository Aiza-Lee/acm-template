#include "aizalib.h"
/*
 * 树直径路径（树形 DP） (Tree Diameter Path - Tree DP)
 *
 * Overview:
 *     基于树形动态规划自底向上求解树的直径，并同步记录转移前驱链，
 *     在求得最大路径长度的同时在 O(N) 时间内完整回溯重构直径上的有序节点序列。
 *     - 动态规划与转移追踪：
 *       1. 向下最长链：以 root 定根，按 BFS 序倒序自底向上迭代。
 *          每个节点记录其子树内的最长向下链与其末端叶节点。
 *       2. 全局拼合与折返点：在每个节点 u 处合并两条来自不同子分支的最长链，
 *          若优于当前全局最优，则记录直径端点及折返点 LCA（best_lca = u）。
 *       3. 线性路径重构：通过父节点指针数组 fa，从两端点分别向上遍历至 best_lca，
 *          翻转后拼接生成严格按树上移动顺序排列的节点序列 path（u -> ... -> v）。
 *     - 工具：Graph 结构、TreeDiameterDPPath 求解器、solve。
 *
 * API:
 *     struct Graph<T>(n)           — 树的带权邻接表表示（1-based，边权需非负）。
 *     Graph::add_edge(u, v, w = 1) — 添加无向树边 (u, v, w)。
 *     TreeDiameterDPPath<T>(G)     — 构造树直径与路径求解器。
 *     solve(root = 1)              — 求解直径与有序点序列，返回 Result {len, u, v,
 *                                     edges, path}。
 *
 * Notes:
 *     1. 下标统一为 1-based。
 *     2. 边权需非负；path.front() == u, path.back() == v，相邻节点间必有树边。
 *     3. Time: O(N)；Space: O(N)。
 *
 * Related:
 *     TreeDiameter·树直径(树形DP).cpp: 仅求端点与长度的轻量版。
 *     TreeDiameter·树直径(两次DFS).cpp: 两次遍历实现。
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
struct TreeDiameterDPPath {
    struct Result {
        T len{};
        int u = 1, v = 1;
        int edges = 0;
        std::vector<int> path{1};
    };

    const Graph<T>& G;

    TreeDiameterDPPath(const Graph<T>& G) : G(G) {}

    Result solve(int root = 1) const {
        AST(1 <= root && root <= G.n);
        std::vector<int> fa(G.n + 1, 0), ord;
        ord.reserve(G.n);
        ord.emplace_back(root);
        rep(i, 0, (int)ord.size() - 1) {
            int u = ord[i];
            for (auto [v, w] : G.adj[u]) if (v != fa[u]) {
                fa[v] = u;
                ord.emplace_back(v);
            }
        }

        std::vector<Down> down(G.n + 1);
        Result ans{T{}, root, root, 0, {root}};
        int best_lca = root;
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
            if (_better_result(mx1.len + mx2.len, mx1.end, mx2.end,
                    mx1.cnt + mx2.cnt, ans)) {
                ans.len = mx1.len + mx2.len;
                ans.u = mx1.end;
                ans.v = mx2.end;
                ans.edges = mx1.cnt + mx2.cnt;
                best_lca = u;
            }
        }

        if (ans.u > ans.v) std::swap(ans.u, ans.v);
        ans.path = _rebuild_path(ans.u, ans.v, best_lca, fa);
        ans.edges = (int)ans.path.size() - 1;
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

    static bool _better_result(T len, int u, int v, int edges, const Result& ans) {
        if (ans.len < len) return true;
        if (len < ans.len) return false;
        if (edges != ans.edges) return edges > ans.edges;
        if (u > v) std::swap(u, v);
        int a = ans.u, b = ans.v;
        if (a > b) std::swap(a, b);
        return std::pair(u, v) < std::pair(a, b);
    }

    std::vector<int> _rebuild_path(int pu, int pv, int lca,
        const std::vector<int>& fa) const {
        std::vector<int> asc, desc;
        for (int x = pu; x != lca; x = fa[x]) asc.emplace_back(x);
        asc.emplace_back(lca);
        for (int x = pv; x != lca; x = fa[x]) desc.emplace_back(x);
        desc.emplace_back(lca);
        std::reverse(desc.begin(), desc.end());
        desc.erase(desc.begin());

        std::vector<int> path;
        path.reserve(asc.size() + desc.size());
        for (int x : asc) path.emplace_back(x);
        for (int x : desc) path.emplace_back(x);
        return path;
    }
};