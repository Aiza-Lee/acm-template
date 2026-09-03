#include "aizalib.h"

/*
 * 树直径（树形 DP，输出有序点集）
 *
 * Overview:
 *      树直径：在带权树中找最长路径（点权和 = 边权和累计）。
 *      本模板沿用树形 DP：以 root 定根，按 BFS 序倒序处理。
 *      每个节点记录最长向下链与最远叶端点，全局最优在每个节点处
 *      合并两条子链取最大。同步记录 parent 链，最终可在 O(N) 重建
 *      直径路径上的有序点集。
 *
 *      相较 TreeDiameter·树直径(树形DP).cpp，本模板额外输出
 *      std::vector<int> path，便于在直径上做点级操作
 *      （点权和、区间染色、统计等）。
 *
 * API:
 *      Graph<T>(n) / add_edge(u, v, w = 1)
 *          建无向树，1-based，边权需非负。
 *      TreeDiameterDPPath<T>(G).solve(root = 1)
 *          求直径与路径，返回
 *              Result { len, u, v, edges, path }
 *          len: 直径长度（边权和）；
 *          u, v: 两端点，u <= v；
 *          edges: 路径边数 == path.size() - 1；
 *          path: u -> v 的有序点集。
 *          复杂度 O(N)。
 *
 * Notes:
 *      1. 1-based indexing；n >= 1。
 *      2. 边权需非负；权重零时退化为最长（边数最多）路径。
 *      3. path.front() == u，path.back() == v；相邻点之间必有边。
 *      4. 平局规则：len 较大者优先；同 len 时边数多者优先；再相同按
 *          端点字典序。路径在所有 tie-break 后做单次重建。
 *
 * Related:
 *      TreeDiameter·树直径(树形DP).cpp: 不输出路径的简化版。
 *      TreeDiameter·树直径(两次DFS).cpp: 两次 BFS/DFS 实现，同样只输出端点。
 */

template<typename T>
concept TreeDiameterWeight = std::default_initializable<T> && std::totally_ordered<T> && requires(T a, T b) {
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