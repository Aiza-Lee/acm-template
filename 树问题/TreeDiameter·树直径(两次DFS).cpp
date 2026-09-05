#include "aizalib.h"

/*
 * Tree Diameter (树直径，两次 DFS)
 *
 * Overview:
 *     在带权树上任取起点找最远点 a，再从 a 找最远点 b，a-b 即为一条直径。
 *     适用于非负边权树。
 *
 * API:
 *     struct Graph<T>(n)               — 建无向树，1-based，边权需非负
 *     Graph::add_edge(u, v, w = 1)     — 添加无向边 (u, v, w)
 *     struct TreeDiameterTwoDFS<T>(G)  — 树直径求解器
 *     TreeDiameterTwoDFS::solve(s = 1) — 求直径，返回 Result {len, u, v, edges}
 *
 * Notes:
 *     1. 1-based indexing，边权必须非负（若存在负权边需使用树形 DP）。
 *     2. Time: O(N)；Space: O(N)。
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
struct TreeDiameterTwoDFS {
    struct Result {
        T len{};
        int u = 1, v = 1, edges = 0;
    };

    const Graph<T>& G;

    TreeDiameterTwoDFS(const Graph<T>& G) : G(G) {}

    Result solve(int s = 1) const {
        AST(1 <= s && s <= G.n);
        auto a = _farthest(s);
        auto b = _farthest(a.u);
        return {b.dis, a.u, b.u, b.dep};
    }

private:
    struct Far {
        int u, dep;
        T dis;
    };

    static bool _better_far(int u, T du, int dep, const Far& best) {
        if (best.dis < du) return true;
        if (du < best.dis) return false;
        if (dep != best.dep) return dep > best.dep;
        return u < best.u;
    }

    Far _farthest(int s) const {
        std::vector<int> fa(G.n + 1), dep(G.n + 1), stk{s};
        std::vector<T> dis(G.n + 1);
        Far best{s, 0, T{}};
        while (!stk.empty()) {
            int u = stk.back();
            stk.pop_back();
            if (_better_far(u, dis[u], dep[u], best)) best = {u, dep[u], dis[u]};
            for (auto [v, w] : G.adj[u]) if (v != fa[u]) {
                fa[v] = u;
                dep[v] = dep[u] + 1;
                dis[v] = dis[u] + w;
                stk.emplace_back(v);
            }
        }
        return best;
    }
};
