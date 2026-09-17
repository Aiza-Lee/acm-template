#include "aizalib.h"
/*
 * 树直径（两次 DFS） (Tree Diameter - Two DFS)
 *
 * Overview:
 *     在非负边权无向树上利用两阶段最远点搜索快速求出树的直径（最长简单路径）。
 *     - 几何贪心性质：在非负权树中，从任意起点出发能够到达的最远节点必定是某条直径的
 *       端点之一。
 *     - 两阶段遍历结构：
 *       1. 第一阶段：从任意起点 s 出发遍历，找到距离 s 最远的节点 a。
 *       2. 第二阶段：从节点 a 出发再次遍历，找到距离 a 最远的节点 b。
 *       3. 则路径 a -> b 即为树上的一条直径，其长度为 dis(a, b)。
 *     - 工具：Graph 结构、TreeDiameterTwoDFS 求解器、solve。
 *
 * API:
 *     struct Graph<T>(n)           — 树的带权邻接表表示（1-based，边权需非负）。
 *     Graph::add_edge(u, v, w = 1) — 添加非负权无向边 (u, v, w)。
 *     TreeDiameterTwoDFS<T>(G)     — 构造树直径求解器。
 *     solve(s = 1)                 — 从 s 出发两阶段求解直径，返回 Result {len, u,
 *                                     v, edges}。
 *
 * Notes:
 *     1. 下标统一为 1-based。
 *     2. 边权必须非负；若存在负边权，必须改用树形 DP 算法。
 *     3. Time: O(N)；Space: O(N)。
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
            if (_better_far(u, dis[u], dep[u], best)) {
                best = {u, dep[u], dis[u]};
            }
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
