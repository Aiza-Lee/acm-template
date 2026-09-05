#include "aizalib.h"
/*
 * Centroid Decomposition (点分治框架)
 *
 * Overview:
 *     树上点分治通用框架。递归寻找当前连通块重心，通过回调函数统计跨重心的路径/子树贡献，
 *     随后标记重心已访问并递归分治各个子连通块。
 *
 * API:
 *     CentroidDecomposition(n)   — 初始化包含 n 个节点的点分治结构体 (1-based)
 *     add_edge(u, v)             — 添加无向边 (u, v)
 *     solve(root, calc)          — 从 root 出发执行点分治，对每层重心调用 calc(c)
 *     get_size(u, fa = 0)        — 计算 u 所在未删除连通块大小
 *     get_centroid(u, fa, total) — 获取包含 u 的连通块重心
 *
 * Notes:
 *     1. 1-based indexing。
 *     2. Time: 框架分治深度 O(log N)，总框架开销 O(N log N)；Space: O(N)。
 *     3. 回调函数 calc(int c) 中可遍历 c 的出边访问未删除邻点（!done[v]），统计跨越 c 的答案。
 */

struct CentroidDecomposition {
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<int> siz;
    std::vector<bool> done;

    CentroidDecomposition(int n) : n(n), adj(n + 1), siz(n + 1, 0), done(n + 1, false) {}

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    int get_size(int u, int fa = 0) {
        siz[u] = 1;
        for (int v : adj[u]) {
            if (v == fa || done[v]) continue;
            siz[u] += get_size(v, u);
        }
        return siz[u];
    }

    int get_centroid(int u, int fa, int total) {
        for (int v : adj[u]) {
            if (v == fa || done[v]) continue;
            if (siz[v] > total / 2) return get_centroid(v, u, total);
        }
        return u;
    }

    template<class Callback>
    void solve(int root, Callback&& calc) {
        auto decomp = [&](auto&& self, int u) -> void {
            int total = get_size(u, 0);
            int c = get_centroid(u, 0, total);
            calc(c);
            done[c] = true;
            for (int v : adj[c]) {
                if (!done[v]) {
                    self(self, v);
                }
            }
        };
        decomp(decomp, root);
    }
};
