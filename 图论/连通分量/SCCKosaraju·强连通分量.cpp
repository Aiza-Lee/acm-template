#include "aizalib.h"
#include "../0-base[ignore]/Graph·图.cpp"

/*
 * SCCKosaraju·强连通分量
 *
 * Overview:
 *     先在原图上按后序完成第一遍 DFS，再在反图上按逆后序做第二遍 DFS 求 SCC。
 *
 * API:
 *     SCCKosaraju(g) — 构造时通过两遍 DFS 完成求解，复杂度 O(V + E)
 *     scc_cnt        — 强连通分量个数
 *     scc            — scc[u] 为点 u 所属分量编号 (1-based)
 *     sz             — sz[c] 为分量 c 包含的顶点数
 *     comp           — comp[c] 为分量 c 包含的顶点列表
 *
 * Notes:
 *     1. 1-based indexing. 对不连通图会自动逐块处理。
 *     2. 分量编号顺序取决于第二遍 DFS 的展开顺序。
 */

struct SCCKosaraju {
    Graph& g;                           // 原图引用
    int scc_cnt;                        // 强连通分量个数
    std::vector<int> scc;               // scc[u] 为点 u 所属分量编号
    std::vector<int> sz;                // sz[c] 为分量 c 包含的顶点数
    std::vector<std::vector<int>> comp; // comp[c] 为分量 c 包含的顶点列表
    std::vector<std::vector<int>> radj; // 反图
    std::vector<int> order;             // 第一遍 DFS 后序
    std::vector<char> vis;              // DFS 访问标记

    SCCKosaraju(Graph& g)
        : g(g), scc_cnt(0), scc(g.n + 1), sz(1, 0), comp(1),
          radj(g.n + 1), order(), vis(g.n + 1) {
        order.reserve(g.n);
        rep(u, 1, g.n) {
            for (int v : g.adj[u]) radj[v].emplace_back(u);
        }
        rep(u, 1, g.n) if (!vis[u]) _dfs1(u);
        std::fill(vis.begin(), vis.end(), 0);
        std::reverse(order.begin(), order.end());
        for (int u : order) {
            if (vis[u]) continue;
            ++scc_cnt;
            comp.emplace_back();
            _dfs2(u, scc_cnt);
            sz.emplace_back((int)comp[scc_cnt].size());
        }
    }

private:
    void _dfs1(int u) {
        vis[u] = 1;
        for (int v : g.adj[u]) {
            if (!vis[v]) _dfs1(v);
        }
        order.emplace_back(u);
    }

    void _dfs2(int u, int id) {
        vis[u] = 1;
        scc[u] = id;
        comp[id].emplace_back(u);
        for (int v : radj[u]) {
            if (!vis[v]) _dfs2(v, id);
        }
    }
};
