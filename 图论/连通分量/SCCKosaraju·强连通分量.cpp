#pragma once
#include "aizalib.h"
#include "../0-base[ignore]/Graph·图.cpp"
/*
 * SCCKosaraju·强连通分量
 *
 * Overview:
 *     先在原图上按后序完成第一遍 DFS，再在反图上按逆后序做第二遍 DFS 求强连通分量。
 *     - 正拓扑序特性: 若存在跨分量边 u -> v，必有 scc[u] < scc[v]。
 *
 * API:
 *     SCCKosaraju(g) — 构造并求解强连通分量，复杂度 O(V + E)
 *     scc_cnt        — 强连通分量总数
 *     scc            — scc[u] 为节点 u 所属分量编号 (1-based)
 *     operator[]     — 快捷访问 scc[u]
 *
 * Notes:
 *     1. 1-based indexing.
 *     2. 编号顺序天然为正拓扑序（与 Tarjan 逆拓扑序相反）。
 *
 * Related:
 *     SCCTarjan·强连通分量.cpp::SCC — Tarjan 单遍 DFS 求解器
 */

struct SCCKosaraju {
    int scc_cnt;          // 强连通分量总数
    std::vector<int> scc; // scc[u]: 节点所属分量编号 (1-based)

    SCCKosaraju(const Graph& g) : scc_cnt(0), scc(g.n + 1, 0) {
        std::vector<std::vector<int>> radj(g.n + 1);
        rep(u, 1, g.n) {
            for (int v : g.adj[u]) radj[v].emplace_back(u);
        }

        std::vector<int> order;
        std::vector<char> vis(g.n + 1, 0);
        order.reserve(g.n);

        auto dfs1 = [&](auto&& self, int u) -> void {
            vis[u] = 1;
            for (int v : g.adj[u]) {
                if (!vis[v]) self(self, v);
            }
            order.emplace_back(u);
        };
        auto dfs2 = [&](auto&& self, int u, int id) -> void {
            scc[u] = id;
            for (int v : radj[u]) {
                if (!scc[v]) self(self, v, id);
            }
        };

        rep(u, 1, g.n) {
            if (!vis[u]) dfs1(dfs1, u);
        }
        per(i, (int)order.size() - 1, 0) {
            int u = order[i];
            if (!scc[u]) dfs2(dfs2, u, ++scc_cnt);
        }
    }

    int operator[](int u) const { return scc[u]; }
};
