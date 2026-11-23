#pragma once
#include "aizalib.h"
#include "../0-base[ignore]/Graph·图.cpp"
/*
 * SCCTarjan·强连通分量
 *
 * Overview:
 *     使用 Tarjan 算法求解有向图的强连通分量，仅保留最小拓扑输出 {scc_cnt, scc}，
 *     不驻留工作数组与各分量点集。
 *     - 天然逆拓扑序: 若存在跨分量边 u -> v，必有 scc[u] > scc[v]；按 per(c,
 *       scc_cnt, 1) 遍历即为正拓扑序。
 *
 * API:
 *     SCC(g)     — 构造并求解强连通分量，复杂度 O(V + E)
 *     scc_cnt    — 强连通分量总数
 *     scc        — scc[u] 为点 u 所属分量编号 (1-based)
 *     operator[] — 快捷访问 scc[u]
 *
 * Notes:
 *     1. 1-based indexing. 适用于任意有向图（含重边、自环与非连通图）。
 *     2. 逆拓扑序特性: 拓扑 DP 可直接按 per(c, scc_cnt, 1) 遍历。
 *
 * Related:
 *     - SCCCondensation·强连通缩点.cpp::SCCCondensation — 缩点去重 DAG
 *     - SCCKosaraju·强连通分量.cpp::SCCKosaraju — 两遍 DFS 求解
 *     - TwoSAT·2-SAT可满足性.cpp::TwoSAT — 复用 SCC 求解 2-SAT
 */

struct SCC {
    int scc_cnt;          // 强连通分量总数
    std::vector<int> scc; // scc[u]: 节点所属分量 (1-based)

    SCC(const Graph& g) : scc_cnt(0), scc(g.n + 1, 0) {
        std::vector<int> dfn(g.n + 1, 0), low(g.n + 1, 0), stk;
        std::vector<char> in_stk(g.n + 1, 0);
        stk.reserve(g.n);
        int dfn_cnt = 0;

        auto tarjan = [&](auto&& self, int u) -> void {
            dfn[u] = low[u] = ++dfn_cnt;
            stk.emplace_back(u);
            in_stk[u] = 1;
            for (int v : g.adj[u]) {
                if (!dfn[v]) {
                    self(self, v);
                    low[u] = std::min(low[u], low[v]);
                } else if (in_stk[v]) {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }
            if (low[u] != dfn[u]) return;
            ++scc_cnt;
            while (true) {
                int x = stk.back();
                stk.pop_back();
                in_stk[x] = 0;
                scc[x] = scc_cnt;
                if (x == u) break;
            }
        };

        rep(i, 1, g.n) {
            if (!dfn[i]) tarjan(tarjan, i);
        }
    }

    int operator[](int u) const { return scc[u]; }
};
