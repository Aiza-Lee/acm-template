#pragma once
#include "aizalib.h"
#include "../0-base[ignore]/GraphForwardStar·前向星.cpp"
/*
 * VDCC·点双连通分量
 *
 * Overview:
 *     使用 Tarjan 算法求解无向图的割点与点双连通分量，仅保留最小顶点结构输出
 *     {cut_cnt, is_cut, cuts, v_dcc_cnt, v_dcc}。
 *
 * API:
 *     VDCC(g)   — 构造函数，求解割点与点双连通分量，复杂度 O(V + E)
 *     cut_cnt   — 割点总数
 *     is_cut    — is_cut[u] 标记节点 u 是否为割点
 *     cuts      — 割点列表
 *     v_dcc_cnt — 点双连通分量总数
 *     v_dcc     — v_dcc[c] 为分量 c 包含的节点列表 (1-based)
 *
 * Notes:
 *     1. 1-based indexing. 需使用 GraphFS 链式前向星成对加无向边。
 *     2. 割点可同时属于多个点双连通分量；孤立点自成一个仅含该点的单点分量。
 *
 * Related:
 *     BlockCutTree·圆方树.cpp::BlockCutTree — 圆方树构建与点属性判断
 */

struct VDCC {
    int cut_cnt;                         // 割点总数
    std::vector<bool> is_cut;            // is_cut[u]: 是否为割点
    std::vector<int> cuts;               // 割点列表
    int v_dcc_cnt;                       // 点双连通分量总数
    std::vector<std::vector<int>> v_dcc; // v_dcc[c]: 分量点集 (1-based)

    VDCC(const GraphFS& g)
        : cut_cnt(0), is_cut(g.n + 1, false), cuts(),
          v_dcc_cnt(0), v_dcc(1) {
        std::vector<int> dfn(g.n + 1, 0), low(g.n + 1, 0), stk_v;
        stk_v.reserve(g.n);
        int dfn_cnt = 0, root = 0;

        auto tarjan = [&](auto&& self, int u, int fa_edge) -> void {
            low[u] = dfn[u] = ++dfn_cnt;
            stk_v.emplace_back(u);
            int son_cnt = 0;

            for (int e = g.head[u]; e; e = g.nxt[e]) {
                int v = g.to[e];
                if (!dfn[v]) {
                    ++son_cnt;
                    self(self, v, e);
                    low[u] = std::min(low[u], low[v]);
                    if (low[v] >= dfn[u]) {
                        if ((u != root || son_cnt > 1) && !is_cut[u]) {
                            is_cut[u] = true;
                            cuts.emplace_back(u);
                        }
                        ++v_dcc_cnt;
                        std::vector<int> comp;
                        int x;
                        do {
                            x = stk_v.back();
                            stk_v.pop_back();
                            comp.emplace_back(x);
                        } while (x != v);
                        comp.emplace_back(u);
                        v_dcc.emplace_back(std::move(comp));
                    }
                } else if ((e ^ 1) != fa_edge && dfn[v] < dfn[u]) {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }
            if (u == root && son_cnt == 0) {
                stk_v.pop_back();
                ++v_dcc_cnt;
                v_dcc.emplace_back(std::vector<int>{u});
            }
        };

        rep(i, 1, g.n) {
            if (!dfn[i]) {
                root = i;
                tarjan(tarjan, i, 0);
            }
        }
        cut_cnt = (int)cuts.size();
    }
};
