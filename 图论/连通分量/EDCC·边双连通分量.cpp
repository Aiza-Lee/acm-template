#pragma once
#include "aizalib.h"
#include "../0-base[ignore]/GraphForwardStar·前向星.cpp"
/*
 * EDCC·边双连通分量
 *
 * Overview:
 *     使用 Tarjan 算法求解无向图的割边（桥）与边双连通分量，仅保留最小结构输出
 *     {bridge_cnt, is_bridge, bridges, dcc_cnt, dcc}。
 *
 * API:
 *     EDCC(g)    — 构造函数，求解割边与边双连通分量，复杂度 O(V + E)
 *     bridge_cnt — 桥的总数
 *     is_bridge  — is_bridge[e] 标记边表槽位 e 是否为桥
 *     bridges    — 桥列表，元素为 tuple<int, int, int>(u, v, edge_id)
 *     dcc_cnt    — 边双连通分量个数
 *     dcc        — dcc[u] 为节点 u 所属边双连通分量编号 (1-based)
 *     operator[] — 快捷访问 dcc[u]
 *
 * Notes:
 *     1. 1-based indexing. 需使用 GraphFS 链式前向星成对加无向边 (tot 从 1 开始)。
 *     2. edge_id 为 1-based 的无向边插入次序 (e / 2)。
 *     3. 天然支持多重图（重边与自环），重边不会被误判为桥。
 *     4. 缩点建桥树: 遍历 bridges 中的每条边 (u, v, eid)，连接 dcc[u] 与 dcc[v]
 *        即可；叶子节点 (度数 1) 个数 leaf_cnt 对应边双连通化最少加边数 (leaf_cnt +
 *        1) / 2。
 */

struct EDCC {
    int bridge_cnt;                                 // 桥的总数
    std::vector<bool> is_bridge;                    // 槽位 e 是否为桥
    std::vector<std::tuple<int, int, int>> bridges; // 桥 (u, v, edge_id)
    int dcc_cnt;                                    // 边双分量个数
    std::vector<int> dcc;                           // dcc[u]: 节点所属分量

    EDCC(const GraphFS& g)
        : bridge_cnt(0), is_bridge(g.to.size(), false), bridges(),
          dcc_cnt(0), dcc(g.n + 1, 0) {
        std::vector<int> dfn(g.n + 1, 0), low(g.n + 1, 0);
        int dfn_cnt = 0;

        auto tarjan = [&](auto&& self, int u, int from_edge) -> void {
            low[u] = dfn[u] = ++dfn_cnt;
            for (int e = g.head[u]; e; e = g.nxt[e]) {
                if (e == (from_edge ^ 1)) continue;
                int v = g.to[e];
                if (!dfn[v]) {
                    self(self, v, e);
                    low[u] = std::min(low[u], low[v]);
                    if (low[v] > dfn[u]) {
                        is_bridge[e] = is_bridge[e ^ 1] = true;
                        ++bridge_cnt;
                        bridges.emplace_back(u, v, e / 2);
                    }
                } else {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }
        };

        auto dfs_dcc = [&](auto&& self, int u, int id) -> void {
            dcc[u] = id;
            for (int e = g.head[u]; e; e = g.nxt[e]) {
                if (is_bridge[e]) continue;
                int v = g.to[e];
                if (!dcc[v]) self(self, v, id);
            }
        };

        rep(i, 1, g.n) {
            if (!dfn[i]) tarjan(tarjan, i, 0);
        }
        rep(i, 1, g.n) {
            if (!dcc[i]) dfs_dcc(dfs_dcc, i, ++dcc_cnt);
        }
    }

    int operator[](int u) const { return dcc[u]; }
};
