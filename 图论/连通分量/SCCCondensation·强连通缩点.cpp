#pragma once
#include "aizalib.h"
#include "SCCTarjan·强连通分量.cpp"
/*
 * SCCCondensation·强连通缩点
 *
 * Overview:
 *     复用 SCCTarjan 求解强连通分量，并构建缩点后的去重有向无环图（DAG），
 *     提供缩点邻接表、边列表、出入度及正拓扑序。
 *
 * API:
 *     SCCCondensation(g) — 构造并求解 SCC 与缩点 DAG，复杂度 O(V + E log E)
 *     scc_cnt            — 缩点后的分量个数（DAG 顶点数）
 *     scc                — scc[u] 为原图节点 u 所属分量编号 (1-based)
 *     dag                — 缩点 DAG 邻接表（去重边与自环，1-based）
 *     dag_edges          — 缩点图的有向边列表，元素为 pair<int, int>(cu, cv)
 *     ind                — ind[c] 为缩点图节点 c 的入度
 *     outd               — outd[c] 为缩点图节点 c 的出度
 *     topo               — 缩点图的正拓扑序节点序列（长度为 scc_cnt）
 *
 * Notes:
 *     1. 1-based indexing.
 *     2. 有向图强连通化最少加边: scc_cnt <= 1 时为 0，否则为 max(count(ind == 0),
 *        count(outd == 0))。
 *
 * Related:
 *     SCCTarjan·强连通分量.cpp::SCC — 底层强连通分量求解器
 */

struct SCCCondensation {
    int scc_cnt;                                // 缩点 DAG 顶点数
    std::vector<int> scc;                       // 原图节点所属分量
    std::vector<std::vector<int>> dag;          // 缩点 DAG 邻接表
    std::vector<std::pair<int, int>> dag_edges; // 边列表
    std::vector<int> ind;                       // 各节点入度
    std::vector<int> outd;                      // 各节点出度
    std::vector<int> topo;                      // 正拓扑序序列

    SCCCondensation(const Graph& g) {
        SCC core(g);
        scc_cnt = core.scc_cnt;
        scc = std::move(core.scc);
        dag.assign(scc_cnt + 1, {});
        ind.assign(scc_cnt + 1, 0);
        outd.assign(scc_cnt + 1, 0);

        rep(u, 1, g.n) {
            int cu = scc[u];
            for (int v : g.adj[u]) {
                int cv = scc[v];
                if (cu != cv) dag[cu].emplace_back(cv);
            }
        }
        rep(c, 1, scc_cnt) {
            std::sort(dag[c].begin(), dag[c].end());
            dag[c].erase(std::unique(dag[c].begin(), dag[c].end()), dag[c].end());
            outd[c] = (int)dag[c].size();
            for (int to : dag[c]) {
                ++ind[to];
                dag_edges.emplace_back(c, to);
            }
        }
        topo.reserve(scc_cnt);
        per(c, scc_cnt, 1) topo.emplace_back(c);
    }
};
