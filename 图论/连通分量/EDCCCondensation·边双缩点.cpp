#include "aizalib.h"
#include "../0-base[ignore]/GraphForwardStar·前向星.cpp"

/*
 * EDCCCondensation·边双缩点
 *
 * Overview:
 *     使用 Tarjan 算法求解边双连通分量，并将每个边双缩为一个点，
 *     割边（桥）作为树边构建缩点树/森林（Bridge Tree / Bridge Forest）。
 *
 * API:
 *     EDCCCondensation(g) — 构造并求解边双与缩点树，复杂度 O(V + E)
 *     dcc_cnt             — 边双分量总数（缩点树的顶点数）
 *     dcc                 — dcc[u] 为原图节点 u 所属分量编号 (1-based)
 *     sz                  — sz[c] 为分量 c 包含的原图节点数
 *     comp                — comp[c] 为分量 c 包含的原图节点列表
 *     bridge_cnt          — 桥的总数（缩点森林的边数）
 *     is_bridge           — is_bridge[e] 标记边槽位 e 是否为桥
 *     bridges             — 桥列表，元素为 tuple<int, int, int>(u, v, edge_id)
 *     tree_adj            — 缩点树邻接表，元素为 pair<int, int>(to_dcc, edge_id)
 *     tree_edges          — 缩点树边列表，元素为 tuple<int, int, int>(cu, cv,
 *                            edge_id)
 *     deg                 — deg[c] 为缩点树中节点 c 的度数
 *
 * Notes:
 *     1. 1-based indexing. 原图连通时缩点图为单棵树，非连通时为森林。
 *     2. tree_adj 存储出边目标分量及对应原图无向边插入序号 edge_id (1-based)。
 *     3. 缩点树叶子节点数 (deg[c] == 1) 对应常见加边问题：使图成为边双连通图所需最少
 *        添加边数为 (leaf_cnt + 1) / 2。
 *     4. 若仅需割边与分量标记而无需构建缩点树，请使用同目录下的
 *        EDCC·边双连通分量.cpp。
 */

struct EDCCCondensation {
    GraphFS& g;                                     // 原图引用
    int dcc_cnt;                                    // 边双分量个数
    std::vector<int> dcc;                           // dcc[u]: 节点所属分量
    std::vector<int> sz;                            // sz[c]: 分量大小
    std::vector<std::vector<int>> comp;             // comp[c]: 分量点集
    int bridge_cnt;                                 // 桥的总数
    std::vector<bool> is_bridge;                    // 槽位是否为桥
    std::vector<std::tuple<int, int, int>> bridges; // 桥 (u, v, edge_id)
    // 缩点树邻接表 (to_dcc, edge_id) 与边列表 (cu, cv, edge_id)
    std::vector<std::vector<std::pair<int, int>>> tree_adj;
    std::vector<std::tuple<int, int, int>> tree_edges;
    std::vector<int> deg;                           // 缩点树节点度数

    EDCCCondensation(GraphFS& g)
        : g(g), dcc_cnt(0), dcc(g.n + 1, 0), sz(1, 0), comp(1),
          bridge_cnt(0), is_bridge(g.to.size(), false), bridges(),
          tree_adj(), tree_edges(), deg() {
        _find_bridges();
        _find_dcc();
        _build_tree();
    }

private:
    std::vector<int> _dfn, _low;
    int _dfn_cnt = 0;

    void _find_bridges() {
        _dfn.assign(g.n + 1, 0);
        _low.assign(g.n + 1, 0);
        rep(i, 1, g.n) {
            if (!_dfn[i]) _tarjan(i, 0);
        }
    }

    void _tarjan(int u, int from_edge) {
        _low[u] = _dfn[u] = ++_dfn_cnt;
        for (int e = g.head[u]; e; e = g.nxt[e]) {
            int v = g.to[e];
            if (e == (from_edge ^ 1)) continue;
            if (!_dfn[v]) {
                _tarjan(v, e);
                _low[u] = std::min(_low[u], _low[v]);
                if (_low[v] > _dfn[u]) {
                    is_bridge[e] = true;
                    is_bridge[e ^ 1] = true;
                    ++bridge_cnt;
                    bridges.emplace_back(u, v, e / 2);
                }
            } else {
                _low[u] = std::min(_low[u], _dfn[v]);
            }
        }
    }

    void _find_dcc() {
        rep(i, 1, g.n) {
            if (!dcc[i]) {
                ++dcc_cnt;
                comp.emplace_back();
                _dfs_dcc(i, dcc_cnt);
                sz.emplace_back((int)comp[dcc_cnt].size());
            }
        }
    }

    void _dfs_dcc(int u, int id) {
        dcc[u] = id;
        comp[id].emplace_back(u);
        for (int e = g.head[u]; e; e = g.nxt[e]) {
            if (is_bridge[e]) continue;
            int v = g.to[e];
            if (!dcc[v]) _dfs_dcc(v, id);
        }
    }

    void _build_tree() {
        tree_adj.assign(dcc_cnt + 1, {});
        deg.assign(dcc_cnt + 1, 0);
        tree_edges.reserve(bridge_cnt);
        for (auto [u, v, eid] : bridges) {
            int cu = dcc[u], cv = dcc[v];
            tree_adj[cu].emplace_back(cv, eid);
            tree_adj[cv].emplace_back(cu, eid);
            tree_edges.emplace_back(cu, cv, eid);
            ++deg[cu];
            ++deg[cv];
        }
    }
};
