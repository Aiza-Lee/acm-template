#include "aizalib.h"
#include "../0-base[ignore]/GraphForwardStar·前向星.cpp"

/*
 * EDCC·边双连通分量
 *
 * Overview:
 *     使用 Tarjan 算法求解无向图的割边（桥）与边双连通分量（简单标记版）。
 *
 * API:
 *     EDCC(g)    — 构造函数，求解割边与边双连通分量，复杂度 O(V + E)
 *     bridge_cnt — 桥的总数
 *     is_bridge  — is_bridge[e] 标记边表槽位 e 是否为桥
 *     bridges    — 桥列表，元素为 tuple<int, int, int>(u, v, edge_id)
 *     dcc_cnt    — 边双连通分量个数
 *     dcc        — dcc[u] 为节点 u 所属的边双连通分量编号 (1-based)
 *     sz         — sz[c] 为分量 c 包含的节点数
 *     comp       — comp[c] 为分量 c 包含的节点列表
 *
 * Notes:
 *     1. 1-based indexing. 需使用 GraphFS 链式前向星成对加无向边 (tot 从 1 开始)。
 *     2. edge_id 为 1-based 的无向边插入次序 (e / 2)。
 *     3. 天然支持多重图（重边与自环），重边不会被误判为树边反向边。
 *     4. 若需要构建缩点后的桥森林/缩点树结构，请使用同目录下的
 *        EDCCCondensation·边双缩点.cpp。
 */

struct EDCC {
    GraphFS& g;                                     // 图的引用
    int bridge_cnt;                                 // 桥的总数
    std::vector<bool> is_bridge;                    // 槽位 e 是否为桥
    std::vector<std::tuple<int, int, int>> bridges; // 桥 (u, v, edge_id)
    int dcc_cnt;                                    // 边双分量个数
    std::vector<int> dcc;                           // dcc[u]: 节点所属分量
    std::vector<int> sz;                            // sz[c]: 分量节点数
    std::vector<std::vector<int>> comp;             // comp[c]: 分量点集

    EDCC(GraphFS& g)
        : g(g), bridge_cnt(0), is_bridge(g.to.size(), false),
          bridges(), dcc_cnt(0), dcc(g.n + 1, 0), sz(1, 0), comp(1),
          _dfn(g.n + 1, 0), _low(g.n + 1, 0), _dfn_cnt(0) {
        rep(i, 1, g.n) {
            if (!_dfn[i]) _tarjan(i, 0);
        }
        rep(i, 1, g.n) {
            if (!dcc[i]) {
                ++dcc_cnt;
                comp.emplace_back();
                _dfs_dcc(i, dcc_cnt);
                sz.emplace_back((int)comp[dcc_cnt].size());
            }
        }
    }

private:
    std::vector<int> _dfn, _low;
    int _dfn_cnt;

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

    void _dfs_dcc(int u, int id) {
        dcc[u] = id;
        comp[id].emplace_back(u);
        for (int e = g.head[u]; e; e = g.nxt[e]) {
            if (is_bridge[e]) continue;
            int v = g.to[e];
            if (!dcc[v]) _dfs_dcc(v, id);
        }
    }
};
