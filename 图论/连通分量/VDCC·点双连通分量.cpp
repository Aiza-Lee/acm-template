#include "aizalib.h"
#include "../0-base[ignore]/GraphForwardStar·前向星.cpp"

/*
 * VDCC·点双连通分量
 *
 * Overview:
 *     使用 Tarjan 算法求解无向图的割点与点双连通分量（简单标记版）。
 *
 * API:
 *     VDCC(g)     — 构造函数，求解割点与点双连通分量，复杂度 O(V + E)
 *     cut_cnt     — 割点总数
 *     is_cut      — is_cut[u] 标记节点 u 是否为割点
 *     cuts        — 割点列表
 *     v_dcc_cnt   — 点双连通分量总数
 *     v_dcc       — v_dcc[c] 为分量 c 包含的节点列表 (1-based)
 *     v_dcc_edges — v_dcc_edges[c] 为分量 c 包含的原图边插入序号列表 (1-based)
 *
 * Notes:
 *     1. 1-based indexing. 需使用 GraphFS 链式前向星成对加无向边。
 *     2. 割点可同时属于多个点双连通分量；孤立点自成一个仅含该点的单点分量。
 *     3. 若需要构建圆方树（Block-Cut Tree）以便树剖/LCA/DP，请使用同目录下的
 *        BlockCutTree·圆方树.cpp。
 */

struct VDCC {
    GraphFS& g;                                 // 图的引用
    int cut_cnt;                                // 割点总数
    std::vector<bool> is_cut;                   // is_cut[u]: 是否为割点
    std::vector<int> cuts;                      // 割点列表
    int v_dcc_cnt;                              // 点双连通分量总数
    std::vector<std::vector<int>> v_dcc;        // v_dcc[c]: 分量点集 (1-based)
    std::vector<std::vector<int>> v_dcc_edges;  // v_dcc_edges[c]: 分量边集 (1-based)

    VDCC(GraphFS& g)
        : g(g), cut_cnt(0), is_cut(g.n + 1, false), cuts(),
          v_dcc_cnt(0), v_dcc(1), v_dcc_edges(1),
          _dfn(g.n + 1, 0), _low(g.n + 1, 0), _dfn_cnt(0), _root(0) {
        rep(i, 1, g.n) {
            if (!_dfn[i]) {
                _root = i;
                _tarjan(i, 0);
            }
        }
        cut_cnt = (int)cuts.size();
    }

private:
    std::vector<int> _dfn, _low;
    std::vector<int> _stk_v, _stk_e;
    int _dfn_cnt;
    int _root;

    void _tarjan(int u, int fa_edge) {
        _low[u] = _dfn[u] = ++_dfn_cnt;
        _stk_v.emplace_back(u);
        int son_cnt = 0;

        for (int e = g.head[u]; e; e = g.nxt[e]) {
            int v = g.to[e];
            if (!_dfn[v]) {
                _stk_e.emplace_back(e);
                ++son_cnt;
                _tarjan(v, e);
                _low[u] = std::min(_low[u], _low[v]);
                if (_low[v] >= _dfn[u]) {
                    if (u != _root || son_cnt > 1) {
                        if (!is_cut[u]) {
                            is_cut[u] = true;
                            cuts.emplace_back(u);
                        }
                    }
                    ++v_dcc_cnt;
                    std::vector<int> comp;
                    std::vector<int> comp_edges;

                    int x;
                    do {
                        x = _stk_v.back();
                        _stk_v.pop_back();
                        comp.emplace_back(x);
                    } while (x != v);
                    comp.emplace_back(u);

                    int edge_id;
                    do {
                        edge_id = _stk_e.back();
                        _stk_e.pop_back();
                        comp_edges.emplace_back(edge_id / 2);
                    } while (edge_id != e);

                    v_dcc.emplace_back(std::move(comp));
                    v_dcc_edges.emplace_back(std::move(comp_edges));
                }
            } else if ((e ^ 1) != fa_edge && _dfn[v] < _dfn[u]) {
                _stk_e.emplace_back(e);
                _low[u] = std::min(_low[u], _dfn[v]);
            }
        }
        if (u == _root && son_cnt == 0) {
            _stk_v.pop_back();
            ++v_dcc_cnt;
            v_dcc.emplace_back(std::vector<int>{u});
            v_dcc_edges.emplace_back(std::vector<int>{});
        }
    }
};
