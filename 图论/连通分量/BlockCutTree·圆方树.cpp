#include "aizalib.h"
#include "../0-base[ignore]/GraphForwardStar·前向星.cpp"

/*
 * BlockCutTree·圆方树
 *
 * Overview:
 *     使用 Tarjan 算法求解点双连通分量，并构建圆方树/森林（Block-Cut Tree）。
 *     原图节点作为圆点 (1..n)，各个点双连通分量作为方点 (n+1..n+v_dcc_cnt)。
 *
 * API:
 *     BlockCutTree(g)  — 构造函数，求解点双并建立圆方树，复杂度 O(V + E)
 *     n                — 原图节点数（圆点个数）
 *     v_dcc_cnt        — 点双连通分量总数（方点个数）
 *     bct_n            — 圆方树总节点数 (n + v_dcc_cnt)
 *     is_cut           — is_cut[u] 原图节点 u 是否为割点
 *     cuts             — 割点列表
 *     v_dcc            — v_dcc[c] 点双连通分量 c 包含的原图节点列表 (1-based)
 *     tree_adj         — 圆方树邻接表 (大小 bct_n + 1，圆点与方点无向连边)
 *     tree_edges       — 圆方树无向边列表，元素为 pair<int, int>(round_u, square_v)
 *     is_round(u)      — 判断节点 u 是否为圆点 (1 <= u <= n)
 *     is_square(u)     — 判断节点 u 是否为方点 (n < u <= bct_n)
 *     square_id(c)     — 获取点双分量 c 对应的方点编号 (n + c)
 *     dcc_of_square(s) — 获取方点 s 对应的点双分量编号 (s - n)
 *
 * Notes:
 *     1. 1-based indexing. 原图连通时圆方树为树，非连通时为森林。
 *     2. 圆方树是二分图：边仅存在于圆点和方点之间，不存在圆-圆或方-方边。
 *     3. 割点判定：度数大于 1 的圆点即为原图中的割点。
 *     4. 路径性质：原图两点 u, v 间所有简单路径的必须经过点（割点），对应圆方树 u
 *        到 v 树上路径中的所有圆点；路径覆盖的所有点双对应路径中的所有方点。
 *     5. 若仅需割点与点双分量标记而无需构建圆方树，请使用同目录下的
 *        VDCC·点双连通分量.cpp。
 */

struct BlockCutTree {
    GraphFS& g;                                 // 原图引用
    int n;                                      // 圆点个数
    int v_dcc_cnt;                              // 方点个数（点双分量数）
    int bct_n;                                  // 圆方树总点数
    std::vector<bool> is_cut;                   // is_cut[u]: 原图点是否为割点
    std::vector<int> cuts;                      // 割点列表
    std::vector<std::vector<int>> v_dcc;        // 各点双包含的原图点 (1-based)
    std::vector<std::vector<int>> tree_adj;     // 圆方树邻接表 (1..bct_n)
    std::vector<std::pair<int, int>> tree_edges; // 圆方树边列表 (round, square)

    BlockCutTree(GraphFS& g)
        : g(g), n(g.n), v_dcc_cnt(0), bct_n(0),
          is_cut(g.n + 1, false), cuts(), v_dcc(1),
          tree_adj(), tree_edges() {
        _find_v_dcc();
        _build_tree();
    }

    bool is_round(int u) const {
        return 1 <= u && u <= n;
    }

    bool is_square(int u) const {
        return n < u && u <= bct_n;
    }

    int square_id(int c) const {
        return n + c;
    }

    int dcc_of_square(int s) const {
        return s - n;
    }

private:
    std::vector<int> _dfn, _low, _stk_v;
    int _dfn_cnt = 0;
    int _root = 0;

    void _find_v_dcc() {
        _dfn.assign(n + 1, 0);
        _low.assign(n + 1, 0);
        rep(i, 1, n) {
            if (!_dfn[i]) {
                _root = i;
                _tarjan(i, 0);
            }
        }
    }

    void _tarjan(int u, int fa_edge) {
        _low[u] = _dfn[u] = ++_dfn_cnt;
        _stk_v.emplace_back(u);
        int son_cnt = 0;

        for (int e = g.head[u]; e; e = g.nxt[e]) {
            int v = g.to[e];
            if (!_dfn[v]) {
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
                    int x;
                    do {
                        x = _stk_v.back();
                        _stk_v.pop_back();
                        comp.emplace_back(x);
                    } while (x != v);
                    comp.emplace_back(u);
                    v_dcc.emplace_back(std::move(comp));
                }
            } else if ((e ^ 1) != fa_edge && _dfn[v] < _dfn[u]) {
                _low[u] = std::min(_low[u], _dfn[v]);
            }
        }
        if (u == _root && son_cnt == 0) {
            _stk_v.pop_back();
            ++v_dcc_cnt;
            v_dcc.emplace_back(std::vector<int>{u});
        }
    }

    void _build_tree() {
        bct_n = n + v_dcc_cnt;
        tree_adj.assign(bct_n + 1, {});
        rep(c, 1, v_dcc_cnt) {
            int sq = square_id(c);
            for (int u : v_dcc[c]) {
                tree_adj[u].emplace_back(sq);
                tree_adj[sq].emplace_back(u);
                tree_edges.emplace_back(u, sq);
            }
        }
    }
};
