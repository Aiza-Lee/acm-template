#pragma once
#include "aizalib.h"
#include "VDCC·点双连通分量.cpp"
/*
 * BlockCutTree·圆方树
 *
 * Overview:
 *     复用 VDCC 求解点双连通分量，并构建圆方树/森林（Block-Cut Tree）。
 *     原图节点作为圆点 (1..n)，各个点双连通分量作为方点 (n+1..n+v_dcc_cnt)。
 *     二分图结构：仅在方点 square_id(c) 与该分量包含的每个圆点 u 间连无向边。
 *
 * API:
 *     BlockCutTree(g)  — 构造函数，求解点双并建立圆方树，复杂度 O(V + E)
 *     n                — 原图节点数（圆点个数）
 *     v_dcc_cnt        — 点双连通分量总数（方点个数）
 *     bct_n            — 圆方树总节点数 (n + v_dcc_cnt)
 *     is_cut           — is_cut[u] 原图节点 u 是否为割点
 *     cuts             — 割点列表
 *     v_dcc            — 各点双包含的原图节点列表 (1-based)
 *     tree_adj         — 圆方树邻接表 (大小 bct_n + 1，圆方点无向连边)
 *     tree_edges       — 圆方树无向边列表，元素为 pair<int, int>(round_u, square_v)
 *     is_round(u)      — 判断节点 u 是否为圆点 (1 <= u <= n)
 *     is_square(u)     — 判断节点 u 是否为方点 (n < u <= bct_n)
 *     square_id(c)     — 获取点双分量 c 对应的方点编号 (n + c)
 *     dcc_of_square(s) — 获取方点 s 对应的点双分量编号 (s - n)
 *
 * Notes:
 *     1. 1-based indexing. 原图连通时圆方树为树，非连通时为森林。
 *     2. 割点判定等价性: 圆点 u 为割点 <=> 圆方树上 deg(u) > 1。
 *     3. 路径性质: 原图 u 到 v 所有简单路径的必须经过点，对应圆方树 u 到 v
 *        树上路径中的所有圆点；路径覆盖的所有点双对应路径中的所有方点。
 *     4. 路径点双赋权统计（铁人两项模型）: 圆点权 w = -1，方点权 w = |v_dcc[c]|，
 *        任意两圆点 s != t 间简单路径并集的必选中转点数即为圆方树上 (s, t)
 *        路径点权和。
 *
 * Related:
 *     VDCC·点双连通分量.cpp::VDCC — 底层点双连通分量求解器
 */

struct BlockCutTree {
    int n;                                       // 圆点个数
    int v_dcc_cnt;                               // 方点个数（点双分量数）
    int bct_n;                                   // 圆方树总点数
    std::vector<bool> is_cut;                    // 原图点是否为割点
    std::vector<int> cuts;                       // 割点列表
    std::vector<std::vector<int>> v_dcc;         // 各点双包含的原图点 (1-based)
    std::vector<std::vector<int>> tree_adj;      // 圆方树邻接表 (1..bct_n)
    std::vector<std::pair<int, int>> tree_edges; // 边列表 (round, square)

    BlockCutTree(const GraphFS& g) {
        VDCC core(g);
        n = g.n;
        v_dcc_cnt = core.v_dcc_cnt;
        bct_n = n + v_dcc_cnt;
        is_cut = std::move(core.is_cut);
        cuts = std::move(core.cuts);
        v_dcc = std::move(core.v_dcc);
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

    bool is_round(int u) const { return 1 <= u && u <= n; }
    bool is_square(int u) const { return n < u && u <= bct_n; }
    int square_id(int c) const { return n + c; }
    int dcc_of_square(int s) const { return s - n; }
};
