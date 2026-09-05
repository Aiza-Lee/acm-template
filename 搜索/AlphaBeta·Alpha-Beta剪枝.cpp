#include "aizalib.h"
/*
 * Alpha-Beta 剪枝 (Alpha-Beta Pruning)
 *
 * Overview:
 *      用于优化极大极小（Minimax）对抗树搜索的剪枝算法。
 *
 * API:
 *      AlphaBeta(n)                  — 构造包含 n 个节点的博弈树。
 *      add_child(u, v)               — 添加从节点 u 到子节点 v 的转移。
 *      set_val(u, v)                 — 设置叶节点 u 的评估值。
 *      search(u, alph, beta, is_max) — 从节点 u 出发执行 Alpha-Beta 搜索，返回评估值。
 *
 * Notes:
 *      1. Time: 最优 O(b^{d/2})，平均 O(b^{0.75d})，最劣 O(b^d)。
 *      2. Space: O(n)。
 *      3. 节点采用 1-based 下标。
 */

struct AlphaBeta {
    int n;
    std::vector<int> val;
    std::vector<std::vector<int>> son;

    AlphaBeta(int n = 0) : n(n), val(n + 1, 0), son(n + 1) {}

    void init(int num_nodes) {
        n = num_nodes;
        val.assign(n + 1, 0);
        son.assign(n + 1, {});
    }

    void add_child(int u, int v) {
        AST(1 <= u && u <= n && 1 <= v && v <= n);
        son[u].push_back(v);
    }

    void set_val(int u, int v) {
        AST(1 <= u && u <= n);
        val[u] = v;
    }

    int search(int u, int alph, int beta, bool is_max) {
        AST(1 <= u && u <= n);
        if (son[u].empty()) return val[u];
        if (is_max) {
            for (int v : son[u]) {
                alph = std::max(alph, search(v, alph, beta, !is_max));
                if (alph >= beta) break;
            }
            return alph;
        } else {
            for (int v : son[u]) {
                beta = std::min(beta, search(v, alph, beta, !is_max));
                if (alph >= beta) break;
            }
            return beta;
        }
    }
};