#include "aizalib.h"
/*
 * Alpha-Beta Pruning (Alpha-Beta 剪枝)
 *
 * Overview:
 *     双人完全信息零和博弈下优化极大极小 (Minimax) 对抗搜索的经典剪枝算法：
 *     - Minimax 博弈树: 节点交替分为 MAX 层 (己方争取最大收益) 与 MIN 层
 *       (对手力图最小化己方收益)。
 *     - Alpha-Beta 剪枝准则: 维护祖先链上己方已知最大下界 alph 与对手允许的最小上界
 *       beta；当搜索导致 alph >= beta 时发生截断，舍弃后续无效分支。
 *     - 结构与工具: 配合静态局面估值函数，在博弈树展开中快速决出最优策略。
 *
 * API:
 *     AlphaBeta(n)                  — 构造包含 n 个节点的博弈树，O(N)
 *     init(num_nodes)               — 重新初始化包含 num_nodes 个节点的博弈树
 *     add_child(u, v)               — 添加博弈状态转移有向边 u -> v
 *     set_val(u, v)                 — 设置叶节点 u 的静态局面评估值
 *     search(u, alph, beta, is_max) — 从状态 u 启动 Alpha-Beta 剪枝搜索，
 *                                      返回博弈评估值
 *
 * Notes:
 *     1. Time: 最优节点排序下 O(b^(d/2))，随机顺序下约 O(b^(0.75d))，最劣 O(b^d)。
 *     2. Space: O(n)。
 *     3. 节点采用 1-based 下标。
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