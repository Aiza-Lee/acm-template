#include "aizalib.h"
/*
 * 2D Divide and Conquer DP Optimization (高维决策单调性分治优化)
 *
 * Overview:
 *     求解具有决策单调性的分层动态规划方程：
 *     dp[t][i] = min_{j < i}(dp[t - 1][j] + cost(t, j, i))。
 *     - 决策单调性与四边形不等式: 若代价函数满足四边形不等式，则最优决策点满足
 *       opt[t][i] <= opt[t][i + 1]。
 *     - 分层分治区间缩减: 在求解第 t 层的状态区间 [l, r] 时，先计算中点 mid = (l +
 *       r) / 2 在决策区间 [ql, qr] 内的最优决策 opt[t][mid]；
 *       由单调性将左侧决策区间收敛到 [ql, opt[t][mid]]，右侧收敛到 [opt[t][mid],
 *       qr]。
 *
 * API:
 *     DCDP2D(m, n, inf, cost)      — 构造 m 层、每层 n 个状态的优化器，
 *                                     传入转移代价函数
 *     set_base(base)               — 设置第 0 层基础状态数组 dp[0] (大小为 n + 1)
 *     solve_layer(t, l, r, ql, qr) — 分治求解第 t 层状态在 [l, r]
 *                                     区间的值与最优决策
 *     solve_all(l, r, ql, qr)      — 顺序分治求解全部 1..m 层
 *     value(t, i)                  — 查询第 t 层状态 i 的最优值 dp[t][i]
 *     decision(t, i)               — 查询第 t 层状态 i 的最优转移点 opt[t][i]
 *
 * Notes:
 *     1. Time: 单层 O(N log N)，全部 m 层总复杂度 O(M * N log N)。
 *     2. Space: O(M * N)。
 *     3. 状态与层数均采用 1-based 下标，决策点通常为 j in [0, i - 1]。
 */
template<class T, class F>
struct DCDP2D {
    int m, n;
    T inf;
    F cost;
    std::vector<std::vector<T>> dp;     // dp[t][i]: 第 t 层到状态 i 的最优值
    std::vector<std::vector<int>> opt;  // opt[t][i]: 第 t 层状态 i 的最优决策

    DCDP2D(int m, int n, T inf, F cost)
        : m(m), n(n), inf(inf), cost(cost), dp(m + 1, std::vector<T>(n + 1, inf)),
          opt(m + 1, std::vector<int>(n + 1, -1)) {}

    void set_base(const std::vector<T>& base) {
        AST((int)base.size() == n + 1);
        dp[0] = base;
    }

    void _solve(int t, int l, int r, int ql, int qr) {
        if (l > r) return;
        int mid = (l + r) >> 1;
        int rr = std::min(qr, mid - 1);
        AST(ql <= rr);

        T best_val = inf;
        int best_pos = -1;
        int split_pos = ql;
        rep(j, ql, rr) {
            if (dp[t - 1][j] == inf) continue;
            T cur = dp[t - 1][j] + cost(t, j, mid);
            if (cur < best_val) {
                best_val = cur;
                best_pos = j;
                split_pos = j;
            }
        }
        dp[t][mid] = best_val;
        opt[t][mid] = best_pos;
        _solve(t, l, mid - 1, ql, split_pos);
        _solve(t, mid + 1, r, split_pos, qr);
    }

    void solve_layer(int t, int l = 1, int r = -1, int ql = 0, int qr = -1) {
        AST(1 <= t && t <= m);
        if (r == -1) r = n;
        if (qr == -1) qr = n - 1;
        AST(1 <= l && l <= r && r <= n);
        AST(0 <= ql && ql <= qr && qr < n);
        AST(ql <= l - 1);
        std::fill(dp[t].begin(), dp[t].end(), inf);
        std::fill(opt[t].begin(), opt[t].end(), -1);
        _solve(t, l, r, ql, qr);
    }

    void solve_all(int l = 1, int r = -1, int ql = 0, int qr = -1) {
        if (r == -1) r = n;
        if (qr == -1) qr = n - 1;
        rep(t, 1, m) solve_layer(t, l, r, ql, qr);
    }

    T value(int t, int i) const {
        AST(1 <= t && t <= m);
        AST(1 <= i && i <= n);
        return dp[t][i];
    }

    int decision(int t, int i) const {
        AST(1 <= t && t <= m);
        AST(1 <= i && i <= n);
        return opt[t][i];
    }
};

template<class T, class F>
DCDP2D(int, int, T, F) -> DCDP2D<T, F>;
