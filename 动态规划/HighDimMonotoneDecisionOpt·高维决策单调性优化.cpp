#include "aizalib.h"
/*
 * 高维决策单调性分治优化 (2D Divide and Conquer DP Optimization)
 *
 * Overview:
 *      多层 DP 形式：dp[t][i] = min_{j < i}(dp[t - 1][j] + cost(t, j, i))。
 *      若每层最优决策单调，即 opt[t][i] <= opt[t][i + 1]，可逐层分治优化。
 *
 * API:
 *      DCDP2D(m, n, inf, cost)      — 初始化多层优化器，cost(t, j, i) 返回从 j 到 i 的代价。
 *      set_base(base)               — 设置基础层 dp[0]。
 *      solve_layer(t, l, r, ql, qr) — 求解第 t 层。
 *      solve_all(l, r, ql, qr)      — 求解全部 1..m 层。
 *      value(t, i)                  — 获取 dp[t][i]。
 *      decision(t, i)               — 获取 opt[t][i]。
 *
 * Notes:
 *      1. Time: O(m * n log n)。
 *      2. Space: O(m * n)。
 *      3. 状态采用 1-based 下标，决策点通常为 j in [0, i - 1]。
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
