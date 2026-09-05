#include "aizalib.h"
/*
 * Knuth优化 (Knuth Optimization)
 *
 * Overview:
 *      对于区间 DP：dp[i][j] = min_{i <= k < j}(dp[i][k] + dp[k + 1][j] + cost(i, j))。
 *      若 cost 满足四边形不等式与区间单调性，最优分割点满足 opt[i][j - 1] <= opt[i][j] <= opt[i + 1][j]，
 *      可将枚举范围降为均摊 O(1)，总体时间降至 O(n^2)。
 *
 * API:
 *      KnuthDP::solve(n, cost)      — 返回 dp[1][n]，cost(i, j) 为合并区间 [i, j] 的代价。
 *      KnuthDP::solve_full(n, cost) — 返回 pair(dp 表, opt 表)，opt[i][j] 记录最优分割点。
 *
 * Notes:
 *      1. Time: O(n^2)。
 *      2. Space: O(n^2)。
 *      3. 下标采用 1-based，dp[i][i] = 0。
 *      4. 适用条件：cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c) 且 cost(b, c) <= cost(a, d) (a <= b <= c <= d)。
 */

struct KnuthDP {
    /**
     * dp[i][j] = min_{i≤k<j} dp[i][k] + dp[k+1][j] + cost(i,j)
     * cost must satisfy quadrangle inequality + monotonicity
     */
    template<typename CostFn>
    static i64 solve(int n, CostFn cost) {
        if (n <= 0) return 0;
        std::vector<std::vector<i64>> dp(n + 2, std::vector<i64>(n + 2, 0));
        std::vector<std::vector<int>> opt(n + 2, std::vector<int>(n + 2));
        rep(i, 1, n) opt[i][i] = i;
        rep(len, 2, n) {
            rep(i, 1, n - len + 1) {
                int j = i + len - 1;
                dp[i][j] = std::numeric_limits<i64>::max();
                rep(k, opt[i][j - 1], std::min(j - 1, opt[i + 1][j])) {
                    i64 val = dp[i][k] + dp[k + 1][j] + cost(i, j);
                    if (val < dp[i][j]) {
                        dp[i][j] = val;
                        opt[i][j] = k;
                    }
                }
            }
        }
        return dp[1][n];
    }

    /**
     * Returns (dp table, opt table) for reconstruction.
     * dp[i][j]: optimal cost for interval [i,j].
     * opt[i][j]: optimal split point k (i ≤ k < j).
     */
    template<typename CostFn>
    static std::pair<std::vector<std::vector<i64>>, std::vector<std::vector<int>>>
    solve_full(int n, CostFn cost) {
        std::vector<std::vector<i64>> dp(n + 2, std::vector<i64>(n + 2, 0));
        std::vector<std::vector<int>> opt(n + 2, std::vector<int>(n + 2));
        if (n <= 0) return {std::move(dp), std::move(opt)};
        rep(i, 1, n) opt[i][i] = i;
        rep(len, 2, n) {
            rep(i, 1, n - len + 1) {
                int j = i + len - 1;
                dp[i][j] = std::numeric_limits<i64>::max();
                rep(k, opt[i][j - 1], std::min(j - 1, opt[i + 1][j])) {
                    i64 val = dp[i][k] + dp[k + 1][j] + cost(i, j);
                    if (val < dp[i][j]) {
                        dp[i][j] = val;
                        opt[i][j] = k;
                    }
                }
            }
        }
        return {std::move(dp), std::move(opt)};
    }
};
