#include "aizalib.h"
/*
 * Knuth Optimization (Knuth 决策单调性优化)
 *
 * Overview:
 *     区间动态规划的决策单调性双向夹逼优化算法：
 *     dp[i][j] = min_{i <= k < j}(dp[i][k] + dp[k + 1][j] + cost(i, j))。
 *     - 四边形不等式与区间单调性: 若代价函数满足四边形不等式 cost(a, c) + cost(b,
 *       d) <= cost(a, d) + cost(b, c) 且满足区间单调包含 cost(b, c) <= cost(a, d)
 *       (对任意 a <= b <= c <= d 成立)。
 *     - Knuth 夹逼定理: 最优分割点满足 opt[i][j - 1] <= opt[i][j] <= opt[i + 1][j]。
 *
 * API:
 *     KnuthDP::solve(n, cost)      — 求解区间 [1, n] 的最小合并代价 dp[1][n]，
 *                                     O(N^2)
 *     KnuthDP::solve_full(n, cost) — 求解并返回完整 dp 表与最优分割点 opt 表，
 *                                     O(N^2)
 *
 * Notes:
 *     1. Time: O(N^2)。
 *     2. Space: O(N^2)。
 *     3. 下标采用 1-based，基础状态 dp[i][i] = 0。
 *     4. 适用经典问题: 石子合并 (Garsia-Wachs 替代解)、最优二叉搜索树构建。
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
