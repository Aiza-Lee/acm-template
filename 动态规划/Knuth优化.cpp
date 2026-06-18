#include "aizalib.h"

/**
 * Knuth优化 (Knuth Optimization)
 * 算法介绍:
 *     对于区间DP形式: dp[i][j] = min_{i≤k<j} dp[i][k] + dp[k+1][j] + cost(i,j)
 *     若cost满足四边形不等式与区间包含单调性, 最优决策opt[i][j]满足:
 *     opt[i][j-1] ≤ opt[i][j] ≤ opt[i+1][j]
 *     将内层k的枚举范围从O(n)降至均摊O(1), 总体O(n^2)
 * 适用条件:
 *     1. cost(a,c)+cost(b,d) ≤ cost(a,d)+cost(b,c),  a≤b≤c≤d (四边形不等式)
 *     2. cost(b,c) ≤ cost(a,d),  a≤b≤c≤d (区间包含单调性)
 * 模板参数: None
 * Interface:
 *     solve(n, cost):  返回dp[1][n], cost(i,j)为合并区间[i,j]的代价
 *     solve_full(n, cost):  返回(dp表, opt表), opt[i][j]存最优分割点
 * Note:
 *     1. Time: O(n^2)
 *     2. Space: O(n^2)
 *     3. 下标1-based, dp[i][i]=0
 *     4. 常见场景: 石子合并
 *     5. 使用前请自行验证cost是否满足两个不等式条件
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
