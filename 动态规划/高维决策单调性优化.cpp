#include "aizalib.h"

/**
 * 高维决策单调性分治优化
 * 算法介绍:
 * 		多层形式: `dp[t][i] = min_{j < i}(dp[t - 1][j] + cost(t, j, i))`。
 * 		若每层最优决策单调，即 `opt[t][i] <= opt[t][i + 1]`，可逐层分治优化。
 * 适用条件:
 * 		1. 每一层都需要最优决策单调。
 * 		2. 四边形不等式是常见充分条件，见文字资料。
 * 模板参数:
 * 		T: DP 值类型；F: 代价函数类型
 * Interface:
 * 		DCDP2D(m, n, inf, cost): 多层优化，`cost(t, j, i)` 返回第 t 层从 j 转移到 i 的附加代价
 * 		DCDP2D::set_base(base): 设置 `dp[0]`
 * 		DCDP2D::solve_layer(t, l = 1, r = n, ql = 0, qr = n - 1)
 * 		DCDP2D::solve_all(l = 1, r = n, ql = 0, qr = n - 1)
 * Note:
 * 		1. Time: O(m (r - l + 1) log (qr - ql + 2))
 * 		2. Space: O(mn)
 * 		3. 默认使用 1-based 状态下标，决策点通常为 `j in [0, i - 1]`
 * 		4. `cost(...)` 应能较快计算
 */
template<class T, class F>
struct DCDP2D {
	int m, n;
	T inf;
	F cost;
	std::vector<std::vector<T>> dp;		// dp[t][i]: 第 t 层到状态 i 的最优值
	std::vector<std::vector<int>> opt;	// opt[t][i]: 第 t 层状态 i 的最优决策

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
