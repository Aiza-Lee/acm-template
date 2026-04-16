#include "aizalib.h"

/**
 * 决策单调性分治优化
 * 算法介绍:
 * 		一层形式: `f[i] = min_{j < i} cost(j, i)`。
 * 		若最优决策单调，即 `opt[i] <= opt[i + 1]`，可用分治优化。
 * 适用条件:
 * 		1. 本模板需要最优决策单调。
 * 		2. 四边形不等式是常见充分条件，见文字资料。
 * 模板参数:
 * 		T: DP 值类型；F: 代价函数类型
 * Interface:
 * 		DCDP(n, inf, cost): 一层优化，`cost(j, i)` 返回从 j 转移到 i 的代价
 * 		DCDP::solve(l = 1, r = n, ql = 0, qr = n - 1)
 * Note:
 * 		1. Time: O((r - l + 1) log (qr - ql + 2))
 * 		2. Space: O(N)
 * 		3. 默认使用 1-based 状态下标，决策点通常为 `j in [0, i - 1]`
 * 		4. `cost(...)` 应能较快计算
 */
template<class T, class F>
struct DCDP {
	int n;
	T inf;
	F cost;
	std::vector<T> f;		// f[i]: 状态 i 的最优值
	std::vector<int> opt;	// opt[i]: 状态 i 的最优决策

	DCDP(int n, T inf, F cost)
		: n(n), inf(inf), cost(cost), f(n + 1, inf), opt(n + 1, -1) {}

	void _solve(int l, int r, int ql, int qr) {
		if (l > r) return;
		int mid = (l + r) >> 1;
		int rr = std::min(qr, mid - 1);
		AST(ql <= rr);

		T best_val = inf;
		int best_pos = ql;
		rep(j, ql, rr) {
			T cur = cost(j, mid);
			if (cur < best_val) {
				best_val = cur;
				best_pos = j;
			}
		}
		f[mid] = best_val;
		opt[mid] = best_pos;
		_solve(l, mid - 1, ql, best_pos);
		_solve(mid + 1, r, best_pos, qr);
	}

	void solve(int l = 1, int r = -1, int ql = 0, int qr = -1) {
		if (r == -1) r = n;
		if (qr == -1) qr = n - 1;
		AST(1 <= l && l <= r && r <= n);
		AST(0 <= ql && ql <= qr && qr < n);
		AST(ql <= l - 1);
		std::fill(f.begin(), f.end(), inf);
		std::fill(opt.begin(), opt.end(), -1);
		_solve(l, r, ql, qr);
	}

	T value(int i) const {
		AST(1 <= i && i <= n);
		return f[i];
	}

	int decision(int i) const {
		AST(1 <= i && i <= n);
		return opt[i];
	}
};

template<class T, class F>
DCDP(int, T, F) -> DCDP<T, F>;
