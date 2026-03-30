#include "aizalib.h"

/**
 * CRT (中国剩余定理)
 * 算法介绍: 求解模数两两互质的同余方程组 x = r[i] (mod m[i])。
 * 模板参数: None
 * Interface:
 * 		CRT::solve(r, m) // 解同余方程组并返回最小非负解，Time: O(n log max(m[i]))
 * Note:
 * 		1. Time: O(n log max(m[i]))
 * 		2. Space: O(1)
 * 		3. 返回 [0, M) 内最小解，其中 M = m[0] * m[1] * ... * m[n-1]
 * 		4. 用法/技巧: 若模数不互质请改用 ExtendedCRT
 */
struct CRT {
	static i64 _exgcd(i64 a, i64 b, i64 &x, i64 &y) {
		if (!b) {
			x = a >= 0 ? 1 : -1;
			y = 0;
			return std::abs(a);
		}
		i64 x1, y1, g = _exgcd(b, a % b, x1, y1);
		x = y1;
		y = x1 - a / b * y1;
		return g;
	}

	static i64 _norm(i64 x, i64 mod) {
		x %= mod;
		return x < 0 ? x + mod : x;
	}

	static i64 solve(const std::vector<i64> &r, const std::vector<i64> &m) {
		AST(r.size() == m.size());
		int n = (int)r.size();
		if (!n) return 0;

		i128 M = 1;
		for (i64 v : m) {
			AST(v > 0);
			M *= v;
			AST(M <= std::numeric_limits<i64>::max());
		}
		i64 mod = (i64)M, ans = 0;
		rep(i, 0, n - 1) {
			i64 ri = _norm(r[i], m[i]), Mi = mod / m[i], x, y;
			i64 g = _exgcd(Mi, m[i], x, y);
			AST(g == 1);
			x = _norm(x, m[i]);
			ans = (ans + (i128)ri * Mi % mod * x) % mod;
		}
		return _norm(ans, mod);
	}
};
