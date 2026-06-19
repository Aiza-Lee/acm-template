#include "aizalib.h"

/**
 * Extended CRT (扩展中国剩余定理)
 * 算法介绍: 逐个合并同余方程，处理模数不一定互质的情形。
 * 模板参数: None
 * Interface:
 * 		ExtendedCRT::merge(r0, m0, r1, m1) // 合并两个同余方程，Time: O(log max(m0, m1))
 * 		ExtendedCRT::solve(r, m)           // 合并整个方程组，Time: O(n log lcm)
 * Note:
 * 		1. Time: O(n log lcm)
 * 		2. Space: O(1)
 * 		3. 返回 {最小非负解, lcm}，无解返回 {-1, -1}
 * 		4. 用法/技巧: 输入余数可不先归一化，模板内部会转成 [0, m) 范围
 */
struct ExtendedCRT {
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

	static std::pair<i64, i64> merge(i64 r0, i64 m0, i64 r1, i64 m1) {
		AST(m0 > 0 && m1 > 0);
		r0 = _norm(r0, m0), r1 = _norm(r1, m1);
		i64 x, y, g = _exgcd(m0, m1, x, y), d = r1 - r0;
		if (d % g) return {-1, -1};

		i64 mod = m1 / g;
		i64 k = (i128)(d / g) * x % mod;
		i128 lcm = (i128)m0 / g * m1;
		AST(lcm <= std::numeric_limits<i64>::max());
		i64 nxt_m = (i64)lcm;
		i64 nxt_r = _norm(r0 + (i128)m0 * k % nxt_m, nxt_m);
		return {nxt_r, nxt_m};
	}

	static std::pair<i64, i64> solve(const std::vector<i64> &r, const std::vector<i64> &m) {
		AST(r.size() == m.size());
		int n = (int)r.size();
		if (!n) return {0, 1};
		AST(m[0] > 0);
		i64 r0 = _norm(r[0], m[0]), m0 = m[0];
		rep(i, 1, n - 1) {
			auto [nr, nm] = merge(r0, m0, r[i], m[i]);
			if (nm == -1) return {-1, -1};
			r0 = nr, m0 = nm;
		}
		return {r0, m0};
	}
};
