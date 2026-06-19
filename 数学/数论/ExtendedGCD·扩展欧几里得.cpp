#include "aizalib.h"

/**
 * Extended GCD (扩展欧几里得)
 * 算法介绍: 求 ax + by = gcd(a, b) 的一组特解，并可解模逆元与线性不定方程。
 * 模板参数: None
 * Interface:
 * 		ExtendedGCD::exgcd(a, b, x, y)                 // 求 ax + by = gcd(a, b) 的一组特解，Time: O(log min(|a|, |b|))
 * 		ExtendedGCD::inv_mod(a, mod, inv)             // 求 a 在 mod 下的逆元，不存在则返回 false，Time: O(log min(|a|, |mod|))
 * 		ExtendedGCD::solve_linear(a, b, c, x, y, dx, dy) // 解 ax + by = c 并给出通解步长，Time: O(log min(|a|, |b|))
 * Note:
 * 		1. Time: O(log min(|a|, |b|))
 * 		2. Space: O(log min(|a|, |b|))
 * 		3. 若 g = gcd(a, b)，则 ax + by = c 有解当且仅当 c % g == 0
 * 		4. 用法/技巧: solve_linear 返回一组特解，通解为 x = x0 + k * dx, y = y0 - k * dy
 */
struct ExtendedGCD {
	static i64 exgcd(i64 a, i64 b, i64 &x, i64 &y) {
		if (!b) {
			x = a >= 0 ? 1 : -1;
			y = 0;
			return std::abs(a);
		}
		i64 x1, y1, g = exgcd(b, a % b, x1, y1);
		x = y1;
		y = x1 - a / b * y1;
		return g;
	}

	static bool inv_mod(i64 a, i64 mod, i64 &inv) {
		AST(mod > 0);
		i64 x, y;
		if (exgcd(a, mod, x, y) != 1) return false;
		inv = (x % mod + mod) % mod;
		return true;
	}

	static bool solve_linear(i64 a, i64 b, i64 c, i64 &x, i64 &y, i64 &dx, i64 &dy) {
		if (!a && !b) {
			if (c) return false;
			x = y = dx = dy = 0;
			return true;
		}
		i64 g = exgcd(a, b, x, y);
		if (c % g) return false;
		i128 k = c / g;
		i128 xx = (i128)x * k, yy = (i128)y * k;
		AST(std::numeric_limits<i64>::min() <= xx && xx <= std::numeric_limits<i64>::max());
		AST(std::numeric_limits<i64>::min() <= yy && yy <= std::numeric_limits<i64>::max());
		x = (i64)xx, y = (i64)yy;
		dx = b / g, dy = a / g;
		return true;
	}
};
