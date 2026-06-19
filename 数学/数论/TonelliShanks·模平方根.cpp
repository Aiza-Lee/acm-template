#include "aizalib.h"

/**
 * Tonelli-Shanks (模平方根)
 * 算法介绍: 在奇素数模 p 下求 x^2 = a (mod p) 的解。
 * 模板参数: None
 * Interface:
 * 		TonelliShanks::legendre(a, p) // 求勒让德符号；返回 0 / 1 / -1，Time: O(log p)
 * 		TonelliShanks::solve(a, p)    // 返回两个平方根 {x, y}，无解返回 {-1, -1}，Time: O(log^2 p)
 * Note:
 * 		1. Time: solve 为 O(log^2 p)
 * 		2. Space: O(1)
 * 		3. 仅适用于奇素数模；p = 2 时单独处理
 * 		4. 用法/技巧: 若得到一个解 x，另一个解就是 p - x；返回值会按从小到大排好
 */
struct TonelliShanks {
	static i64 _norm(i64 x, i64 mod) {
		x %= mod;
		return x < 0 ? x + mod : x;
	}

	static i64 _pow_mod(i64 a, i64 b, i64 mod) {
		i64 res = 1;
		for (a = _norm(a, mod); b; b >>= 1, a = (i128)a * a % mod)
			if (b & 1) res = (i128)res * a % mod;
		return res;
	}

	static int legendre(i64 a, i64 p) {
		a = _norm(a, p);
		if (!a) return 0;
		i64 x = _pow_mod(a, (p - 1) / 2, p);
		return x == 1 ? 1 : -1;
	}

	static std::pair<i64, i64> solve(i64 a, i64 p) {
		AST(p >= 2);
		a = _norm(a, p);
		if (p == 2) return {a, a};
		if (!a) return {0, 0};
		if (legendre(a, p) != 1) return {-1, -1};
		if (p % 4 == 3) {
			i64 x = _pow_mod(a, (p + 1) / 4, p);
			i64 y = p - x;
			if (x > y) std::swap(x, y);
			return {x, y};
		}

		i64 q = p - 1, s = 0;
		while (!(q & 1)) q >>= 1, ++s;
		i64 z = 2;
		while (legendre(z, p) != -1) ++z;
		i64 c = _pow_mod(z, q, p);
		i64 t = _pow_mod(a, q, p);
		i64 r = _pow_mod(a, (q + 1) / 2, p);
		i64 m = s;
		while (t != 1) {
			i64 i = 1, tt = (i128)t * t % p;
			while (tt != 1) tt = (i128)tt * tt % p, ++i;
			i64 b = _pow_mod(c, 1LL << (m - i - 1), p);
			r = (i128)r * b % p;
			c = (i128)b * b % p;
			t = (i128)t * c % p;
			m = i;
		}
		i64 x = r, y = p - r;
		if (x > y) std::swap(x, y);
		return {x, y};
	}
};
