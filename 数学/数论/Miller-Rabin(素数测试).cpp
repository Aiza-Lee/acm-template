#include "aizalib.h"

/**
 * Miller-Rabin (素数测试)
 * 算法介绍: 用固定底数集合在 64 位范围内做确定性素性判定。
 * 模板参数: None
 * Interface:
 * 		MillerRabin::is_prime(n) // 判定 n 是否为素数，Time: O(log^3 n)
 * Note:
 * 		1. Time: O(log^3 n)
 * 		2. Space: O(1)
 * 		3. 适用于正 64 位整数判素
 * 		4. 用法/技巧: 若还要分解质因数，可直接改用 PollardRho
 */
struct MillerRabin {
	static i64 _mul_mod(i64 a, i64 b, i64 mod) {
		return (i128)a * b % mod;
	}

	static i64 _pow_mod(i64 a, i64 b, i64 mod) {
		i64 res = 1;
		for (a %= mod; b; b >>= 1, a = _mul_mod(a, a, mod))
			if (b & 1) res = _mul_mod(res, a, mod);
		return res;
	}

	static bool _check(i64 a, i64 s, i64 d, i64 n) {
		if (a % n == 0) return true;
		i64 x = _pow_mod(a, d, n);
		if (x == 1 || x == n - 1) return true;
		rep(i, 1, (int)s - 1) {
			x = _mul_mod(x, x, n);
			if (x == n - 1) return true;
		}
		return false;
	}

	static bool is_prime(i64 n) {
		if (n < 2) return false;
		for (i64 p : {2LL, 3LL, 5LL, 7LL, 11LL, 13LL, 17LL, 19LL, 23LL, 29LL, 31LL, 37LL}) {
			if (n % p == 0) return n == p;
		}
		i64 d = n - 1, s = 0;
		while (!(d & 1)) d >>= 1, ++s;
		for (i64 a : {2LL, 325LL, 9375LL, 28178LL, 450775LL, 9780504LL, 1795265022LL}) {
			if (!_check(a, s, d, n)) return false;
		}
		return true;
	}
};
