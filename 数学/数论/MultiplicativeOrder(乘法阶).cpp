#include "aizalib.h"

/**
 * Multiplicative Order (乘法阶)
 * 算法介绍: 在 gcd(a, mod) = 1 时求 a 在模 mod 意义下的最小正阶。
 * 模板参数: None
 * Interface:
 * 		MultiplicativeOrder::order(a, mod) // 求最小正整数 k 使 a^k = 1 (mod mod)；若不存在返回 -1
 * Note:
 * 		1. Time: 试除版约 O(sqrt(mod) + sqrt(phi(mod)) log mod)
 * 		2. Space: O(质因子个数)
 * 		3. 阶一定整除 phi(mod)
 * 		4. 用法/技巧: 若 mod 很大，可把内部试除分解替换成 PollardRho
 */
struct MultiplicativeOrder {
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

	static std::vector<i64> _factor_distinct(i64 n) {
		std::vector<i64> res;
		for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
			if (n % p) continue;
			res.emplace_back(p);
			while (n % p == 0) n /= p;
		}
		if (n > 1) res.emplace_back(n);
		return res;
	}

	static i64 _phi(i64 n) {
		i64 res = n;
		for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
			if (n % p) continue;
			res = res / p * (p - 1);
			while (n % p == 0) n /= p;
		}
		if (n > 1) res = res / n * (n - 1);
		return res;
	}

	static i64 order(i64 a, i64 mod) {
		if (mod <= 1 || std::gcd(a, mod) != 1) return -1;
		i64 ord = _phi(mod);
		for (i64 p : _factor_distinct(ord)) {
			while (ord % p == 0 && _pow_mod(a, ord / p, mod) == 1) ord /= p;
		}
		return ord;
	}
};
