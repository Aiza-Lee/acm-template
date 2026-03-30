#include "aizalib.h"

/**
 * Pollard-Rho (整数分解)
 * 算法介绍: 配合 Miller-Rabin 在 64 位范围内快速分解质因数。
 * 模板参数: None
 * Interface:
 * 		PollardRho::is_prime(n)      // 判定 n 是否为素数，Time: O(log^3 n)
 * 		PollardRho::get_one_factor(n) // 返回一个非平凡因子，Time: 期望 O(n^(1/4) log n)
 * 		PollardRho::factorize(n)     // 返回全部质因子并排序，Time: 期望 O(n^(1/4) log n)
 * 		PollardRho::factor_count(n)  // 返回 {质因子, 次数}，Time: 同 factorize
 * Note:
 * 		1. Time: 期望 O(n^(1/4) log n)
 * 		2. Space: O(log n)
 * 		3. 适用于正 64 位整数，n=1 时分解结果为空
 * 		4. 用法/技巧: 若只需判素数，直接调 is_prime 即可
 */
struct PollardRho {
	static inline std::mt19937_64 _rng{
		(u64)std::chrono::steady_clock::now().time_since_epoch().count()
	};

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

	static i64 _rand(i64 l, i64 r) {
		return std::uniform_int_distribution<i64>(l, r)(_rng);
	}

	static i64 _abs_diff(i64 x, i64 y) {
		return x >= y ? x - y : y - x;
	}

	static i64 _rho(i64 n) {
		if (!(n & 1)) return 2;
		if (n % 3 == 0) return 3;
		while (true) {
			i64 c = _rand(1, n - 1);
			i64 x = _rand(0, n - 1), y = x, d = 1;
			auto f = [&](i64 v) { return (_mul_mod(v, v, n) + c) % n; };
			while (d == 1) {
				x = f(x);
				y = f(f(y));
				d = std::gcd(_abs_diff(x, y), n);
			}
			if (d != n) return d;
		}
	}

	static void _factor(i64 n, std::vector<i64> &res) {
		if (n == 1) return;
		if (is_prime(n)) {
			res.emplace_back(n);
			return;
		}
		i64 d = _rho(n);
		_factor(d, res);
		_factor(n / d, res);
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

	static i64 get_one_factor(i64 n) {
		AST(n > 1 && !is_prime(n));
		return _rho(n);
	}

	static std::vector<i64> factorize(i64 n) {
		AST(n >= 1);
		std::vector<i64> res;
		_factor(n, res);
		std::sort(res.begin(), res.end());
		return res;
	}

	static std::vector<std::pair<i64, int>> factor_count(i64 n) {
		auto fac = factorize(n);
		std::vector<std::pair<i64, int>> res;
		for (i64 p : fac) {
			if (res.empty() || res.back().first != p) res.emplace_back(p, 1);
			else ++res.back().second;
		}
		return res;
	}
};
