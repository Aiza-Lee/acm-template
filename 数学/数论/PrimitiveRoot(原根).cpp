#include "aizalib.h"

/**
 * Primitive Root (原根)
 * 算法介绍: 判定模数是否存在原根，并求一个最小原根。
 * 模板参数: None
 * Interface:
 * 		PrimitiveRoot::has_primitive_root(mod) // 判断 mod 是否存在原根，Time: 试除版约 O(sqrt(mod))
 * 		PrimitiveRoot::is_primitive_root(g, mod) // 判断 g 是否为 mod 的原根，Time: 试除版约 O(sqrt(mod))
 * 		PrimitiveRoot::find(mod)               // 求最小原根，无原根返回 -1；时间取决于试除分解 + 枚举
 * Note:
 * 		1. Time: 试除分解版本约为 O(sqrt(mod))
 * 		2. Space: O(因子个数(phi(mod)))
 * 		3. 无原根时 find 返回 -1；约定 mod=1 时返回 0
 * 		4. 用法/技巧: 若模数很大，可把内部试除分解替换成 PollardRho
 */
struct PrimitiveRoot {
	static i64 _pow_mod(i64 a, i64 b, i64 mod) {
		i64 res = 1;
		for (a %= mod; b; b >>= 1, a = (i128)a * a % mod)
			if (b & 1) res = (i128)res * a % mod;
		return res;
	}

	static std::vector<std::pair<i64, int>> _factor_count(i64 n) {
		std::vector<std::pair<i64, int>> res;
		if (n <= 1) return res;
		for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
			if (n % p) continue;
			int c = 0;
			while (n % p == 0) n /= p, ++c;
			res.emplace_back(p, c);
		}
		if (n > 1) res.emplace_back(n, 1);
		return res;
	}

	static std::vector<i64> _factor_distinct(i64 n) {
		std::vector<i64> res;
		for (auto [p, _] : _factor_count(n)) res.emplace_back(p);
		return res;
	}

	static i64 _phi(i64 n) {
		i64 res = n;
		for (auto [p, _] : _factor_count(n)) res = res / p * (p - 1);
		return res;
	}

	static bool has_primitive_root(i64 mod) {
		AST(mod >= 1);
		if (mod == 1 || mod == 2 || mod == 4) return true;
		auto fac = _factor_count(mod);
		if (fac.size() == 1) return fac[0].first != 2;
		return fac.size() == 2 && fac[0].first == 2 && fac[0].second == 1
			&& fac[1].first != 2;
	}

	static bool is_primitive_root(i64 g, i64 mod) {
		if (!has_primitive_root(mod) || std::gcd(g, mod) != 1) return false;
		if (mod == 1) return g == 0;
		i64 phi = _phi(mod);
		for (i64 p : _factor_distinct(phi)) {
			if (_pow_mod(g, phi / p, mod) == 1) return false;
		}
		return true;
	}

	static i64 find(i64 mod) {
		if (!has_primitive_root(mod)) return -1;
		if (mod == 1) return 0;
		i64 phi = _phi(mod);
		auto fac = _factor_distinct(phi);
		for (i64 g = 1; g < mod; ++g) {
			if (std::gcd(g, mod) != 1) continue;
			bool ok = true;
			for (i64 p : fac) {
				if (_pow_mod(g, phi / p, mod) == 1) {
					ok = false;
					break;
				}
			}
			if (ok) return g;
		}
		return -1;
	}
};
