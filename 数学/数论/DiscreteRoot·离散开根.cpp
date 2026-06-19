#include "aizalib.h"

/**
 * Discrete Root (离散开根)
 * 算法介绍: 在素数模 p 下求解 x^k = a (mod p)，做法是原根化后转成线性同余。
 * 模板参数: None
 * Interface:
 * 		DiscreteRoot::solve(k, a, p) // 返回全部解并从小到大排序；要求 p 为素数且 k > 0
 * Note:
 * 		1. Time: 试除版约 O(求原根 + sqrt(p) + ans log p)
 * 		2. Space: O(sqrt(p))
 * 		3. 当 a != 0 时，先把 a 写成 g^t，再解 k * y = t (mod p - 1)
 * 		4. 用法/技巧: 当前版本仅处理素数模；一般合数模需要额外讨论
 */
struct DiscreteRoot {
	struct Hash {
		static u64 _splitmix64(u64 x) {
			x += 0x9e3779b97f4a7c15;
			x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
			x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
			return x ^ (x >> 31);
		}
		size_t operator()(u64 x) const {
			static const u64 seed =
				std::chrono::steady_clock::now().time_since_epoch().count();
			return _splitmix64(x + seed);
		}
	};

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

	static i64 _inv_mod(i64 a, i64 mod) {
		i64 x, y, g = _exgcd(a, mod, x, y);
		AST(g == 1);
		return _norm(x, mod);
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

	static i64 _primitive_root_prime(i64 p) {
		AST(p >= 2);
		if (p == 2) return 1;
		auto fac = _factor_distinct(p - 1);
		for (i64 g = 2; g < p; ++g) {
			bool ok = true;
			for (i64 q : fac) {
				if (_pow_mod(g, (p - 1) / q, p) == 1) {
					ok = false;
					break;
				}
			}
			if (ok) return g;
		}
		return -1;
	}

	static i64 _bsgs(i64 a, i64 b, i64 mod) {
		a = _norm(a, mod), b = _norm(b, mod);
		if (mod == 1) return 0;
		if (b == 1) return 0;
		i64 m = (i64)std::sqrt((ld)mod) + 1;
		AST(m <= std::numeric_limits<int>::max());
		int lim = (int)m;
		std::unordered_map<i64, int, Hash> mp;
		mp.reserve((size_t)lim * 2 + 1);
		i64 cur = b;
		rep(j, 0, lim - 1) {
			mp[cur] = j;
			cur = (i128)cur * a % mod;
		}
		i64 step = _pow_mod(a, m, mod);
		cur = step;
		rep(i, 1, lim) {
			auto it = mp.find(cur);
			if (it != mp.end()) return (i64)i * m - it->second;
			cur = (i128)cur * step % mod;
		}
		return -1;
	}

	static std::vector<i64> solve(i64 k, i64 a, i64 p) {
		AST(k > 0 && p >= 2);
		a = _norm(a, p);
		if (p == 2) return {a};
		if (!a) return {0};

		i64 g = _primitive_root_prime(p);
		i64 t = _bsgs(g, a, p);
		if (t == -1) return {};
		i64 mod = p - 1, d = std::gcd(k, mod);
		if (t % d) return {};

		i64 md = mod / d;
		i64 y0 = (i128)(t / d) * _inv_mod(k / d, md) % md;
		std::vector<i64> res;
		res.reserve((size_t)d);
		rep(i, 0, (int)d - 1) {
			i64 y = y0 + (i64)i * md;
			res.emplace_back(_pow_mod(g, y, p));
		}
		std::sort(res.begin(), res.end());
		res.erase(std::unique(res.begin(), res.end()), res.end());
		return res;
	}
};
