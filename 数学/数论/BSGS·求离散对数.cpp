#include "aizalib.h"

/**
 * BSGS (Baby-Step Giant-Step)
 * 算法介绍: 求离散对数 a^x = b (mod p)，并支持 ExBSGS 处理 gcd(a, p) != 1。
 * 模板参数: None
 * Interface:
 * 		BSGS::solve(a, b, p)    // 在 gcd(a, p) = 1 时求最小非负解，Time: O(sqrt(p))
 * 		BSGS::ex_solve(a, b, p) // 一般情形求最小非负解，Time: O(sqrt(p))
 * Note:
 * 		1. Time: O(sqrt(p))
 * 		2. Space: O(sqrt(p))
 * 		3. 返回最小非负整数解，无解返回 -1
 * 		4. 用法/技巧: solve 仅适用于 gcd(a, p) = 1；一般情况直接用 ex_solve
 */
struct BSGS {
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

	static i64 _gcd(i64 a, i64 b) {
		return b ? _gcd(b, a % b) : std::abs(a);
	}

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
		i64 x, y;
		i64 g = _exgcd(a, mod, x, y);
		AST(g == 1);
		return _norm(x, mod);
	}

	static i64 solve(i64 a, i64 b, i64 p) {
		AST(p > 0);
		if (p == 1) return 0;
		a = _norm(a, p), b = _norm(b, p);
		if (b == 1) return 0;
		if (a == 0) return b == 0 ? 1 : -1;
		if (_gcd(a, p) != 1) return -1;
		i64 m = (i64)std::sqrt((ld)p) + 1;
		AST(m <= std::numeric_limits<int>::max());
		int lim = (int)m;

		std::unordered_map<i64, int, Hash> mp;
		mp.reserve((size_t)lim * 2 + 1);
		i64 cur = b;
		rep(j, 0, lim - 1) {
			mp[cur] = j;
			cur = (i128)cur * a % p;
		}

		i64 step = _pow_mod(a, m, p);
		cur = step;
		rep(i, 1, lim) {
			auto it = mp.find(cur);
			if (it != mp.end()) return (i64)i * m - it->second;
			cur = (i128)cur * step % p;
		}
		return -1;
	}

	static i64 ex_solve(i64 a, i64 b, i64 p) {
		AST(p > 0);
		if (p == 1) return 0;
		a = _norm(a, p), b = _norm(b, p);
		if (b == 1) return 0;
		if (a == 0) return b == 0 ? 1 : -1;

		i64 k = 1, d = 0;
		while (true) {
			i64 g = _gcd(a, p);
			if (g == 1) break;
			if (b % g) return -1;
			b /= g, p /= g;
			k = (i128)k * (a / g) % p;
			++d;
			if (k == b) return d;
		}
		if (p == 1) return d;
		i64 invk = _inv_mod(k, p);
		i64 t = solve(a, (i128)b * invk % p, p);
		return t == -1 ? -1 : t + d;
	}
};
