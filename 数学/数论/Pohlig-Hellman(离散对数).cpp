#include "aizalib.h"

/**
 * Pohlig-Hellman (离散对数 - 光滑阶)
 * 算法介绍: 当 ord(a) 的质因子较小（smooth）时，将 a^x = b (mod p) 的离散对数分解到各质数幂子群中分别求解，再通过 CRT 合并。
 * 模板参数: None
 * Interface:
 * 		PohligHellman::solve(a, b, p)                   // 求最小非负整数解，Time: O(sum e_i * q_i)
 * 		PohligHellman::solve(a, b, p, factors_of_pm1)    // 带 p-1 预分解因子，格式 vector<pair<质数, 指数>>
 * Note:
 * 		1. p 需为质数（群阶 = p-1），gcd(a, p) = 1；否则需 ExBSGS
 * 		2. Time: O(sum e_i * q_i)，ord(a) = prod q_i^{e_i}
 * 		   小规模子群使用暴力枚举 O(q)；当单个 q 较大（如 > 1e6）时可替换为 BSGS::solve 降至 O(sqrt(q))
 * 		3. Space: O(质因子个数)
 * 		4. 返回最小非负整数解，无解返回 -1
 * 		5. 用法/技巧: p 很大时需用 PollardRho 分解 p-1，内部自动将阶降至 ord(a)
 * 		6. 常见模数 p-1 分解:
 * 		   NTT 友好质数（smooth，Pohlig-Hellman 高效）:
 * 		     998244353  → 2^23 × 7 × 17
 * 		     469762049  → 2^26 × 7
 * 		     167772161  → 2^25 × 5
 * 		     1004535809 → 2^21 × 479
 * 		     104857601  → 2^22 × 5^2
 * 		   ACM 常用质数（非 smooth，BSGS 更合适）:
 * 		     1000000007 → p-1 = 2 × 500000003（大质因子）
 * 		     1000000009 → p-1 = 2^3 × 125000001（大质因子）
 */
struct PohligHellman {
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

	// 小规模子群 DLP: g 的阶为 q，求 g^x = h (mod p)，x in [0, q)
	// 当 q 较大（如 q > 1e6）时可替换为 BSGS::solve(g, h, p) 降至 O(sqrt(q))
	static i64 _dlp_small(i64 g, i64 h, i64 q, i64 p) {
		if (h == 1) return 0;
		i64 cur = 1;
		for (i64 x = 1; x < q; ++x) {
			cur = (i128)cur * g % p;
			if (cur == h) return static_cast<i64>(x);
		}
		return -1;
	}

	// 求 x mod q^e，逐位提取 base-q 数字
	static i64 _solve_prime_power(i64 a, i64 b, i64 n, i64 q, int e, i64 p) {
		i64 x_q = 0, q_pow = 1;
		i64 g = _pow_mod(a, n / q, p);  // ord(g) = q
		i64 inv_a = _inv_mod(a, p);

		for (int k = 0; k < e; ++k) {
			i64 base = (i128)_pow_mod(inv_a, x_q, p) * b % p;
			i64 h = _pow_mod(base, n / (q_pow * q), p);
			i64 d = _dlp_small(g, h, q, p);
			if (d == -1) return -1;
			x_q += d * q_pow;
			q_pow *= q;
		}
		return x_q;
	}

	static std::vector<std::pair<i64, int>> _factor(i64 n) {
		std::vector<std::pair<i64, int>> res;
		for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
			if (n % p) continue;
			int c = 0;
			while (n % p == 0) n /= p, ++c;
			res.emplace_back(p, c);
		}
		if (n > 1) res.emplace_back(n, 1);
		return res;
	}

	static i64 solve(i64 a, i64 b, i64 p) {
		AST(p > 1);
		a = _norm(a, p), b = _norm(b, p);
		if (b == 1) return 0;
		if (a == 0) return b == 0 ? 1 : -1;
		if (std::gcd(a, p) != 1) return -1;
		return solve(a, b, p, _factor(p - 1));
	}

	static i64 solve(i64 a, i64 b, i64 p, const std::vector<std::pair<i64, int>>& factors_of_pm1) {
		AST(p > 1);
		a = _norm(a, p), b = _norm(b, p);
		if (b == 1) return 0;
		if (a == 0) return b == 0 ? 1 : -1;
		if (std::gcd(a, p) != 1) return -1;

		// Compute ord(a) and its factorization from factors of p-1
		i64 n = p - 1;
		std::vector<std::pair<i64, int>> ord_factors;
		for (auto [q, e] : factors_of_pm1) {
			int cnt = e;
			while (cnt > 0 && _pow_mod(a, n / q, p) == 1) n /= q, --cnt;
			if (cnt > 0) ord_factors.emplace_back(q, cnt);
		}
		if (ord_factors.empty()) return -1;  // a ≡ 1 (mod p) and b ≠ 1, no solution

		i64 x = 0, M = 1;
		for (auto [q, e] : ord_factors) {
			i64 qe = 1;
			for (int k = 0; k < e; ++k) qe *= q;
			i64 x_i = _solve_prime_power(a, b, n, q, e, p);
			if (x_i == -1) return -1;

			// CRT: combine x (mod M) with x_i (mod qe)
			i64 diff = _norm(x_i - x, qe);
			i64 t, y;
			i64 g = _exgcd(M, qe, t, y);
			if (diff % g) return -1;
			t = (i128)_norm(t, qe) * (diff / g) % qe;
			x = x + M * t;
			M = M / g * qe;
			x = _norm(x, M);
		}
		return x;
	}
};
