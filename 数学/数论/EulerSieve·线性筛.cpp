#include "aizalib.h"

/**
 * Euler Sieve (线性筛)
 * 算法介绍: 在线性时间内求质数表、最小质因子、欧拉函数与莫比乌斯函数。
 * 模板参数: None
 * Interface:
 * 		EulerSieve(n)                    // 预处理 [1, n] 的 primes / minp / phi / mu，Time: O(n)
 * 		sieve.primes / minp / phi / mu  // 质数表、最小质因子、欧拉函数、莫比乌斯函数
 * 		sieve.is_prime(x)               // 判断 x 是否为质数，Time: O(1)
 * 		sieve.factorize(x)              // 分解 x 的质因数，Time: O(质因子个数)
 * Note:
 * 		1. Time: O(n)
 * 		2. Space: O(n)
 * 		3. 下标按 1..n，phi[1]=mu[1]=1，minp[x] 为最小质因子
 * 		4. 用法/技巧: 本文件下方保留了通用积性函数线性筛 MultiplicativeSieve<T>
 */
struct EulerSieve {
	int n;
	std::vector<int> primes;
	std::vector<int> minp;
	std::vector<int> phi;
	std::vector<int> mu;

	EulerSieve(int n) : n(n), minp(n + 1), phi(n + 1), mu(n + 1) {
		if (n >= 1) phi[1] = mu[1] = 1;
		rep(i, 2, n) {
			if (!minp[i]) {
				minp[i] = i;
				primes.emplace_back(i);
				phi[i] = i - 1;
				mu[i] = -1;
			}
			for (int p : primes) {
				if ((i64)i * p > n) break;
				minp[i * p] = p;
				if (i % p == 0) {
					phi[i * p] = phi[i] * p;
					mu[i * p] = 0;
					break;
				}
				phi[i * p] = phi[i] * (p - 1);
				mu[i * p] = -mu[i];
			}
		}
	}

	bool is_prime(int x) const {
		return 2 <= x && x <= n && minp[x] == x;
	}

	std::vector<std::pair<int, int>> factorize(int x) const {
		AST(1 <= x && x <= n);
		std::vector<std::pair<int, int>> res;
		while (x > 1) {
			int p = minp[x], c = 0;
			while (x % p == 0) x /= p, ++c;
			res.emplace_back(p, c);
		}
		return res;
	}
};

/**
 * Multiplicative Sieve (通用积性函数线性筛)
 * 算法介绍: 在线性筛过程中维护 i = rest * p^k，从而筛任意积性函数 f。
 * 模板参数: T
 * Interface:
 * 		MultiplicativeSieve<T>(n, calc_pk) // 预处理 [1, n] 的 f[i]，Time: O(n)
 * 		sieve.f[i]                         // f(i)
 * 		sieve.primes / minp / low / cnt   // 质数表、最小质因子、最小质因子最高幂、对应次数
 * 		sieve.is_prime(x)                  // 判断 x 是否为质数，Time: O(1)
 * 		sieve.factorize(x)                 // 分解 x 的质因数，Time: O(质因子个数)
 * 		calc_pk(p, c, pk) -> f(p^c)        // 返回 prime power 的函数值
 * Note:
 * 		1. Time: O(n)
 * 		2. Space: O(n)
 * 		3. 要求 f(1) = 1，且对 gcd(a, b) = 1 有 f(ab) = f(a) * f(b)
 * 		4. 用法/技巧: 只需定义 prime power 的值；例如 phi(p^c) = pk - pk / p，mu(p) = -1，mu(p^c) = 0 (c > 1)
 */
template<class T>
struct MultiplicativeSieve {
	int n;
	std::vector<int> primes;
	std::vector<int> minp;
	std::vector<int> low; // low[i] = i 中最小质因子的最高幂 p^k
	std::vector<int> cnt; // cnt[i] = i 中最小质因子的次数 k
	std::vector<T> f;

	template<class CalcPrimePower>
	MultiplicativeSieve(int n, CalcPrimePower calc_pk)
		: n(n), minp(n + 1), low(n + 1), cnt(n + 1), f(n + 1) {
		if (n >= 1) f[1] = T(1), low[1] = 1;
		rep(i, 2, n) {
			if (!minp[i]) {
				minp[i] = i;
				primes.emplace_back(i);
				low[i] = i;
				cnt[i] = 1;
				f[i] = calc_pk(i, 1, i);
			}
			for (int p : primes) {
				if ((i64)i * p > n) break;
				int x = i * p;
				minp[x] = p;
				if (p == minp[i]) {
					low[x] = low[i] * p;
					cnt[x] = cnt[i] + 1;
					f[x] = f[i / low[i]] * calc_pk(p, cnt[x], low[x]);
					break;
				}
				low[x] = p;
				cnt[x] = 1;
				f[x] = f[i] * calc_pk(p, 1, p);
			}
		}
	}

	bool is_prime(int x) const {
		return 2 <= x && x <= n && minp[x] == x;
	}

	std::vector<std::pair<int, int>> factorize(int x) const {
		AST(1 <= x && x <= n);
		std::vector<std::pair<int, int>> res;
		while (x > 1) {
			int p = minp[x], c = 0;
			while (x % p == 0) x /= p, ++c;
			res.emplace_back(p, c);
		}
		return res;
	}
};
