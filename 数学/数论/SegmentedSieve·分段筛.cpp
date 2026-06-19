#include "aizalib.h"

/**
 * Segmented Sieve (分段筛 / 区间筛)
 * 算法介绍: 求区间 [L, R] 内的质数，先筛出 [2, sqrt(R)] 的质数，再标记区间内合数。
 * 模板参数: None
 * Interface:
 * 		SegmentedSieve(L, R) // 筛区间 [L, R] 内的质数，Time: O(sqrt(R) log log R + (R - L + 1) log log R)
 * 		sieve.is_prime[i]    // 偏移 i 对应的数 L + i 是否为质数
 * 		sieve.primes         // [L, R] 内所有质数
 * 		sieve.check(x)       // 判定区间内单点 x 是否为质数，Time: O(1)
 * Note:
 * 		1. Time: O(sqrt(R) log log R + (R - L + 1) log log R)
 * 		2. Space: O(sqrt(R) + (R - L + 1))
 * 		3. 适用于 R 很大但区间长度 R - L + 1 不大；偏移 i 对应值 L + i
 * 		4. 用法/技巧: 每个质数 p 从 max(p * p, ceil(L / p) * p) 开始标记；0 和 1 需单独去掉
 */
struct SegmentedSieve {
	i64 L, R;
	std::vector<int> base_primes;
	std::vector<char> is_prime;
	std::vector<i64> primes;

	static std::vector<int> _simple_sieve(int n) {
		std::vector<char> isp(n + 1, true);
		std::vector<int> res;
		if (n >= 0) isp[0] = false;
		if (n >= 1) isp[1] = false;
		rep(i, 2, n) {
			if (!isp[i]) continue;
			res.emplace_back(i);
			if ((i64)i * i > n) continue;
			for (int j = i * i; j <= n; j += i) isp[j] = false;
		}
		return res;
	}

	SegmentedSieve(i64 l, i64 r) : L(l), R(r), is_prime((size_t)(r - l + 1), true) {
		AST(0 <= l && l <= r);
		i64 lim64 = (i64)std::sqrt((ld)R);
		while ((lim64 + 1) * (lim64 + 1) <= R) ++lim64;
		while (lim64 * lim64 > R) --lim64;
		AST(lim64 <= std::numeric_limits<int>::max());
		base_primes = _simple_sieve((int)lim64);

		for (int p : base_primes) {
			i64 st = std::max((i64)p * p, (L + p - 1) / p * (i64)p);
			for (i64 x = st; x <= R; x += p) is_prime[(size_t)(x - L)] = false;
		}
		if (L == 0) is_prime[0] = false;
		if (L <= 1 && 1 <= R) is_prime[(size_t)(1 - L)] = false;

		for (size_t i = 0; i < is_prime.size(); ++i) {
			if (is_prime[i]) primes.emplace_back(L + (i64)i);
		}
	}

	bool check(i64 x) const {
		AST(L <= x && x <= R);
		return is_prime[(size_t)(x - L)];
	}
};
