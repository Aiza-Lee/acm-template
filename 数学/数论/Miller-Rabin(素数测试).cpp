#include "aizalib.h"

/**
 * Miller-Rabin 素数测试
 * 算法介绍: 概率性素数测试算法 (在选定特定底数集时变为确定性算法)
 * Interface: 
 *      bool is_prime(i64 n);
 * Note:
 *      1. Time: O(k log^3 n)
 *      2. Deterministic for n < 2^64 using standard bases
 *      3. Uses __int128 for modular arithmetic to prevent overflow
 */

namespace NumberTheory {

// 快速幂 (With Mod)
static i64 qpow(i64 a, i64 b, i64 m) {
	i64 res = 1;
	a %= m;
	while (b) {
		if (b & 1) res = (__int128)res * a % m;
		a = (__int128)a * a % m;
		b >>= 1;
	}
	return res;
}

// 单次测试: true 可能为素数, false 必定合数
static bool miller_rabin(i64 n, i64 a, i64 d, int s) {
	if (n <= a) return true;
	i64 x = qpow(a, d, n);
	if (x == 1 || x == n - 1) return true;
	while (--s) {
		x = (__int128)x * x % n;
		if (x == n - 1) return true;
	}
	return false;
}

// 确定性底数集 (n < 4,759,123,141)
const std::vector<i64> small_bases = { 2, 7, 61 };

// 确定性底数集 (n < 2^64)
const std::vector<i64> large_bases = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };

bool is_prime(i64 n) {
	if (n < 4) return n == 2 || n == 3;
	if (n % 2 == 0 || n % 3 == 0) return false;

	i64 d = n - 1; 
	int s = 0;
	while (!(d & 1)) d >>= 1, s++;

	// 选择合适的底数集
	const auto& bases = (n < 4759123141LL) ? small_bases : large_bases;

	for (i64 a : bases) {
		if (!miller_rabin(n, a, d, s)) return false;
	}
	return true;
}

} // namespace NumberTheory
