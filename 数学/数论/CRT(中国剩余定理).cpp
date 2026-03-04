#include "aizalib.h"

/**
 * 中国剩余定理 (CRT)
 * 算法介绍: 求解模数两两互质的线性同余方程组: x = r[i] (mod m[i])
 * Template Args: None
 * Interface: 
 *      i64 crt(const vector<i64>& r, const vector<i64>& m);
 * Note:
 *      1. Time: O(N log(max_m))
 *      2. Requirement: gcd(m[i], m[j]) = 1 for i \neq j
 *      3. Overflow protection using __int128
 */

namespace NumberTheory {

static i64 exgcd_crt(i64 a, i64 b, i64 &x, i64 &y) {
	if (!b) { x = 1; y = 0; return a; }
	i64 d = exgcd_crt(b, a % b, y, x);
	y -= a / b * x;
	return d;
}

// 求解 x \equiv r[i] \pmod{m[i]}
// 要求: m[i] 两两互质
i64 crt(const std::vector<i64>& r, const std::vector<i64>& m) {
	assert(r.size() == m.size());
	int n = r.size();
	i64 M = 1, ans = 0;
	
	for (i64 v : m) M *= v;
	
	rep(i, 0, n - 1) {
		i64 Mi = M / m[i], x, y;
		exgcd_crt(Mi, m[i], x, y);
		x = (x % m[i] + m[i]) % m[i]; // 保证逆元为正
		
		// term = r[i] * Mi * x (mod M)
		i64 term = (__int128)r[i] * Mi % M * x % M;
		ans = (ans + term) % M;
	}
	return (ans + M) % M;
}

} // namespace NumberTheory
