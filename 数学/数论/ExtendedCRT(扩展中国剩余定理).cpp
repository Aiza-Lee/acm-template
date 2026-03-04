#include "aizalib.h"

/**
 * 扩展中国剩余定理 (Extended CRT)
 * 算法介绍: 求解模数不互质的线性同余方程组: x = r[i] (mod m[i])
 * Template Args: None
 * Interface: 
 *      pair<i64, i64> excrt(const vector<i64>& r, const vector<i64>& m);
 * Note:
 *      1. Time: O(N log(LCM))
 *      2. Return: {最小非负解 x, LCM}, 无解返回 {-1, -1}
 */

namespace NumberTheory {

i64 exgcd(i64 a, i64 b, i64 &x, i64 &y) {
	if (!b) { x = 1; y = 0; return a; }
	i64 d = exgcd(b, a % b, y, x);
	y -= a / b * x;
	return d;
}

// 求解 x \equiv r[i] \pmod{m[i]}
std::pair<i64, i64> excrt(const std::vector<i64>& r, const std::vector<i64>& m) {
	assert(r.size() == m.size());
	int n = r.size();
	i64 m0 = 1, r0 = 0;
	
	rep(i, 0, n - 1) {
		i64 m1 = m[i], r1 = r[i], k, y;
		i64 g = exgcd(m0, m1, k, y);
		// Sol: m0 * k + m1 * y = r1 - r0
		if ((r1 - r0) % g != 0) return {-1, -1};
		
		k = (__int128)((r1 - r0) / g) * k % (m1 / g);
		i64 lcm = m0 / g * m1;
		
		r0 = (r0 + (__int128)k * m0) % lcm;
		if (r0 < 0) r0 += lcm;
		m0 = lcm;
	}
	return {r0, m0};
}

} // namespace NumberTheory
