#include "aizalib.h"

/**
 * Floor Sum (整除分块)
 * 算法介绍: 利用 ⌊n/i⌋ 取值只有 O(√n) 种的性质，将求和按值分块计算。
 * Interface:
 *      FloorSum::sum(n)        — ∑_{i=1}^{n} ⌊n/i⌋
 *      FloorSum::sum(n, k)     — ∑_{i=1}^{n} ⌊k/i⌋
 *      FloorSum::sum_w(n, k, pref) — ∑_{i=1}^{n} f(i)⋅⌊k/i⌋ (pref 为 f 前缀和)
 *      FloorSum::sum2(n, m)    — ∑_{i=1}^{min(n,m)} ⌊n/i⌋⋅⌊m/i⌋
 * Note:
 *      1. Time: O(√n) per call
 *      2. Space: O(1)
 *      3. 核心思想：对于 i ∈ [l, r], ⌊n/l⌋ = ⌊n/r⌋, 其中 r = n / (n / l)
 *      4. 原理: floor_sum_w 需要 f 的前缀和数组 (1-indexed), pref[0] = 0
 */
struct FloorSum {
	static i64 sum(i64 n) {
		i64 ans = 0;
		for (i64 l = 1, r; l <= n; l = r + 1) {
			r = n / (n / l);
			ans += (n / l) * (r - l + 1);
		}
		return ans;
	}

	static i64 sum(i64 n, i64 k) {
		i64 ans = 0;
		for (i64 l = 1, r; l <= n; l = r + 1) {
			if (k / l == 0) break;
			r = std::min(k / (k / l), n);
			ans += (k / l) * (r - l + 1);
		}
		return ans;
	}

	static i64 sum_w(i64 n, i64 k, const std::vector<i64>& pref) {
		i64 ans = 0;
		for (i64 l = 1, r; l <= n; l = r + 1) {
			if (k / l == 0) break;
			r = std::min(k / (k / l), n);
			ans += (pref[r] - pref[l - 1]) * (k / l);
		}
		return ans;
	}

	static i64 sum2(i64 n, i64 m) {
		i64 ans = 0;
		i64 lim = std::min(n, m);
		for (i64 l = 1, r; l <= lim; l = r + 1) {
			r = std::min(n / (n / l), m / (m / l));
			ans += (n / l) * (m / l) * (r - l + 1);
		}
		return ans;
	}
};
