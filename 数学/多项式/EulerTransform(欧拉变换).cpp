#include "0-base/Poly.hpp"

namespace poly_ext {

/**
 * 快速计算 \prod (1 - x^{a_i})
 * @param a 指数数组 {a_1, a_2, ...}
 * @param n 结果多项式的最高次项 (通常为 \sum a_i，如果只需要前 k 项可以传 k)
 * @return 结果多项式
 */
Poly euler_transform(const Poly& a, int n) {
	std::vector<int> sum(n + 1, 0);
	rep(i, 1, n) {
		if (i < (int)a.size() && a[i]) {
			int val = mul(a[i], i);
			for (int j = i; j <= n; j += i) inc(sum[j], val);
		}
	}

	Poly ln_b(n + 1);
	rep(i, 1, n) {
		ln_b[i] = mul(sum[i], inv(i));
	}

	return ln_b.exp(n + 1);
}

/**
 * 快速计算 \prod (1 - x^{a_i})
 * @param a 指数数组 {a_1, a_2, ...}
 * @param n 结果多项式的最高次项 (通常为 \sum a_i，如果只需要前 k 项可以传 k)
 * @return 结果多项式
 */
Poly fast_multiply_one_minus_xa(const std::vector<int>& a, int n) {
	Poly cnt(n + 1);
	for (int x : a) {
		if (x <= n) dec(cnt[x], 1);
	}
	return euler_transform(cnt, n);
}

} // namespace poly_ext
