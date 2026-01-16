#include "0-base/Poly.hpp"
namespace poly_ext {

// 模意义下的循环卷积 (长度为 n)
// res[k] = sum_{i+j = k (mod n)} a[i] * b[j]
Poly cyclic_convolution(const Poly& a, const Poly& b, int n) {
	auto c = a * b;
	Poly res(n);
	rep(i, 0, (int)c.size() - 1) {
		inc(res[i % n], c[i]);
	}
	return res;
}

// 模意义下的循环差卷积 (长度为 n)
// res[k] = sum_{i-j = k (mod n)} a[i] * b[j]
Poly cyclic_difference_convolution(const Poly& a, const Poly& b, int n) {
	Poly b_rev(n);
	rep(i, 0, (int)b.size() - 1) {
		int idx = (n - i % n) % n;
		inc(b_rev[idx], b[i]);
	}
	return cyclic_convolution(a, b_rev, n);
}

// 循环移位 (右移 k 位), 相当于乘 x^k mod (x^n - 1)
Poly cyclic_shift(const Poly& a, int k, int n) {
	Poly res(n);
	k = (k % n + n) % n;
	rep(i, 0, (int)a.size() - 1) {
		inc(res[(i + k) % n], a[i]);
	}
	return res;
}

} // namespace poly_ext