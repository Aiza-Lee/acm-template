#include "0-base/Poly.hpp"
namespace poly_ext {
// Lagrange Interpolation
// Given points (x_i, y_i), return P(x)
// O(n^2)
inline Poly lagrange_interpolation(const std::vector<int>& x, const std::vector<int>& y) {
	int n = x.size();
	Poly res(n);

	// Calculate M(x) = product(x - x_i)
	Poly M(1);
	M[0] = 1;
	rep(i, 0, n - 1) {
		Poly new_M(i + 2);
		rep(j, 0, i) {
			new_M[j] = sub(new_M[j], mul(M[j], x[i]));
			new_M[j + 1] = add(new_M[j + 1], M[j]);
		}
		M = std::move(new_M);
	}

	rep(j, 0, n - 1) {
		int den = 1;
		rep(i, 0, n - 1) {
			if (i == j) continue;
			den = mul(den, sub(x[j], x[i]));
		}
		int coef = mul(y[j], inv(den));

		Poly q(n);
		int carry = 0;
		per(i, n, 1) {
			q[i - 1] = add(M[i], mul(carry, x[j]));
			carry = q[i - 1];
		}

		rep(i, 0, n - 1) {
			res[i] = add(res[i], mul(q[i], coef));
		}
	}
	return res;
}

// Consecutive Lagrange Interpolation
// Given y_0, y_1, ..., y_n at x = 0, 1, ..., n
// Evaluate P(xi) in O(n)
inline int lagrange_consecutive(const std::vector<int>& y, int xi) {
	int n = y.size() - 1;
	if (xi <= n && xi >= 0) return y[xi];

	std::vector<int> pre(n + 1), suf(n + 1);
	pre[0] = xi;
	rep(i, 1, n) pre[i] = mul(pre[i - 1], sub(xi, i));
	suf[n] = sub(xi, n);
	per(i, n - 1, 0) suf[i] = mul(suf[i + 1], sub(xi, i));

	std::vector<int> fact(n + 1), inv_fact(n + 1);
	fact[0] = 1;
	rep(i, 1, n) fact[i] = mul(fact[i - 1], i);
	inv_fact[n] = inv(fact[n]);
	per(i, n - 1, 0) inv_fact[i] = mul(inv_fact[i + 1], i + 1);

	int ans = 0;
	rep(i, 0, n) {
		int num = mul((i == 0 ? 1 : pre[i - 1]), (i == n ? 1 : suf[i + 1]));
		int den = mul(inv_fact[i], inv_fact[n - i]);
		int term = mul(y[i], mul(num, den));
		if ((n - i) & 1) ans = sub(ans, term);
		else ans = add(ans, term);
	}
	return ans;
}
} // namespace poly_ext