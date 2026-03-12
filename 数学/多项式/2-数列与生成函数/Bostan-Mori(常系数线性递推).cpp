#include "0-base/Poly.hpp"

namespace poly_ext {

/**
 * Bostan-Mori 算法求 [x^n] P(x)/Q(x)
 * 时间复杂度 O(d log d log n), d = deg(Q)
 */
int bostan_mori(i64 n, Poly P, Poly Q) {
	while (n > 0) {
		auto Q_neg = Q;
		for (int i = 1; i < (int)Q_neg.size(); i += 2) {
			Q_neg[i] = sub(0, Q_neg[i]);
		}

		auto U = P * Q_neg;
		auto V = Q * Q_neg;

		// P(x) <- U_odd(x) if n is odd, else U_even(x)
		// Note: U(x) = U_even(x^2) + x * U_odd(x^2)
		// If n is even, we want [x^{n/2}] U_even(x) / V_even(x)
		// If n is odd, we want [x^{(n-1)/2}] U_odd(x) / V_even(x)
		
		P.clear();
		for (int i = (n & 1); i < (int)U.size(); i += 2) {
			P.push_back(U[i]);
		}

		Q.clear();
		for (int i = 0; i < (int)V.size(); i += 2) {
			Q.push_back(V[i]);
		}

		n >>= 1;
	}
	if (Q.empty() || Q[0] == 0) return 0;
	return mul(P.empty() ? 0 : P[0], inv(Q[0]));
}

/**
 * 求线性递推数列第 n 项
 * a_n = \sum_{i=1}^d c_i a_{n-i}
 * @param n 要求第 n 项 (0-indexed)
 * @param c 递推系数 {c_1, c_2, ..., c_d}
 * @param a 初值 {a_0, a_1, ..., a_{d-1}}
 * @return a_n
 */
int linear_recurrence(long long n, const std::vector<int>& c, const std::vector<int>& a) {
	int d = c.size();
	if (n < d) return a[n];

	// Q(x) = 1 - \sum c_i x^i
	Poly Q(d + 1);
	Q[0] = 1;
	rep(i, 0, d - 1) {
		Q[i + 1] = sub(0, c[i]);
	}

	// P(x) = A(x) * Q(x) mod x^d
	Poly A(a);
	auto P = A * Q;
	P.resize(d);

	return bostan_mori(n, P, Q);
}

} // namespace poly_ext
