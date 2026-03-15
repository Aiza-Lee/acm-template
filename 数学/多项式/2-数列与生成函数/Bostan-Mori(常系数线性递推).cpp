#include "0-base/Poly.hpp"

namespace poly_ext {

/**
 * Bostan-Mori(常系数线性递推)
 * 算法介绍: 用于求解有理分式 [x^n] P(x)/Q(x) 和常系数线性递推数列第 n 项。数学原理见文字部分。
 * 模板参数: 无
 * Interface: 
 *  - bostan_mori(n, P, Q): 求 [x^n] P(x)/Q(x)
 *  - linear_recurrence(n, c, a): 求常系数线性递推数列第 n 项
 * Note:
 * 		1. Time: O(d log d log n), d = deg(Q)
 * 		2. Space: O(d)
 */

/**
 * @brief Bostan-Mori 算法求 [x^n] P(x)/Q(x)
 * @param n 求第 n 项的系数 (0-indexed)
 * @param P 分子多项式
 * @param Q 分母多项式 (需要保证 Q[0] != 0)
 * @return [x^n] P(x)/Q(x) 取模后的值
 */
int bostan_mori(i64 n, Poly P, Poly Q) {
	while (n > 0) {
		auto Q_neg = Q;
		for (int i = 1; i < (int)Q_neg.size(); i += 2)
			Q_neg[i] = sub(0, Q_neg[i]);

		auto U = P * Q_neg, V = Q * Q_neg;
		
		P.clear();
		for (int i = (n & 1); i < (int)U.size(); i += 2)
			P.push_back(U[i]);

		Q.clear();
		for (int i = 0; i < (int)V.size(); i += 2)
			Q.push_back(V[i]);

		n >>= 1;
	}
	if (Q.empty() || Q[0] == 0) return 0;
	return mul(P.empty() ? 0 : P[0], inv(Q[0]));
}

/**
 * @brief 求线性递推数列第 n 项
 * 满足 a_n = \sum_{i=1}^d c_i a_{n-i}
 * @param n 要求第 n 项 (0-indexed)
 * @param c 递推系数 {c_1, c_2, ..., c_d}
 * @param a 初值 {a_0, a_1, ..., a_{d-1}}
 * @return a_n 取模后的值
 */
int linear_recurrence(long long n, const std::vector<int>& c, const std::vector<int>& a) {
	int d = c.size();
	if (n < d) return a[n];

	Poly Q(d + 1);
	Q[0] = 1;
	rep(i, 0, d - 1)
		Q[i + 1] = sub(0, c[i]);

	Poly A(a);
	auto P = A * Q;
	P.resize(d);

	return bostan_mori(n, P, Q);
}

} // namespace poly_ext
