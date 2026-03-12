#include "aizalib.h"
#include "0-base/Poly.hpp"

/**
 * 集合元素连续幂和 (Power Sum)
 * 算法介绍:
 * 给定集合 (或多重集) A = {a_0, a_1, ..., a_{n-1}}，
 * 对所有的 i 在 [0, m-1] 内，快速求出 g(i) = (a_0)^i + (a_1)^i + ... + (a_{n-1})^i
 * 
 * 推导过程:
 * 构造多项式: P(x) = (1 - a_0*x) * (1 - a_1*x) * ... * (1 - a_{n-1}*x)
 * 对 P(x) 取自然对数:
 * ln P(x) = sum( ln(1 - a_k*x) )
 * 由于 ln(1 - z) = - (z + z^2/2 + z^3/3 + ...)，代入得:
 * ln P(x) = - sum_k sum_i (a_k^i * x^i / i)
 *         = - sum_i (x^i / i * sum_k a_k^i)
 *         = - sum_i (g(i) / i * x^i)
 * 所以对于 i >= 1，有 g(i) = -i * (ln P(x) 第 i 次项的系数)。
 * 特别地，g(0) = n。
 * 
 * 模板参数:
 * 
 * Interface:
 * Poly power_sums(const std::vector<int>& a, int m);
 * 
 * Note:
 * 		1. Time: O(n log^2 n + m log m)
 * 		2. Space: O(n + m)
 * 		3. O(n log^2 n) 部分对应分治乘法，O(m log m) 对应多项式求导积分及求逆。
 */

namespace poly_ext {

Poly power_sums(const std::vector<int>& a, int m) {
	if (a.empty()) {
		return Poly(m, 0);
	}
	
	auto solve = [&](auto self, int l, int r) -> Poly {
		if (l == r) {
			// 构造 1 - a[l]*x
			Poly res(2);
			res[0] = 1;
			res[1] = sub(0, a[l]); 
			return res;
		}
		int mid = l + (r - l) / 2;
		return self(self, l, mid) * self(self, mid + 1, r);
	};
	
	Poly P = solve(solve, 0, (int)a.size() - 1);
	
	P.resize(m);
	
	Poly lnP = P.ln(m);
	
	Poly res(m);
	res[0] = a.size();
	rep(i, 1, m - 1) {
		res[i] = mul(lnP[i], sub(0, i)); 
	}
	
	return res;
}

} // namespace poly_ext
