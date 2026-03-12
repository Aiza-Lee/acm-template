#include "aizalib.h"
#include "0-base/Poly.hpp"

/**
 * 第二类斯特林数行 (Stirling Numbers of the Second Kind)
 * 算法介绍:
 * 		计算第二类斯特林数第 n 行，即 S(n, x) 对于 x = 0...n。
 * 		容斥原理: S(n, k) = sum_{i=0}^{k} ( (-1)^{k-i} / (k-i)! ) * ( i^n / i! )
 * 模板参数:
 * 
 * Interface:
 * 		Poly stirling2_row(int n);
 * 
 * Note:
 * 		1. Time: O(n log n)
 * 		2. Space: O(n)
 * 		3. 0-based 结果多项式，第 i 项即为 S(n, i)。
 */

namespace poly_ext {

Poly stirling2_row(int n) {
	if (n == 0) {
		Poly f(1);
		f[0] = 1;
		return f;
	}
	
	std::vector<int> ifac(n + 1);
	int fac_n = 1;
	rep(i, 1, n) fac_n = mul(fac_n, i);
	ifac[n] = inv(fac_n);
	per(i, n - 1, 0) ifac[i] = mul(ifac[i + 1], i + 1);
	
	Poly a(n + 1), b(n + 1);
	rep(i, 0, n) {
		a[i] = (i % 2 == 1) ? sub(0, ifac[i]) : ifac[i];
		b[i] = mul(fp(i, n), ifac[i]);
	}
	
	Poly ans = a * b;
	ans.resize(n + 1);
	return ans;
}

} // namespace poly_ext
