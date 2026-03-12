#include "aizalib.h"
#include "0-base/Poly.hpp"

/**
 * 多项式平移 (Taylor Shift)
 * 算法介绍:
 * 		计算平移多项式 g(x) = f(x+c)。
 * 		利用差卷积求解: g[i] = sum_{j=i}^{n-1} f[j] * C(j, i) * c^{j-i}
 * 模板参数:
 * 
 * Interface:
 * 		Poly taylor_shift(const Poly& f, i64 c);
 * 
 * Note:
 * 		1. Time: O(n log n)
 * 		2. Space: O(n)
 */

namespace poly_ext {

Poly taylor_shift(const Poly& f, i64 c_in) {
	int n = f.size();
	if (n == 0) return f;
	
	int c = mod(c_in);
	std::vector<int> fac(n + 1, 1), ifac(n + 1);
	rep(i, 1, n) fac[i] = mul(fac[i - 1], i);
	ifac[n] = inv(fac[n]);
	per(i, n - 1, 0) ifac[i] = mul(ifac[i + 1], i + 1);
	
	Poly a(n), b(n);
	rep(i, 0, n - 1) {
		a[n - 1 - i] = mul(f[i], fac[i]); 
	}
	
	int pw = 1;
	rep(i, 0, n - 1) {
		b[i] = mul(pw, ifac[i]);
		pw = mul(pw, c);
	}
	
	Poly res = a * b;
	Poly g(n);
	rep(i, 0, n - 1) {
		g[i] = mul(res[n - 1 - i], ifac[i]);
	}
	
	return g;
}

} // namespace poly_ext
