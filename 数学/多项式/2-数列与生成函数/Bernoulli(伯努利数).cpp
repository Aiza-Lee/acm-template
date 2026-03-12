#include "aizalib.h"
#include "0-base/Poly.hpp"

/**
 * 伯努利数 (Bernoulli Numbers)
 * 算法介绍:
 * 		由生成函数 x / (e^x - 1) 计算前 n 项伯努利数。
 * 		(e^x - 1) / x = sum (x^i / (i+1)!)
 * 		B = ( (e^x - 1) / x )^{-1}
 * 模板参数:
 * 
 * Interface:
 * 		Poly bernoulli(int n);
 * 
 * Note:
 * 		1. Time: O(n log n)
 * 		2. Space: O(n)
 * 		3. 计算出的为0-based结果多项式，第 i 项为 B_i 的模意义下值。
 */

namespace poly_ext {

Poly bernoulli(int n) {
	std::vector<int> fac(n + 2, 1), ifac(n + 2);
	rep(i, 1, n + 1) fac[i] = mul(fac[i - 1], i);
	ifac[n + 1] = inv(fac[n + 1]);
	per(i, n, 0) ifac[i] = mul(ifac[i + 1], i + 1);
	
	Poly a(n + 1);
	rep(i, 0, n) {
		a[i] = ifac[i + 1];
	}
	
	Poly b = a.inverse(n + 1);
	rep(i, 0, n) {
		b[i] = mul(b[i], fac[i]);
	}
	
	return b;
}

} // namespace poly_ext
