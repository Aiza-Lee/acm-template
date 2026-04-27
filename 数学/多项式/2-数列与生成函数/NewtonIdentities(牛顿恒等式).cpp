#include "aizalib.h"
#include "0-base/Poly.hpp"

namespace poly_ext {

/**
 * Newton 恒等式
 * 用途:
 * 		在同一组根的 p、e、f 之间转换。
 * 		公式和应用见文字资料「数学/多项式/牛顿恒等式」。
 *
 * Interface:
 * 		Poly p_to_e(p, n)	// p_1..p_n -> e_0..e_n
 * 		Poly e_to_p(e, m)	// e -> p_0..p_m
 * 		Poly p_to_f(p, n)	// p -> f(x)=prod(x-a_i)
 * 		Poly f_to_p(f, m)	// f -> p_0..p_m
 *
 * Note:
 * 		1. p[0] 表示 p_0，输入时可忽略；返回时 p[0] = 根数。
 * 		2. e[0] = 1，f 为低位在前且 f[n] = 1。
 * 		3. Time: O(n log n) / O(m log m)。
 * 		4. 由 p 恢复 e 时需要 1..n 在模意义下可逆。
 */
struct NewtonIdentities {
	// E(x)=sum e_k x^k = exp(sum (-1)^{k-1} p_k x^k / k)
	static Poly p_to_e(const Poly& p, int n) {
		Poly ln_e(n + 1);
		rep(k, 1, n) if (k < (int)p.size()) {
			ln_e[k] = mul(p[k], inv(k));
			if (!(k & 1)) ln_e[k] = sub(0, ln_e[k]);
		}
		return ln_e.exp(n + 1);
	}

	// E'(x)/E(x)=sum (-1)^{k-1} p_k x^{k-1}
	static Poly e_to_p(const Poly& e, int m) {
		AST(!e.empty() && e[0] == 1);
		int n = (int)e.size() - 1;
		Poly p(m + 1);
		p[0] = n;
		if (m == 0) return p;
		Poly de = e.deriv();
		de.resize(m);
		Poly q = de * e.inverse(m);
		q.resize(m);
		rep(k, 1, m) p[k] = (k & 1) ? q[k - 1] : sub(0, q[k - 1]);
		return p;
	}

	// 由幂和构造低位在前的首一特征多项式 prod(x-a_i)
	static Poly p_to_f(const Poly& p, int n) {
		Poly e = p_to_e(p, n), f(n + 1);
		rep(k, 0, n) {
			f[n - k] = (k & 1) ? sub(0, e[k]) : e[k];
		}
		return f;
	}

	// 由低位在前的首一特征多项式恢复幂和
	static Poly f_to_p(const Poly& f, int m) {
		int n = (int)f.size() - 1;
		AST(n >= 0 && f[n] == 1);
		Poly e(n + 1);
		rep(k, 0, n) {
			e[k] = (k & 1) ? sub(0, f[n - k]) : f[n - k];
		}
		return e_to_p(e, m);
	}
};

} // namespace poly_ext
