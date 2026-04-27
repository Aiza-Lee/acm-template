#include "aizalib.h"
#include "0-base/Poly.hpp"

namespace poly_ext {

/**
 * Newton 恒等式
 * 用途:
 * 		在同一组根的幂和 p、基本对称 e、首一特征多项式 f 之间转换。
 * 		公式和应用见文字资料「数学/多项式/牛顿恒等式」。
 *
 * Interface:
 * 		Poly elementary_from_power(p, n)		// 幂和 p_1..p_n -> 基本对称 e_0..e_n
 * 		Poly power_from_elementary(e, m)		// 基本对称 e -> 幂和 p_0..p_m
 * 		Poly characteristic_from_power(p, n)	// 幂和 -> 特征多项式 prod(x-a_i)
 * 		Poly power_from_characteristic(f, m)	// 特征多项式 -> 幂和 p_0..p_m
 *
 * Note:
 * 		1. p[0] 表示 p_0，输入时可忽略；返回时 p[0] = 根数。
 * 		2. e[0] = 1，f 为低位在前且 f[n] = 1。
 * 		3. Time: O(n^2) 或 O(nm)，适合手写和中等规模。
 * 		4. elementary_from_power 需要 1..n 在模意义下可逆。
 */
struct NewtonIdentities {
	// 由幂和递推基本对称多项式
	static Poly elementary_from_power(const Poly& p, int n) {
		Poly e(n + 1);
		e[0] = 1;
		rep(k, 1, n) {
			int s = 0;
			rep(i, 1, k) {
				int v = (i < (int)p.size()) ? mul(e[k - i], p[i]) : 0;
				if (i & 1) inc(s, v);
				else dec(s, v);
			}
			e[k] = mul(s, inv(k));
		}
		return e;
	}

	// 由基本对称多项式递推幂和
	static Poly power_from_elementary(const Poly& e, int m) {
		int n = (int)e.size() - 1;
		Poly p(m + 1);
		p[0] = n;
		rep(k, 1, m) {
			int s = 0;
			rep(i, 1, std::min(k - 1, n)) {
				int v = mul(e[i], p[k - i]);
				if (i & 1) inc(s, v);
				else dec(s, v);
			}
			if (k <= n) {
				int v = mul(k, e[k]);
				if (k & 1) inc(s, v);
				else dec(s, v);
			}
			p[k] = s;
		}
		return p;
	}

	// 由幂和构造低位在前的首一特征多项式 prod(x-a_i)
	static Poly characteristic_from_power(const Poly& p, int n) {
		Poly e = elementary_from_power(p, n), f(n + 1);
		rep(k, 0, n) {
			f[n - k] = (k & 1) ? sub(0, e[k]) : e[k];
		}
		return f;
	}

	// 由低位在前的首一特征多项式恢复幂和
	static Poly power_from_characteristic(const Poly& f, int m) {
		int n = (int)f.size() - 1;
		AST(n >= 0 && f[n] == 1);
		Poly e(n + 1);
		rep(k, 0, n) {
			e[k] = (k & 1) ? sub(0, f[n - k]) : f[n - k];
		}
		return power_from_elementary(e, m);
	}
};

} // namespace poly_ext
