#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"

namespace poly_ext {

/*
 * Newton's Identities (牛顿恒等式与对称多项式变换)
 *
 * Overview:
 *      建立同一组变量根集合的幂和对称多项式 p_k = sum a_i^k、初等对称多项式 e_k =
 *      sum a_{i_1}...a_{i_k} 以及首一特征多项式 f(x) = prod (x - a_i)
 *      之间的快速双向转换。
 *      利用母函数对数导数关系 E(x) = sum e_k x^k = exp(sum (-1)^{k-1} p_k x^k / k)，
 *      结合多项式 Exp 与求逆在 O(n log n) 时间内完成三者互化。
 *
 * API:
 *     p_to_e(p, n) — 幂和 p_1..p_n 转化为初等对称多项式 e_0..e_n，复杂度 O(n log n)。
 *     e_to_p(e, m) — 初等对称多项式 e 转化为前 m 项幂和 p_0..p_m，复杂度 O(m log m)。
 *     p_to_f(p, n) — 幂和 p 转化为特征多项式 f(x) = prod(x - a_i)，复杂度 O(n log n)。
 *     f_to_p(f, m) — 特征多项式 f 转化为前 m 项幂和 p_0..p_m，复杂度 O(m log m)。
 *
 * Notes:
 *      1. 要求 1..n 在模意义下存在逆元。
 *      2. p[0] 返回根的总数 n，e[0] = 1，f 为升幂排列且最高次项系数 f[n] = 1。
 *
 * Related:
 *      数学/多项式/0-base/Poly·多项式全家桶.hpp: 多项式 Exp 与求逆。
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
