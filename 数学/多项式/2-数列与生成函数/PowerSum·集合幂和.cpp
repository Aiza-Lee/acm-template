#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"

/*
 * Set Power Sum (集合元素连续幂和)
 *
 * Overview:
 *      给定多重集合 A = {a_0, a_1, ..., a_{n-1}}，对所有 0 <= i < m，求出集合幂和
 *      g(i) = sum_{k=0}^{n-1} (a_k)^i。
 *      构造生成多项式 P(x) = prod (1 - a_k * x)，两边取对数展开得 ln P(x) = -
 *      sum_{i=1}^inf (g(i) / i) * x^i。
 *      利用分治 NTT 展开 P(x) 后求 Ln，即可在 O(n log^2 n + m log m)
 *      内求出全部幂和。
 *
 * API:
 *     power_sums(a, m) — 计算集合 a 的 0 到 m-1 次幂和，返回长为 m 的多项式。
 *                         复杂度 O(n log^2 n + m log m) 时间，O(n + m) 空间。
 *
 * Notes:
 *      1. 特别地，常数项 g(0) = n。
 *
 * Related:
 *      数学/多项式/2-数列与生成函数/NewtonIdentities·牛顿恒等式.cpp:
 *      对称多项式体系。
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
