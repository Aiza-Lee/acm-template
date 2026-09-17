#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"

/*
 * Bernoulli Numbers (伯努利数)
 *
 * Overview:
 *      基于指数生成函数 B(x) = x / (e^x - 1) = sum_{i=0}^inf B_i * x^i / i! 通过对
 *      A(x) = (e^x - 1) / x = sum_{i=0}^inf x^i / (i+1)! 进行多项式求逆，在 O(n log
 *      n) 时间内计算出前 n+1 项伯努利数 B_0, B_1, ..., B_n。
 *
 * API:
 *     bernoulli(n) — 计算 B_0 到 B_n 的值，返回长为 n+1 的多项式。复杂度 O(n log
 *                     n) 时间，O(n) 空间。
 *
 * Notes:
 *      1. 结果以 0-based 给出，b[i] 即为 B_i mod 998244353。
 *      2. 约定 B_1 = -1/2 (取模为 (md-1)/2)。
 *
 * Related:
 *      数学/多项式/0-base/Poly·多项式全家桶.hpp: 多项式求逆支持。
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
