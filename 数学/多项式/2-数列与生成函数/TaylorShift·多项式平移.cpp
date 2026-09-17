#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"

/*
 * Taylor Shift (多项式泰勒平移)
 *
 * Overview:
 *      计算多项式平移 g(x) = f(x + c)。
 *      根据泰勒展开 g_i = sum_{j=i}^{n-1} f_j * C(j, i) * c^{j-i}，整理得 i! * g_i
 *      = sum_{j=i}^{n-1} (j! * f_j) * (c^{j-i} / (j-i)!)。
 *      此式为标准下标差卷积，令 A_{n-1-j} = j! * f_j 与 B_k = c^k / k! 卷积，即可在
 *      O(n log n) 时间内求出 g(x)。
 *
 * API:
 *     taylor_shift(f, c) — 计算多项式 f(x + c)，返回与 f 同长的多项式。复杂度 O(n
 *                           log n) 时间，O(n) 空间。
 *
 * Notes:
 *      1. c 可为任意整数，内部自动对模数取模。
 *
 * Related:
 *      数学/多项式/1-卷积与变换/Convolution·卷积.cpp: 差卷积。
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
