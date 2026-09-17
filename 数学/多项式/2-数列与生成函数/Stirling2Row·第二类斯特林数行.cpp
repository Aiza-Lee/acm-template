#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"

/*
 * Stirling Numbers of the Second Kind Row (第二类斯特林数行)
 *
 * Overview:
 *      计算第二类斯特林数第 n 行的所有值 S(n, 0), S(n, 1), ..., S(n, n)。
 *      利用容斥卷积展开式 S(n, k) = sum_{i=0}^k ((-1)^{k-i} / (k-i)!) * (i^n / i!)，
 *      构造 A_i = (-1)^i / i! 与 B_i = i^n / i!，二者直接进行 NTT 多项式卷积即可在
 *      O(n log n) 时间内求出整行。
 *
 * API:
 *     stirling2_row(n) — 计算第二类斯特林数第 n 行，返回长为 n+1 的多项式。复杂度
 *                         O(n log n) 时间，O(n) 空间。
 *
 * Notes:
 *      1. 结果以 0-based 多项式给出，第 i 项为 S(n, i) mod 998244353。
 *      2. 特别地 S(0, 0) = 1。
 *
 * Related:
 *      数学/多项式/0-base/Poly·多项式全家桶.hpp: 底层 NTT 卷积。
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
