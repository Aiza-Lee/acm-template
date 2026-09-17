#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"
/*
 * Cyclic Convolution (多项式循环卷积与循环移位)
 *
 * Overview:
 *      计算长度为 n 的模 (x^n - 1) 意义下的循环加法卷积与循环差卷积。
 *      利用普通 NTT 卷积后将所有高次项按模 n 累加折叠；循环移位对应在环 R[x] /
 *      (x^n - 1) 下乘单项式 x^k。
 *
 * API:
 *     cyclic_convolution(a, b, n)            — 计算循环加法卷积 res[k] =
 *                                               sum_{i+j=k (mod n)} a[i] * b[j]。
 *                                               复杂度 O((|a|+|b|) log(|a|+|b|))。
 *     cyclic_difference_convolution(a, b, n) — 计算循环差卷积 res[k] = sum_{i-j=k
 *                                               (mod n)} a[i] * b[j]。复杂度
 *                                               O((|a|+|b|) log(|a|+|b|))。
 *     cyclic_shift(a, k, n)                  — 计算多项式 a 模 (x^n - 1) 循环右移
 *                                               k 位。复杂度 O(|a|)。
 *
 * Notes:
 *      1. 要求 n > 0。
 *
 * Related:
 *      数学/多项式/1-卷积与变换/Convolution·卷积.cpp: 普通卷积与差卷积。
 */

namespace poly_ext {

// 模意义下的循环卷积 (长度为 n)
// res[k] = sum_{i+j = k (mod n)} a[i] * b[j]
Poly cyclic_convolution(const Poly& a, const Poly& b, int n) {
    auto c = a * b;
    Poly res(n);
    rep(i, 0, (int)c.size() - 1) {
        inc(res[i % n], c[i]);
    }
    return res;
}

// 模意义下的循环差卷积 (长度为 n)
// res[k] = sum_{i-j = k (mod n)} a[i] * b[j]
Poly cyclic_difference_convolution(const Poly& a, const Poly& b, int n) {
    Poly b_rev(n);
    rep(i, 0, (int)b.size() - 1) {
        int idx = (n - i % n) % n;
        inc(b_rev[idx], b[i]);
    }
    return cyclic_convolution(a, b_rev, n);
}

// 循环移位 (右移 k 位), 相当于乘 x^k mod (x^n - 1)
Poly cyclic_shift(const Poly& a, int k, int n) {
    Poly res(n);
    k = (k % n + n) % n;
    rep(i, 0, (int)a.size() - 1) {
        inc(res[(i + k) % n], a[i]);
    }
    return res;
}

} // namespace poly_ext