#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"
/*
 * Polynomial Convolutions (多项式普通卷积与差卷积)
 *
 * Overview:
 *      提供普通加法卷积与差卷积（下标差卷积）。
 *      加法卷积 res[k] = sum_{i+j=k} a[i] * b[j] 直接映射为多项式乘法；差卷积
 *      res[k] = sum_{i-j=k} a[i] * b[j] 通过反转序列 b 转化为标准加法卷积，卷积第 p
 *      项对应差值 k = p - (|b| - 1)。
 *
 * API:
 *     convolution(a, b)            — 计算普通加法卷积 a * b，复杂度 O((n+m)
 *                                     log(n+m))。
 *     difference_convolution(a, b) — 计算差卷积，返回 a 与 reverse(b) 的卷积结果。
 *                                     复杂度 O((n+m) log(n+m))。
 *
 * Notes:
 *      1. 差卷积结果中，下标 |b| - 1 对应 i - j = 0，下标 |b| - 1 + d 对应 i - j =
 *         d。
 *
 * Related:
 *      数学/多项式/0-base/Poly·多项式全家桶.hpp: 底层多项式类。
 */

namespace poly_ext {

// 普通卷积 (普通多项式乘法)
// res[k] = sum_{i+j = k} a[i] * b[j]
Poly convolution(const Poly& a, const Poly& b) {
    return a * b;
}

/**
 * 差卷积
 * 计算 res[k] = sum_{i-j = k} a[i] * b[j]
 * 
 * 原理:
 * 令 j' = m - 1 - j (其中 m = b.size())
 * 则 j = m - 1 - j'
 * i - j = k  =>  i - (m - 1 - j') = k  =>  i + j' = k + m - 1
 * 
 * 也就是计算 a 和 reverse(b) 的卷积。
 * 卷积结果的第 p 项 (p = i + j') 对应原式中 k = p - (m - 1) 的系数。
 * 
 * 下标对应关系:
 * res[p] 对应差值为 k = p - (b.size() - 1) 的项
 * 特别地:
 * - res[b.size() - 1] 对应 k = 0 (i = j) 的项
 * - res[b.size() - 1 + d] 对应 k = d (i = j + d) 的项
 * - res[b.size() - 1 - d] 对应 k = -d (i = j - d) 的项
 */
Poly difference_convolution(const Poly& a, const Poly& b) {
    Poly b_rev = b;
    std::reverse(b_rev.begin(), b_rev.end());
    return a * b_rev;
}

} // namespace poly_ext
