#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"
/*
 * Euler Transform (多项式欧拉变换与无限乘积展开)
 *
 * Overview:
 *      计算无标号组合对象的欧拉变换 B(x) = prod_{i=1}^inf (1 - x^i)^{-a_i}。
 *      两边取对数 ln B(x) = sum_{i=1}^inf a_i * sum_{k=1}^inf (x^{ik} / k)，
 *      通过调和级数枚举倍数在 O(n log n) 时间内构造出 ln B(x)，再通过多项式 Exp
 *      还原 B(x)。
 *      同理可快速计算型如 prod (1 - x^{a_i}) 的多项式连乘积。
 *
 * API:
 *     euler_transform(a, n)            — 计算 a 的欧拉变换，返回长为 n+1 的多项式，
 *                                         复杂度 O(n log n)。
 *     fast_multiply_one_minus_xa(a, n) — 快速展开 prod_{v in a} (1 - x^v) mod
 *                                         x^{n+1}，复杂度 O(n log n)。
 *
 * Notes:
 *      1. 要求 a[0] = 0。
 *      2. 输出多项式截断到 x^n（长度为 n + 1）。
 *
 * Related:
 *      数学/多项式/0-base/Poly·多项式全家桶.hpp: 多项式 Exp 算子。
 */

namespace poly_ext {

Poly euler_transform(const Poly& a, int n) {
    std::vector<int> sum(n + 1, 0);
    rep(i, 1, n) {
        if (i < (int)a.size() && a[i]) {
            int val = mul(a[i], i);
            for (int j = i; j <= n; j += i) inc(sum[j], val);
        }
    }

    Poly ln_b(n + 1);
    rep(i, 1, n) {
        ln_b[i] = mul(sum[i], inv(i));
    }

    return ln_b.exp(n + 1);
}

Poly fast_multiply_one_minus_xa(const std::vector<int>& a, int n) {
    Poly cnt(n + 1);
    for (int x : a) {
        if (x <= n) dec(cnt[x], 1);
    }
    return euler_transform(cnt, n);
}

} // namespace poly_ext
