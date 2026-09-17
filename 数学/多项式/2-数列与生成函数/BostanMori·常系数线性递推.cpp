#include "aizalib.h"
#include "0-base/Poly·多项式全家桶.hpp"
/*
 * Bostan-Mori Algorithm (常系数线性递推与有理分式展开)
 *
 * Overview:
 *      求解有理分式系数 [x^n] P(x)/Q(x) 以及 d 阶常系数线性递推的第 n 项。
 *      利用倍增转换 P(x)/Q(x) = (P(x)*Q(-x)) / (Q(x)*Q(-x))，分母变为关于 x^2
 *      的偶多项式，根据 n 的奇偶性提取分子的奇偶项，每次将求解规模 n 折半，在 O(d
 *      log d log n) 时间内求出远项 a_n，大幅优于传统多项式取模。
 *
 * API:
 *     bostan_mori(n, P, Q)       — 求解 [x^n] P(x)/Q(x)，复杂度 O(d log d log n)。
 *     linear_recurrence(n, c, a) — 求解满足 a_n = sum_{i=1..d} c_i * a_{n-i} 的第
 *                                   n 项，其中 c = {c1, ..., cd}, a = {a0, ...,
 *                                   a_{d-1}}。复杂度 O(d log d log n)。
 *
 * Notes:
 *      1. 要求 Q[0] != 0。
 *      2. n 为 64 位整数 (i64)，下标 0-based。
 *
 * Related:
 *      数学/多项式/2-数列与生成函数/BerlekampMassey·最短线性递推.cpp:
 *      求解最短递推式。
 */

namespace poly_ext {

int bostan_mori(i64 n, Poly P, Poly Q) {
    while (n > 0) {
        auto Q_neg = Q;
        for (int i = 1; i < (int)Q_neg.size(); i += 2)
            Q_neg[i] = sub(0, Q_neg[i]);

        auto U = P * Q_neg, V = Q * Q_neg;
        
        P.clear();
        for (int i = (n & 1); i < (int)U.size(); i += 2)
            P.push_back(U[i]);

        Q.clear();
        for (int i = 0; i < (int)V.size(); i += 2)
            Q.push_back(V[i]);

        n >>= 1;
    }
    if (Q.empty() || Q[0] == 0) return 0;
    return mul(P.empty() ? 0 : P[0], inv(Q[0]));
}

int linear_recurrence(
    i64 n, const std::vector<int>& c, const std::vector<int>& a
) {
    int d = c.size();
    if (n < d) return a[n];

    Poly Q(d + 1);
    Q[0] = 1;
    rep(i, 0, d - 1)
        Q[i + 1] = sub(0, c[i]);

    Poly A(a);
    auto P = A * Q;
    P.resize(d);

    return bostan_mori(n, P, Q);
}

} // namespace poly_ext
