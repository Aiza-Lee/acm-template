#pragma once
#include "../0-base/Poly·多项式全家桶.hpp"
/*
 * LagrangeInversion·拉格朗日反演
 *
 * Overview:
 *     给定多项式 F(x) 满足 F(0) = 0, F'(0) != 0，求其复合逆 G(x) (满足 G(F(x)) = x) 的第 n 项系数 [x^n] G(x)，
 *     以及扩展拉格朗日反演 [x^n] H(G(x))。
 *
 * API:
 *     coeff(Poly F, int n)                     — 计算复合逆 G(x) 的第 n 项系数 [x^n] G(x)
 *     generalized_coeff(Poly F, Poly H, int n) — 计算复合函数 H(G(x)) 的第 n 项系数 [x^n] H(G(x))
 *
 * Notes:
 *     1. 要求 F[0] == 0 且 F[1] != 0。当 F[1] != 1 时已内置归一化系数缩放。
 *     2. 复杂度为 O(n log n)，依托多项式全家桶求逆与幂函数实现。
 */

namespace Lagrange {

int coeff(Poly F, int n) {
    if (n == 0) return 0;
    if (F.empty()) return 0;
    // F(x) / x
    F.erase(F.begin());
    int f1 = (F.empty() ? 0 : F[0]);
    if (f1 == 0) return 0;
    int inv_f1 = inv(f1);
    for (int& x : F) x = mul(x, inv_f1);

    // (F(x)/x)^(-n)
    int k = mod(-n);
    F = F.pow(k, n);

    if (n - 1 >= (int)F.size()) return 0;
    int factor = mul(inv(n), Poly::Core::fp(inv_f1, n));
    return mul(F[n - 1], factor);
}

int generalized_coeff(Poly F, Poly H, int n) {
    if (n == 0) return H.empty() ? 0 : H[0];
    if (F.empty()) return 0;

    F.erase(F.begin());
    int f1 = (F.empty() ? 0 : F[0]);
    if (f1 == 0) return 0;
    int inv_f1 = inv(f1);
    for (int& x : F) x = mul(x, inv_f1);

    int k = mod(-n);
    F = F.pow(k, n);

    H = H.deriv();
    F = F * H;

    if (n - 1 >= (int)F.size()) return 0;
    int factor = mul(inv(n), Poly::Core::fp(inv_f1, n));
    return mul(F[n - 1], factor);
}

} // namespace Lagrange
