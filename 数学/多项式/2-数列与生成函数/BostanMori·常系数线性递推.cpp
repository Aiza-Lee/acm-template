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
 *      同时提供转置 Bostan-Mori（SliceCoeff）与区间 Bostan-Mori 扩展，支持在 O(d
 *      log d log n + (d + m) log(d + m)) 时间内批量求解区间 [n, n+m-1] 的连续项。
 *
 * API:
 *     bostan_mori(n, P, Q)                — 求解单项 [x^n] P(x)/Q(x)，复杂度 O(d
 *                                            log d log n)。
 *     linear_recurrence(n, c, a)          — 求解满足 a_n = sum_{i=1..d} c_i *
 *                                            a_{n-i} 的单项 a_n，复杂度 O(d log d
 *                                            log n)。
 *     transposed_bostan_mori(n, Q)        — 转置 Bostan-Mori 算法，求解 [x^n ...
 *                                            x^{n+d-1}] 1/Q(x)，复杂度 O(d log d
 *                                            log n)。
 *     bostan_mori_range(n, m, P, Q)       — 区间 Bostan-Mori，求解 [x^n ...
 *                                            x^{n+m-1}] P(x)/Q(x)，复杂度 O(d log d
 *                                            log n + (d + m) log(d + m))。
 *     linear_recurrence_range(n, m, c, a) — 求解满足递推关系的连续区间项 a_n, ...,
 *                                            a_{n+m-1}，复杂度 O(d log d log n + (d
 *                                            + m) log(d + m))。
 *
 * Notes:
 *      1. 要求 Q[0] != 0。
 *      2. n 为 64 位整数 (i64)，下标 0-based；m 为求解区间长度。
 *      3. 当 deg(P) >= deg(Q) 时，区间算法自动通过多项式带余除法处理多项式整除部分。
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

Poly transposed_bostan_mori(i64 n, Poly Q) {
    while (!Q.empty() && Q.back() == 0) Q.pop_back();
    if (Q.empty() || Q[0] == 0) return {};
    int d = (int)Q.size() - 1;
    if (d <= 0) return {};
    i64 N = n + d - 1;
    if (N < 0) return Poly(d, 0);

    std::vector<std::pair<i64, Poly>> stack;
    i64 cur_N = N;
    Poly cur_Q = Q;
    while (cur_N > 0) {
        Poly Q_neg = cur_Q;
        for (int i = 1; i < (int)Q_neg.size(); i += 2)
            Q_neg[i] = sub(0, Q_neg[i]);
        auto A = cur_Q * Q_neg;
        Poly V(d + 1);
        rep(i, 0, d) V[i] = (2 * i < (int)A.size()) ? A[2 * i] : 0;
        stack.emplace_back(cur_N, std::move(Q_neg));
        cur_N >>= 1;
        cur_Q = std::move(V);
    }

    Poly W(d, 0);
    W[d - 1] = inv(cur_Q[0]);

    while (!stack.empty()) {
        auto [n_val, Q_neg] = std::move(stack.back());
        stack.pop_back();

        Poly S(2 * d, 0);
        if (n_val % 2 == 0) {
            rep(i, 0, d - 1) S[2 * i + 1] = W[i];
        } else {
            rep(i, 0, d - 1) S[2 * i] = W[i];
        }

        auto B = Q_neg * S;
        W.assign(d, 0);
        rep(i, 0, d - 1) {
            if (d + i < (int)B.size()) W[i] = B[d + i];
        }
    }
    return W;
}

Poly bostan_mori_range(i64 n, int m, Poly P, Poly Q) {
    if (m <= 0) return {};
    while (!Q.empty() && Q.back() == 0) Q.pop_back();
    if (Q.empty() || Q[0] == 0) return Poly(m, 0);

    int d = (int)Q.size() - 1;
    if (d == 0) {
        int inv_q0 = inv(Q[0]);
        Poly ans(m, 0);
        rep(i, 0, m - 1) {
            if (n + i < (i64)P.size()) ans[i] = mul(P[n + i], inv_q0);
        }
        return ans;
    }

    Poly ans(m, 0);
    Poly Prem = P;
    if ((int)P.size() >= (int)Q.size()) {
        auto [R, rem] = P.div_mod(Q);
        Prem = rem;
        rep(i, 0, m - 1) {
            if (n + i < (i64)R.size()) ans[i] = R[n + i];
        }
    }
    while (!Prem.empty() && Prem.back() == 0) Prem.pop_back();
    if (Prem.empty()) return ans;

    if (n == 0) {
        auto invQ = Q.inverse(m);
        auto series = Prem * invQ;
        rep(i, 0, m - 1) {
            if (i < (int)series.size()) inc(ans[i], series[i]);
        }
        return ans;
    }

    auto invQ_slice = transposed_bostan_mori(n, Q);
    auto prod = Q * invQ_slice;
    prod.resize(d);
    auto Pk = (prod * Prem) % Q;

    auto invQ = Q.inverse(m);
    auto series = Pk * invQ;
    rep(i, 0, m - 1) {
        if (i < (int)series.size()) inc(ans[i], series[i]);
    }
    return ans;
}

std::vector<int> linear_recurrence_range(
    i64 n, int m, const std::vector<int>& c, const std::vector<int>& a
) {
    if (m <= 0) return {};
    int d = c.size();
    if (d == 0) return std::vector<int>(m, 0);

    Poly Q(d + 1);
    Q[0] = 1;
    rep(i, 0, d - 1)
        Q[i + 1] = sub(0, c[i]);

    Poly A(a);
    if ((int)A.size() < d) A.resize(d, 0);
    auto P = A * Q;
    P.resize(d);

    Poly res = bostan_mori_range(n, m, P, Q);
    return std::vector<int>(res.begin(), res.end());
}

} // namespace poly_ext
