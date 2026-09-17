#include "aizalib.h"

/*
 * Extended Lucas Theorem (扩展 Lucas 定理)
 *
 * Overview:
 *      计算任意正整数模数 mod 下的组合数 C(n, m) mod mod。
 *      将模数质因数分解为 mod = prod(p_i^{k_i})，分别计算每个质数幂模数下的 C(n, m)
 *      mod p_i^{k_i}，最后利用扩展中国剩余定理 (ExCRT) 合并出模 mod 的解。
 *      计算 n! mod p^k 时提取所有质因子 p，利用阶乘循环节和递归分治计算与 p
 *      互质部分。
 *
 * API:
 *     C(n, m, mod) — 计算组合数 C(n, m) mod mod。复杂度 O(sum(p_i^{k_i}) + t log
 *                     n) 时间，O(max(p_i^{k_i})) 空间，其中 t 为 mod
 *                     的不同质因子数。
 *
 * Notes:
 *      1. 要求 mod > 0。若 m < 0 或 m > n 返回 0。若 mod = 1 返回 0。
 *      2. 适用于各 p_i^{k_i} 较小（通常 <= 10^6）的一般模数组合数计算。
 *
 * Related:
 *      数学/数论/ExtendedCRT·扩展中国剩余定理.cpp: 用于合并各同余方程。
 */
struct ExLucas {
    static i64 _norm(i64 x, i64 mod) {
        x %= mod;
        return x < 0 ? x + mod : x;
    }

    static i64 _pow_mod(i64 a, i64 b, i64 mod) {
        i64 res = 1;
        for (a = _norm(a, mod); b; b >>= 1, a = (i128)a * a % mod)
            if (b & 1) res = (i128)res * a % mod;
        return res;
    }

    static i64 _exgcd(i64 a, i64 b, i64 &x, i64 &y) {
        if (!b) {
            x = a >= 0 ? 1 : -1;
            y = 0;
            return std::abs(a);
        }
        i64 x1, y1, g = _exgcd(b, a % b, x1, y1);
        x = y1;
        y = x1 - a / b * y1;
        return g;
    }

    static i64 _inv_mod(i64 a, i64 mod) {
        i64 x, y;
        i64 g = _exgcd(a, mod, x, y);
        AST(g == 1);
        return _norm(x, mod);
    }

    static std::vector<std::pair<i64, int>> _factor_count(i64 n) {
        std::vector<std::pair<i64, int>> res;
        for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
            if (n % p) continue;
            int c = 0;
            while (n % p == 0) n /= p, ++c;
            res.emplace_back(p, c);
        }
        if (n > 1) res.emplace_back(n, 1);
        return res;
    }

    static i64 _count_p(i64 n, i64 p) {
        i64 res = 0;
        while (n) n /= p, res += n;
        return res;
    }

    struct PrimePowerBinom {
        i64 p, pk;
        std::vector<i64> pre;

        PrimePowerBinom(i64 p, int c) : p(p), pk(1) {
            rep(i, 1, c) pk *= p;
            pre.assign((size_t)pk + 1, 1);
            for (i64 i = 1; i <= pk; ++i) {
                pre[(size_t)i] = pre[(size_t)(i - 1)];
                if (i % p) pre[(size_t)i] = (i128)pre[(size_t)i] * i % pk;
            }
        }

        i64 _fact(i64 n) const {
            if (!n) return 1;
            i64 res = ExLucas::_pow_mod(pre[(size_t)pk], n / pk, pk);
            res = (i128)res * pre[(size_t)(n % pk)] % pk;
            return (i128)res * _fact(n / p) % pk;
        }

        i64 C(i64 n, i64 m) const {
            if (m < 0 || m > n) return 0;
            i64 e = ExLucas::_count_p(n, p) - ExLucas::_count_p(m, p)
                - ExLucas::_count_p(n - m, p);
            i64 a = _fact(n);
            i64 b = _fact(m);
            i64 c = _fact(n - m);
            i64 res = a;
            res = (i128)res * ExLucas::_inv_mod(b, pk) % pk;
            res = (i128)res * ExLucas::_inv_mod(c, pk) % pk;
            res = (i128)res * ExLucas::_pow_mod(p, e, pk) % pk;
            return res;
        }
    };

    static std::pair<i64, i64> _merge(i64 r0, i64 m0, i64 r1, i64 m1) {
        r0 = _norm(r0, m0), r1 = _norm(r1, m1);
        i64 x, y, g = _exgcd(m0, m1, x, y), d = r1 - r0;
        if (d % g) return {-1, -1};
        i64 mod = m1 / g;
        i64 k = (i128)(d / g) * x % mod;
        i128 lcm = (i128)m0 / g * m1;
        AST(lcm <= std::numeric_limits<i64>::max());
        i64 nm = (i64)lcm;
        i64 nr = _norm(r0 + (i128)m0 * k % nm, nm);
        return {nr, nm};
    }

    static i64 C(i64 n, i64 m, i64 mod) {
        if (m < 0 || m > n) return 0;
        if (mod == 1) return 0;
        auto fac = _factor_count(mod);
        i64 r0 = 0, m0 = 1;
        for (auto [p, c] : fac) {
            PrimePowerBinom cur(p, c);
            auto [nr, nm] = _merge(r0, m0, cur.C(n, m), cur.pk);
            r0 = nr, m0 = nm;
        }
        return r0;
    }
};
