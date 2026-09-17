#include "aizalib.h"

/*
 * Pohlig-Hellman Algorithm (光滑阶离散对数算法)
 *
 * Overview:
 *      在素数模 p 下求解离散对数 a^x = b (mod p)。
 *      当群阶 p - 1 或 a 的乘法阶 ord(a) 的所有质因子均较小（光滑数 smooth
 *      number）时，将离散对数问题分解为各个素数幂子群 a^(n / q_i^{e_i}) 上的子问题。
 *      在每个素数幂子群中通过 p-adic 逐位提取确定 x mod q_i^{e_i}，
 *      最后利用中国剩余定理 (CRT) 组合得到模 ord(a) 的全局解。
 *
 * API:
 *     solve(a, b, p)                 — 自动试除分解 p-1 求解 a^x = b (mod p)
 *                                       的最小非负解，无解返回 -1。复杂度 O(sqrt(p)
 *                                       + sum(e_i * q_i)) 时间。
 *     solve(a, b, p, factors_of_pm1) — 传入 p-1 的质因数分解 {(q, e)} 求解。复杂度
 *                                       O(sum(e_i * q_i)) 时间。
 *
 * Notes:
 *      1. 要求 p 为质数且 gcd(a, p) = 1。
 *      2. 极适用于 NTT 模数（如 998244353 = 2^23 * 7 * 17 + 1 等）。
 *      3. 若 p - 1 含有大质因子（如 10^9 + 7），建议改用 BSGS。
 *
 * Related:
 *      数学/数论/BSGS·求离散对数.cpp: 大质因子子群或一般离散对数。
 *      数学/数论/CRT·中国剩余定理.cpp: 子群结果合并。
 */
struct PohligHellman {
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

    // 小规模子群 DLP: g 的阶为 q，求 g^x = h (mod p)，x in [0, q)
    // 当 q 较大（如 q > 1e6）时可替换为 BSGS::solve(g, h, p) 降至 O(sqrt(q))
    static i64 _dlp_small(i64 g, i64 h, i64 q, i64 p) {
        if (h == 1) return 0;
        i64 cur = 1;
        for (i64 x = 1; x < q; ++x) {
            cur = (i128)cur * g % p;
            if (cur == h) return static_cast<i64>(x);
        }
        return -1;
    }

    // 求 x mod q^e，逐位提取 base-q 数字
    static i64 _solve_prime_power(i64 a, i64 b, i64 n, i64 q, int e, i64 p) {
        i64 x_q = 0, q_pow = 1;
        i64 g = _pow_mod(a, n / q, p);  // ord(g) = q
        i64 inv_a = _inv_mod(a, p);

        for (int k = 0; k < e; ++k) {
            i64 base = (i128)_pow_mod(inv_a, x_q, p) * b % p;
            i64 h = _pow_mod(base, n / (q_pow * q), p);
            i64 d = _dlp_small(g, h, q, p);
            if (d == -1) return -1;
            x_q += d * q_pow;
            q_pow *= q;
        }
        return x_q;
    }

    static std::vector<std::pair<i64, int>> _factor(i64 n) {
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

    static i64 solve(i64 a, i64 b, i64 p) {
        AST(p > 1);
        a = _norm(a, p), b = _norm(b, p);
        if (b == 1) return 0;
        if (a == 0) return b == 0 ? 1 : -1;
        if (std::gcd(a, p) != 1) return -1;
        return solve(a, b, p, _factor(p - 1));
    }

    static i64 solve(
        i64 a, i64 b, i64 p,
        const std::vector<std::pair<i64, int>> &factors_of_pm1
    ) {
        AST(p > 1);
        a = _norm(a, p), b = _norm(b, p);
        if (b == 1) return 0;
        if (a == 0) return b == 0 ? 1 : -1;
        if (std::gcd(a, p) != 1) return -1;

        // Compute ord(a) and its factorization from factors of p-1
        i64 n = p - 1;
        std::vector<std::pair<i64, int>> ord_factors;
        for (auto [q, e] : factors_of_pm1) {
            int cnt = e;
            while (cnt > 0 && _pow_mod(a, n / q, p) == 1) n /= q, --cnt;
            if (cnt > 0) ord_factors.emplace_back(q, cnt);
        }
        if (ord_factors.empty()) return -1;  // a ≡ 1 (mod p) and b ≠ 1, no solution

        i64 x = 0, M = 1;
        for (auto [q, e] : ord_factors) {
            i64 qe = 1;
            for (int k = 0; k < e; ++k) qe *= q;
            i64 x_i = _solve_prime_power(a, b, n, q, e, p);
            if (x_i == -1) return -1;

            // CRT: combine x (mod M) with x_i (mod qe)
            i64 diff = _norm(x_i - x, qe);
            i64 t, y;
            i64 g = _exgcd(M, qe, t, y);
            if (diff % g) return -1;
            t = (i128)_norm(t, qe) * (diff / g) % qe;
            x = x + M * t;
            M = M / g * qe;
            x = _norm(x, M);
        }
        return x;
    }
};
