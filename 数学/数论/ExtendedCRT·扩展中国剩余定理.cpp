#include "aizalib.h"

/*
 * Extended Chinese Remainder Theorem (ExCRT)
 *
 * Overview:
 *      求解模数不一定两两互质的线性同余方程组 x = r[i] (mod m[i])。
 *      采用增量法，已知前 k-1 个方程的解为 x = r_0 (mod m_0)，与第 k 个方程 x = r_1
 *      (mod m_1) 联立为 r_0 + k_0 * m_0 = r_1 (mod m_1)，即 k_0 * m_0 = r_1 - r_0
 *      (mod m_1)。用扩展欧几里得解线性同余方程，若 (r_1 - r_0) % gcd(m_0, m_1) != 0
 *      则无解；否则合并得到模 lcm(m_0, m_1) 的新同余式。
 *
 * API:
 *     merge(r0, m0, r1, m1) — 合并两个同余方程，返回 {新余数, 新模数}。无解返回
 *                              {-1, -1}。复杂度 O(log(max(m0, m1))) 时间。
 *     solve(r, m)           — 求解整个方程组，返回 {最小非负解, lcm}。无解返回 {-1,
 *                              -1}。复杂度 O(n log(lcm)) 时间，O(1) 额外空间。
 *
 * Notes:
 *      1. 要求各 m[i] > 0。
 *      2. 余数 r[i] 内部会自动正规化到 [0, m[i])。
 *      3. lcm 需在 i64 范围内，若可能溢出需做特别处理。
 *
 * Related:
 *      数学/数论/CRT·中国剩余定理.cpp: 模数两两互质时的经典 CRT。
 */
struct ExtendedCRT {
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

    static i64 _norm(i64 x, i64 mod) {
        x %= mod;
        return x < 0 ? x + mod : x;
    }

    static std::pair<i64, i64> merge(i64 r0, i64 m0, i64 r1, i64 m1) {
        AST(m0 > 0 && m1 > 0);
        r0 = _norm(r0, m0), r1 = _norm(r1, m1);
        i64 x, y, g = _exgcd(m0, m1, x, y), d = r1 - r0;
        if (d % g) return {-1, -1};

        i64 mod = m1 / g;
        i64 k = (i128)(d / g) * x % mod;
        i128 lcm = (i128)m0 / g * m1;
        AST(lcm <= std::numeric_limits<i64>::max());
        i64 nxt_m = (i64)lcm;
        i64 nxt_r = _norm(r0 + (i128)m0 * k % nxt_m, nxt_m);
        return {nxt_r, nxt_m};
    }

    static std::pair<i64, i64> solve(
        const std::vector<i64> &r, const std::vector<i64> &m
    ) {
        AST(r.size() == m.size());
        int n = (int)r.size();
        if (!n) return {0, 1};
        AST(m[0] > 0);
        i64 r0 = _norm(r[0], m[0]), m0 = m[0];
        rep(i, 1, n - 1) {
            auto [nr, nm] = merge(r0, m0, r[i], m[i]);
            if (nm == -1) return {-1, -1};
            r0 = nr, m0 = nm;
        }
        return {r0, m0};
    }
};
