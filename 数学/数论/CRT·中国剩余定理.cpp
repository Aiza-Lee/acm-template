#include "aizalib.h"

/*
 * Chinese Remainder Theorem (CRT)
 *
 * Overview:
 *      求解模数两两互质的线性同余方程组 x = r[i] (mod m[i])。
 *      令总模数 M = prod(m[i]), M[i] = M / m[i]，由于 gcd(M[i], m[i]) = 1，
 *      求出模逆元 t[i] = M[i]^(-1) (mod m[i]) 后，唯一特解为 sum(r[i] * M[i] *
 *      t[i]) (mod M)。
 *
 * API:
 *     solve(r, m) — 求解同余方程组并返回 [0, M) 内的唯一最小非负解。复杂度 O(n
 *                    log(max m)) 时间，O(1) 额外空间。
 *
 * Notes:
 *      1. 要求各 m[i] > 0 且两两互质，且乘积 M = prod(m[i]) 不超过 i64 上限。
 *      2. 若模数不两两互质，应使用 ExtendedCRT。
 *
 * Related:
 *      数学/数论/ExtendedCRT·扩展中国剩余定理.cpp:
 *      模数不保证两两互质时的一般方程组求解。
 */
struct CRT {
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

    static i64 solve(const std::vector<i64> &r, const std::vector<i64> &m) {
        AST(r.size() == m.size());
        int n = (int)r.size();
        if (!n) return 0;

        i128 M = 1;
        for (i64 v : m) {
            AST(v > 0);
            M *= v;
            AST(M <= std::numeric_limits<i64>::max());
        }
        i64 mod = (i64)M, ans = 0;
        rep(i, 0, n - 1) {
            i64 ri = _norm(r[i], m[i]), Mi = mod / m[i], x, y;
            i64 g = _exgcd(Mi, m[i], x, y);
            AST(g == 1);
            x = _norm(x, m[i]);
            ans = (ans + (i128)ri * Mi % mod * x) % mod;
        }
        return _norm(ans, mod);
    }
};
