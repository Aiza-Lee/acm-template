#include "aizalib.h"

/*
 * Multiplicative Order (乘法阶)
 *
 * Overview:
 *      在 gcd(a, mod) = 1 时求 a 在模 mod 意义下的乘法阶 ord_m(a)，即满足 a^k = 1
 *      (mod mod) 的最小正整数 k。
 *      由欧拉定理可知阶必整除 phi(mod)。首先通过试除分解求出 phi(mod)，再分解
 *      phi(mod) 的质因子，从当前可能阶中贪心试除多余质因子，直至无法继续整除。
 *
 * API:
 *     order(a, mod) — 求 a 模 mod 的最小正阶。若 mod <= 1 或 gcd(a, mod) != 1
 *                      则返回 -1。复杂度 O(sqrt(mod) + log(mod) * 质因子个数) 时间，
 *                      O(1) 空间。
 *
 * Notes:
 *      1. 要求 gcd(a, mod) = 1 否则阶不存在。
 *      2. 模数极大时可将试除法替换为 PollardRho 分解。
 *
 * Related:
 *      数学/数论/EulerPhi·单点欧拉函数.cpp: 计算欧拉函数 phi(n)。
 *      数学/数论/PrimitiveRoot·原根.cpp: 原根定义为阶等于 phi(m) 的元素。
 *      数学/数论/ModMultiplicativeGroup·模乘法群.cpp: 任意模数乘法群基底与坐标化求阶。
 */
struct MultiplicativeOrder {
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

    static std::vector<i64> _factor_distinct(i64 n) {
        std::vector<i64> res;
        for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
            if (n % p) continue;
            res.emplace_back(p);
            while (n % p == 0) n /= p;
        }
        if (n > 1) res.emplace_back(n);
        return res;
    }

    static i64 _phi(i64 n) {
        i64 res = n;
        for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
            if (n % p) continue;
            res = res / p * (p - 1);
            while (n % p == 0) n /= p;
        }
        if (n > 1) res = res / n * (n - 1);
        return res;
    }

    static i64 order(i64 a, i64 mod) {
        if (mod <= 1 || std::gcd(a, mod) != 1) return -1;
        i64 ord = _phi(mod);
        for (i64 p : _factor_distinct(ord)) {
            while (ord % p == 0 && _pow_mod(a, ord / p, mod) == 1) ord /= p;
        }
        return ord;
    }
};
