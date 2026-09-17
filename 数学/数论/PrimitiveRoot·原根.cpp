#include "aizalib.h"

/*
 * Primitive Root (原根判定与求解)
 *
 * Overview:
 *      判定模数 m 是否存在原根，并求出其最小原根。
 *      由数论基本定理，模 m 存在原根当且仅当 m in {1, 2, 4, p^k, 2*p^k}，其中 p
 *      为奇素数。
 *      g 为模 m 的原根当且仅当 gcd(g, m) = 1 且对 phi(m) 的每个互异质因子 p，均有
 *      g^(phi(m) / p) != 1 (mod m)。最小原根通常很小（O(m^(1/4))），
 *      通过枚举测试即可快速找到。
 *
 * API:
 *     has_primitive_root(mod)   — 判断模数 mod 是否存在原根。复杂度 O(sqrt(mod))
 *                                  时间。
 *     is_primitive_root(g, mod) — 判断整数 g 是否为模 mod 的原根。复杂度
 *                                  O(sqrt(mod) + log(mod) * 质因子数) 时间。
 *     find(mod)                 — 求模 mod 的最小正原根。若不存在返回 -1；特别地
 *                                  mod=1 时返回 0。复杂度 O(sqrt(mod) + g_min *
 *                                  质因子数 * log mod) 时间。
 *
 * Notes:
 *      1. 要求 mod >= 1。
 *      2. 大模数时可将内部试除分解替换为 PollardRho。
 *
 * Related:
 *      数学/数论/MultiplicativeOrder·乘法阶.cpp: 元素阶的定义与计算。
 *      数学/数论/ModMultiplicativeGroup·模乘法群.cpp: 一般模数下的群结构与正交基底分解。
 */
struct PrimitiveRoot {
    static i64 _pow_mod(i64 a, i64 b, i64 mod) {
        i64 res = 1;
        for (a %= mod; b; b >>= 1, a = (i128)a * a % mod)
            if (b & 1) res = (i128)res * a % mod;
        return res;
    }

    static std::vector<std::pair<i64, int>> _factor_count(i64 n) {
        std::vector<std::pair<i64, int>> res;
        if (n <= 1) return res;
        for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
            if (n % p) continue;
            int c = 0;
            while (n % p == 0) n /= p, ++c;
            res.emplace_back(p, c);
        }
        if (n > 1) res.emplace_back(n, 1);
        return res;
    }

    static std::vector<i64> _factor_distinct(i64 n) {
        std::vector<i64> res;
        for (auto [p, _] : _factor_count(n)) res.emplace_back(p);
        return res;
    }

    static i64 _phi(i64 n) {
        i64 res = n;
        for (auto [p, _] : _factor_count(n)) res = res / p * (p - 1);
        return res;
    }

    static bool has_primitive_root(i64 mod) {
        AST(mod >= 1);
        if (mod == 1 || mod == 2 || mod == 4) return true;
        auto fac = _factor_count(mod);
        if (fac.size() == 1) return fac[0].first != 2;
        return fac.size() == 2 && fac[0].first == 2 && fac[0].second == 1
            && fac[1].first != 2;
    }

    static bool is_primitive_root(i64 g, i64 mod) {
        if (!has_primitive_root(mod) || std::gcd(g, mod) != 1) return false;
        if (mod == 1) return g == 0;
        i64 phi = _phi(mod);
        for (i64 p : _factor_distinct(phi)) {
            if (_pow_mod(g, phi / p, mod) == 1) return false;
        }
        return true;
    }

    static i64 find(i64 mod) {
        if (!has_primitive_root(mod)) return -1;
        if (mod == 1) return 0;
        i64 phi = _phi(mod);
        auto fac = _factor_distinct(phi);
        for (i64 g = 1; g < mod; ++g) {
            if (std::gcd(g, mod) != 1) continue;
            bool ok = true;
            for (i64 p : fac) {
                if (_pow_mod(g, phi / p, mod) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return g;
        }
        return -1;
    }
};
