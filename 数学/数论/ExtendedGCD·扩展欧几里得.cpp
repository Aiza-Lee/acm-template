#include "aizalib.h"

/*
 * Extended Euclidean Algorithm (扩展欧几里得算法)
 *
 * Overview:
 *      求解裴蜀等式 a * x + b * y = gcd(a, b) 的一组整数特解 (x, y)。
 *      在此基础上提供模逆元求解和二元一次不定方程 a * x + b * y = c
 *      的特解与通解步长。
 *
 * API:
 *     exgcd(a, b, x, y)                   — 求 a*x + b*y = gcd(a, b) 的特解，返回
 *                                            gcd(a, b)。复杂度 O(log(min(|a|,
 *                                            |b|))) 时间与递归栈空间。
 *     inv_mod(a, mod, inv)                — 求 a 在模 mod 下的逆元写入 inv，
 *                                            若逆元存在返回 true，否则 false。
 *                                            复杂度 O(log(min(|a|, |mod|)))。
 *     solve_linear(a, b, c, x, y, dx, dy) — 解 a*x + b*y = c，若有解将一组特解写入
 *                                            (x, y)，通解步长写入 (dx, dy)，
 *                                            使得通解为 x = x0 + k*dx, y = y0 -
 *                                            k*dy 并返回 true；无解返回 false。
 *
 * Notes:
 *      1. a*x + b*y = c 有整数解当且仅当 gcd(a, b) | c。
 *      2. 要求特解 x, y 缩放后在 i64 范围内。
 */
struct ExtendedGCD {
    static i64 exgcd(i64 a, i64 b, i64 &x, i64 &y) {
        if (!b) {
            x = a >= 0 ? 1 : -1;
            y = 0;
            return std::abs(a);
        }
        i64 x1, y1, g = exgcd(b, a % b, x1, y1);
        x = y1;
        y = x1 - a / b * y1;
        return g;
    }

    static bool inv_mod(i64 a, i64 mod, i64 &inv) {
        AST(mod > 0);
        i64 x, y;
        if (exgcd(a, mod, x, y) != 1) return false;
        inv = (x % mod + mod) % mod;
        return true;
    }

    static bool solve_linear(i64 a, i64 b, i64 c, i64 &x, i64 &y, i64 &dx, i64 &dy) {
        if (!a && !b) {
            if (c) return false;
            x = y = dx = dy = 0;
            return true;
        }
        i64 g = exgcd(a, b, x, y);
        if (c % g) return false;
        i128 k = c / g;
        i128 xx = (i128)x * k, yy = (i128)y * k;
        AST(std::numeric_limits<i64>::min() <= xx &&
            xx <= std::numeric_limits<i64>::max());
        AST(std::numeric_limits<i64>::min() <= yy &&
            yy <= std::numeric_limits<i64>::max());
        x = (i64)xx, y = (i64)yy;
        dx = b / g, dy = a / g;
        return true;
    }
};
