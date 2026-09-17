#include "aizalib.h"

/*
 * Tonelli-Shanks Algorithm (二次剩余与模平方根)
 *
 * Overview:
 *      在奇素数模 p 下求解二次同余方程 x^2 = a (mod p)。
 *      利用欧拉判别准则 (勒让德符号 a^((p-1)/2) mod p) 判断 a 是否为模 p
 *      的二次剩余。
 *      若 p = 3 (mod 4)，解为直接公式 a^((p+1)/4) mod p；对于一般奇素数，分解 p - 1
 *      = q * 2^s，寻找二次非剩余 z，在 2-Sylow 子群上不断折半消去误差项，
 *      快速逼近精确根。
 *
 * API:
 *     legendre(a, p) — 计算勒让德符号 (a/p)，返回 1 (二次剩余)、-1 (二次非剩余) 或
 *                       0 (a = 0 mod p)。复杂度 O(log p) 时间。
 *     solve(a, p)    — 求解 x^2 = a (mod p)，返回从小到大排序的解对 {x1, x2}。
 *                       无解返回 {-1, -1}。复杂度 O(log^2 p) 时间，O(1) 空间。
 *
 * Notes:
 *      1. 要求 p 为素数。若 p = 2 单独处理；若 a = 0 (mod p) 返回 {0, 0}。
 *      2. 若存在解 x，则另一解为 (p - x) mod p。
 *
 * Related:
 *      数学/数论/DiscreteRoot·离散开根.cpp: 一般高次同余方程求解。
 */
struct TonelliShanks {
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

    static int legendre(i64 a, i64 p) {
        a = _norm(a, p);
        if (!a) return 0;
        i64 x = _pow_mod(a, (p - 1) / 2, p);
        return x == 1 ? 1 : -1;
    }

    static std::pair<i64, i64> solve(i64 a, i64 p) {
        AST(p >= 2);
        a = _norm(a, p);
        if (p == 2) return {a, a};
        if (!a) return {0, 0};
        if (legendre(a, p) != 1) return {-1, -1};
        if (p % 4 == 3) {
            i64 x = _pow_mod(a, (p + 1) / 4, p);
            i64 y = p - x;
            if (x > y) std::swap(x, y);
            return {x, y};
        }

        i64 q = p - 1, s = 0;
        while (!(q & 1)) q >>= 1, ++s;
        i64 z = 2;
        while (legendre(z, p) != -1) ++z;
        i64 c = _pow_mod(z, q, p);
        i64 t = _pow_mod(a, q, p);
        i64 r = _pow_mod(a, (q + 1) / 2, p);
        i64 m = s;
        while (t != 1) {
            i64 i = 1, tt = (i128)t * t % p;
            while (tt != 1) tt = (i128)tt * tt % p, ++i;
            i64 b = _pow_mod(c, 1LL << (m - i - 1), p);
            r = (i128)r * b % p;
            c = (i128)b * b % p;
            t = (i128)t * c % p;
            m = i;
        }
        i64 x = r, y = p - r;
        if (x > y) std::swap(x, y);
        return {x, y};
    }
};
