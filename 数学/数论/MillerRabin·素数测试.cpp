#include "aizalib.h"

/*
 * Miller-Rabin Primality Test (Miller-Rabin 素数测试)
 *
 * Overview:
 *      基于费马小定理与二次探测定理的确定性快速素性测试。
 *      对于 64 位整数范围内的 n，选用固定的 7 个优质测试基底 {2, 325, 9375, 28178,
 *      450775, 9780504, 1795265022}，配合小质数预筛，可在 100% 确定性下判定 n
 *      是否为素数，无需随机化。
 *
 * API:
 *     is_prime(n) — 判定 64 位整数 n 是否为素数。复杂度 O(k log n) 时间，其中 k=7
 *                    为底数个数。
 *
 * Notes:
 *      1. 若 n < 2 返回 false。
 *      2. 内部使用 i128 处理模乘防止溢出。
 *
 * Related:
 *      数学/数论/PollardSRho·寻找质因数.cpp: 基于 Miller-Rabin 的大数质因数分解。
 */
struct MillerRabin {
    static i64 _mul_mod(i64 a, i64 b, i64 mod) {
        return (i128)a * b % mod;
    }

    static i64 _pow_mod(i64 a, i64 b, i64 mod) {
        i64 res = 1;
        for (a %= mod; b; b >>= 1, a = _mul_mod(a, a, mod))
            if (b & 1) res = _mul_mod(res, a, mod);
        return res;
    }

    static bool _check(i64 a, i64 s, i64 d, i64 n) {
        if (a % n == 0) return true;
        i64 x = _pow_mod(a, d, n);
        if (x == 1 || x == n - 1) return true;
        rep(i, 1, (int)s - 1) {
            x = _mul_mod(x, x, n);
            if (x == n - 1) return true;
        }
        return false;
    }

    static bool is_prime(i64 n) {
        if (n < 2) return false;
        for (i64 p : {2LL, 3LL, 5LL, 7LL, 11LL, 13LL, 17LL, 19LL, 23LL, 29LL,
                      31LL, 37LL}) {
            if (n % p == 0) return n == p;
        }
        i64 d = n - 1, s = 0;
        while (!(d & 1)) d >>= 1, ++s;
        for (i64 a : {2LL, 325LL, 9375LL, 28178LL, 450775LL, 9780504LL,
                      1795265022LL}) {
            if (!_check(a, s, d, n)) return false;
        }
        return true;
    }
};
