#include "aizalib.h"

/*
 * Euler's Totient Function (单点欧拉函数)
 *
 * Overview:
 *      计算正整数 n 的欧拉函数 phi(n)，即不超过 n 且与 n 互质的正整数个数。
 *      利用算术基本定理 n = prod(p_i^{k_i})，根据积性性质得计算公式 phi(n) = n *
 *      prod(1 - 1 / p_i) = prod(p_i^{k_i - 1} * (p_i - 1))。
 *      采用试除法在 O(sqrt(n)) 时间内找出全部互异质因子。
 *
 * API:
 *     phi(n) — 计算单个整数 n 的欧拉函数 phi(n)。复杂度 O(sqrt(n)) 时间，O(1)
 *               空间。
 *
 * Notes:
 *      1. 要求 n >= 1。特别地 phi(1) = 1。
 *      2. 多次查询或批量计算请使用线性筛 EulerSieve。
 *
 * Related:
 *      数学/数论/EulerSieve·线性筛.cpp: O(N) 批量预处理积性函数。
  */
struct EulerPhi {
    static i64 phi(i64 n) {
        AST(n >= 1);
        i64 res = n;
        for (i64 p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
            if (n % p) continue;
            res = res / p * (p - 1);
            while (n % p == 0) n /= p;
        }
        if (n > 1) res = res / n * (n - 1);
        return res;
    }
};
