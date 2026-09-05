#include "aizalib.h"
/*
 * EuclideanLike·类欧几里得算法
 *
 * Overview:
 *     类欧几里得算法，用于在对数时间内计算下取整一次函数和 f(a, b, c, n) = ∑_{i=0}^n ⌊(ai + b)/c⌋。
 *
 * API:
 *     sim_euclid::solve(a, b, c, n) — 计算 ∑_{i=0}^n ⌊(ai + b)/c⌋，复杂度 O(log(min(a, c, n)))
 *
 * Notes:
 *     1. 要求 a, b, c >= 0 且 c > 0, n >= 0。
 *     2. 递归通过类辗转相除法将坐标系翻转，复杂度与 gcd 相同。
 */

namespace sim_euclid {
    i64 solve(i64 a, i64 b, i64 c, i64 n) {
        i64 n2 = n * (n + 1) / 2;
        if (a >= c || b >= c) 
            return solve(a % c, b % c, c, n) + (a / c) * n2 + (b / c) * (n + 1);
        i64 m = (a * n + b) / c;
        if (!m) return 0;
        return m * n - solve(c, c - b - 1, a, m - 1);
    }
}
