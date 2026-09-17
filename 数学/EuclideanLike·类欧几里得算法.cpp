#include "aizalib.h"
/*
 * Euclidean-Like Algorithm (类欧几里得算法)
 *
 * Overview:
 *     用于在对数时间内计算下取整一次函数和 f(a, b, c, n) = sum_{i=0}^n floor((a*i +
 *     b) / c)。
 *     利用辗转相除思想在斜率 > 1 时消去整除商、斜率 < 1
 *     时对调主副坐标轴翻转求和区域，在 O(log(min(a, c, n)))
 *     步内完成二维阶梯网格点计数。
 *
 * API:
 *     sim_euclid::solve(a, b, c, n) — 计算 sum_{i=0}^n floor((a*i + b) / c)
 *
 * Notes:
 *     1. 参数要求: a, b, c >= 0 且 c > 0, n >= 0。
 *     2. 时间复杂度: 复杂度与 gcd(a, c) 相同，为 O(log(min(a, c, n)))。
 *     3. 空间复杂度: 递归栈空间 O(log(min(a, c, n)))。
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
