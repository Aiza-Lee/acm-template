#include "aizalib.h"
/**
 * 约瑟夫问题 (Josephus Problem)
 * 算法介绍: 
 *      求解 n 个人（编号 0 到 n-1）围成一圈，每次报数到 k 的人出列，最后剩下的人的编号。
 * 模板参数: 无
 * Interface:
 *     josephus(n, k)           经典线性递推 O(n)
 *     josephus_log(n, k)       优化跳跃版本 O(k log_k n)，适用于 n 极大而 k 较小的情况。
 *     josephus_k2(n)           k = 2 位运算解
 * Note:
 *      1. Time: O(n) 或 O(k log n)
 *      2. Space: O(1)
 *      3. [重要] 所有函数返回的编号是 0-based 索引，即结果域在 [0, n-1]。若题目为 1-based，请在外层结果进行 +1 操作。
 *      4. 递推公式 f(n,k) = (f(n-1,k) + k) % n
 */

namespace Josephus {

int josephus(int n, int k) {
    int res = 0;
    for (int i = 2; i <= n; i++) {
        res = (res + k) % i;
    }
    return res;
}

i64 josephus_log(i64 n, i64 k) {
    if (k == 1) return n - 1;

    i64 res = 0;
    i64 i = 1;

    while (i < n) {
        if (res + k < i) {
            i64 step = (i - res - 1) / (k - 1);
            if (i + step > n) step = n - i;

            res += k * step;
            i += step;
        } else {
            res = (res + k) % (i + 1);
            ++i;
        }
    }
    return res;
}

i64 josephus_k2(i64 n) {
    return ~std::bit_floor(((u64)n) << 1) & ((n << 1) | 1);
}

} // namespace Josephus
