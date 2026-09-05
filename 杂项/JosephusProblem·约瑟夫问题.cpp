#include "aizalib.h"
/*
 * 约瑟夫问题 (Josephus Problem)
 *
 * Overview:
 *      求解 n 个人（编号 0 到 n - 1）围成一圈，每次报数到 k 的人出列，最后剩下的人的编号。
 *
 * API:
 *      josephus(n, k)     — 经典线性递推，时间复杂度 O(n)。
 *      josephus_log(n, k) — 优化跳跃版本，时间复杂度 O(k log_k n)，适用于 n 极大而 k 较小的情况。
 *      josephus_k2(n)     — k = 2 位运算 O(1) 求解。
 *
 * Notes:
 *      1. [重要] 所有函数返回 0-based 编号，结果域在 [0, n - 1]。若题目要求 1-based，请在外层 +1。
 *      2. 递推公式 f(n, k) = (f(n - 1, k) + k) % n。
 *      3. 要求 n >= 1, k >= 1。
 */

namespace Josephus {

inline int josephus(int n, int k) {
    int res = 0;
    for (int i = 2; i <= n; i++) {
        res = (res + k) % i;
    }
    return res;
}

inline i64 josephus_log(i64 n, i64 k) {
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

inline i64 josephus_k2(i64 n) {
    return (n ^ (i64)std::bit_floor((u64)n)) << 1;
}

} // namespace Josephus

using namespace Josephus;
