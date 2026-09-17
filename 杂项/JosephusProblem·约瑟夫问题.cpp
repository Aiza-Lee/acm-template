#include "aizalib.h"
/*
 * 约瑟夫问题 (Josephus Problem)
 *
 * Overview:
 *     求解 n 个人（编号 0 到 n - 1）围成一圈，每次报数到 k 的人出列，
 *     最后剩下的人的初始编号（幸存者编号）。
 *     - 状态转移结构：设 f(i, k) 为 i 个人时的幸存者 0-based 编号。每次第 k
 *       个人出列后，下一次报数起点发生偏移，产生坐标同余映射 f(i, k) = (f(i - 1, k)
 *       + k) mod i。
 *     - 批量跳跃加速：当 k 较小且 res + k < i 时，连续多次转移均不发生取模溢出，
 *       可通过除法批量跳跃 step = (i - res - 1) / (k - 1) 步，将复杂度降至 O(k
 *       log_k n)。
 *     - 位运算闭式解：当 k = 2 时，幸存者编号对应将 n 的二进制最高位消除并左移一位
 *       (n ^ bit_floor(n)) << 1。
 *     - 工具：线性递推 josephus、对数跳跃递推 josephus_log、常数位运算求解
 *       josephus_k2。
 *
 * API:
 *     josephus(n, k)     — 经典线性递推求解幸存者，时间复杂度 O(n)
 *     josephus_log(n, k) — 批量跳跃优化版本，时间复杂度 O(k log_k n)，适用于 n
 *                           极大而 k 较小
 *     josephus_k2(n)     — k = 2 特化位运算求解，时间复杂度 O(1)
 *
 * Notes:
 *     1. 所有函数统一返回 0-based 编号，结果域为 [0, n - 1]。若需 1-based 请在外层
 *        + 1。
 *     2. 要求输入规模满足 n >= 1, k >= 1。
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
