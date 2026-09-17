#include "aizalib.h"
/*
 * 汉明权重枚举 (Hamming Weight / Gosper's Hack)
 *
 * Overview:
 *     利用 Gosper's Hack 算法按字典序生成包含固定数量置位位（即汉明权重固定为
 *     k）的二进制掩码序列，进而支持按汉明权重递增顺序遍历全集状态。
 *     - Gosper's Hack 转移结构：给定当前掩码 x，求严格大于 x
 *       且置位数量相同的最小掩码 next_x：
 *       1. 取最低置位：c = x & -x。
 *       2. 连续进位：r = x + c。
 *       3. 尾部连续 1 归位至最低位：next_x = (((r ^ x) >> 2) / c) | r。
 *     - 工具：next_combination、for_each_k_subset、
 *       for_each_by_hamming_weight。
 *
 * API:
 *     next_combination(x)              — 给定权重为 k 的掩码 x，求相同权重的字典序
 *                                         后继掩码
 *     for_each_k_subset(n, k, f)       — 枚举 n 位中所有权重为 k 的掩码，
 *                                         依次调用回调 f(mask)
 *     for_each_by_hamming_weight(n, f) — 按汉明权重递增顺序枚举全集状态，
 *                                         依次调用回调 f(mask)
 *
 * Notes:
 *     1. Time: 枚举固定 k 个置位的所有状态复杂度为 O(C(n, k))，无多余遍历。
 *     2. Space: O(1)。
 *     3. 当 x 到达上限组合时，next_combination 会产生进位，调用方以 1 << n
 *        作为终止上界。
 */

inline u64 next_combination(u64 x) {
    u64 c = x & -x;
    u64 r = x + c;
    return (((r ^ x) >> 2) / c) | r;
}

template<typename F>
inline void for_each_k_subset(int n, int k, F&& f) {
    if (k == 0) {
        f(0ULL);
        return;
    }
    if (k > n) return;
    u64 limit = 1ULL << n;
    u64 x = (1ULL << k) - 1;
    while (x < limit) {
        f(x);
        x = next_combination(x);
    }
}

template<typename F>
inline void for_each_by_hamming_weight(int n, F&& f) {
    rep(k, 0, n) {
        for_each_k_subset(n, k, f);
    }
}