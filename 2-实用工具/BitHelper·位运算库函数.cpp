#include "aizalib.h"
/*
 * 位运算库函数助手 (Bit Helper)
 *
 * Overview:
 *     整理与封装现代 C++ (C++20/C++23 <bit>) 竞赛常用位运算指令与算术技巧。
 *
 * API:
 *     popcount_u32(x) / popcount_u64(x)       — 统计二进制中置位 1 的个数。
 *     countr_zero_u32(x) / countr_zero_u64(x) — 统计尾零个数，x = 0
 *                                                时安全返回位宽。
 *     countl_zero_u32(x) / countl_zero_u64(x) — 统计前导零个数，x = 0
 *                                                时安全返回位宽。
 *
 *     lowbit(x)                   — 获取最低位 1 对应的数值权重 x & -x。
 *     clear_lowbit(x)             — 清除最低位 1，返回 x & (x - 1)。
 *     is_pow2(x)                  — 判断正整数 x 是否为 2 的正整数次幂。
 *     highest_bit_index_u64(x)    — 获取最高位 1 的下标（0-based，要求 x != 0）。
 *     bit_width_u64(x)            — 存储非零整数 x 所需的最小二进制位数。
 *     bit_floor_u64(x)            — 不大于 x 的最大 2 的整数次幂。
 *     bit_ceil_u64(x)             — 不小于 x 的最小 2 的整数次幂。
 *     rotl64(x, k) / rotr64(x, k) — 64 位循环左移 / 循环右移。
 *     foreach_subset(mask, f)     — 遍历 mask 的全部子集（包含 0 和全集）。
 *
 * Notes:
 *     1. Time: 单个算术操作均为 O(1)；foreach_subset 为 O(2^{popcount(mask)})。
 *     2. Space: O(1)。
 */

inline int popcount_u32(u32 x) { return std::popcount(x); }
inline int popcount_u64(u64 x) { return std::popcount(x); }

inline int countr_zero_u32(u32 x) { return x ? std::countr_zero(x) : 32; }
inline int countr_zero_u64(u64 x) { return x ? std::countr_zero(x) : 64; }
inline int countl_zero_u32(u32 x) { return x ? std::countl_zero(x) : 32; }
inline int countl_zero_u64(u64 x) { return x ? std::countl_zero(x) : 64; }

template<std::integral T>
inline T lowbit(T x) {
    return x & -x;
}

template<std::integral T>
inline T clear_lowbit(T x) {
    return x & (x - 1);
}

template<std::integral T>
inline bool is_pow2(T x) {
    return x > 0 && ((x & (x - 1)) == 0);
}

inline int highest_bit_index_u64(u64 x) {
    AST(x);
    return std::bit_width(x) - 1;
}

inline int bit_width_u64(u64 x) { return std::bit_width(x); }
inline u64 bit_floor_u64(u64 x) { return std::bit_floor(x); }
inline u64 bit_ceil_u64(u64 x) { return x <= 1 ? 1 : std::bit_ceil(x); }

inline u64 rotl64(u64 x, int k) { return std::rotl(x, k); }
inline u64 rotr64(u64 x, int k) { return std::rotr(x, k); }

template<std::integral T, class F>
inline void foreach_subset(T mask, F&& f) {
    for (T s = mask;; s = (s - 1) & mask) {
        f(s);
        if (!s) break;
    }
}
