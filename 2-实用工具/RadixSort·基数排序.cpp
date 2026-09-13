#include "aizalib.h"
/*
 * 基数排序 (Radix Sort - LSD)
 *
 * Overview:
 *     针对 32 位/64 位整数的高性能最低有效位优先（LSD）非比较基数排序。
 *     - 工具：RadixSort 结构体、sort 重载（vector<int> 与 vector<i64>）。
 *
 * API:
 *     RadixSort() — 初始化基数排序桶内存。
 *     sort(a)     — 对 std::vector<int> 或 std::vector<i64> 原地升序排序。
 *
 * Notes:
 *     1. Time: O(P * (n + 2^16))，P 为有效趟数（int 至多 2 趟，i64 至多 4 趟）。
 *     2. Space: O(n + 2^16)，需大小为 n 的辅助缓冲区与 65536 大小的计数组。
 */

struct RadixSort {
    static constexpr int B = 16;
    static constexpr int M = 1 << B;
    std::vector<int> cnt;

    RadixSort() : cnt(M) {}

    void sort(std::vector<int>& a) {
        int n = a.size();
        if (n <= 1) return;
        constexpr u32 MASK = M - 1;
        constexpr u32 SIGN = 1u << 31;
        u32 diff = 0, base = (u32)a[0] ^ SIGN;
        rep(i, 1, n - 1) diff |= ((u32)a[i] ^ SIGN) ^ base;
        std::vector<int> buf(n);
        auto* src = &a;
        auto* dst = &buf;
        rep(seg, 0, 1) {
            int shift = seg * B;
            if (((diff >> shift) & MASK) == 0) continue;
            std::fill(cnt.begin(), cnt.end(), 0);
            for (int x : *src) ++cnt[((u32)x ^ SIGN) >> shift & MASK];
            rep(i, 1, M - 1) cnt[i] += cnt[i - 1];
            per(i, n - 1, 0) {
                int x = (*src)[i];
                (*dst)[--cnt[((u32)x ^ SIGN) >> shift & MASK]] = x;
            }
            std::swap(src, dst);
        }
        if (src != &a) a.swap(*src);
    }

    void sort(std::vector<i64>& a) {
        int n = a.size();
        if (n <= 1) return;
        constexpr u64 MASK = M - 1;
        constexpr u64 SIGN = 1ull << 63;
        u64 diff = 0, base = (u64)a[0] ^ SIGN;
        rep(i, 1, n - 1) diff |= ((u64)a[i] ^ SIGN) ^ base;
        std::vector<i64> buf(n);
        auto* src = &a;
        auto* dst = &buf;
        rep(seg, 0, 3) {
            int shift = seg * B;
            if (((diff >> shift) & MASK) == 0) continue;
            std::fill(cnt.begin(), cnt.end(), 0);
            for (i64 x : *src) ++cnt[((u64)x ^ SIGN) >> shift & MASK];
            rep(i, 1, M - 1) cnt[i] += cnt[i - 1];
            per(i, n - 1, 0) {
                i64 x = (*src)[i];
                (*dst)[--cnt[((u64)x ^ SIGN) >> shift & MASK]] = x;
            }
            std::swap(src, dst);
        }
        if (src != &a) a.swap(*src);
    }
};
