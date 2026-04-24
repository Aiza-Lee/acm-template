#include "aizalib.h"

/**
 * Radix Sort
 * 算法介绍: 对 int / i64 做 LSD 基数排序；每趟按 16 bit 计数，通过翻转符号位转成无符号序。
 * Interface:
 * 		RadixSort rs;
 * 		rs.sort(a);		// 将 vector<int> 或 vector<i64> 升序排序
 * Note:
 * 		1. Time: O(P * (n + 2^16))，P 为有效趟数，32 位至多 2 趟，64 位至多 4 趟
 * 		2. Space: O(n + 2^16)
 * 		3. 常数优化: 16 bit 分桶 + 稳定计数排序 + 若某段 16 bit 全相同则自动跳过该趟
 * 		4. 适合整数排序卡常；模板里只保留 int / i64 两套常用版本
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
