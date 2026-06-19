#include "aizalib.h"
/**
 * Bit Helper
 * 算法介绍: 整理竞赛里常用的位运算库函数速查，并提供少量可直接复用的轻量包装
 * 模板参数: None
 * Interface:
 * 		popcount_u32/u64(x)						: 统计二进制 1 的个数
 * 		countr_zero_u32/u64(x)					: 统计尾零个数，x = 0 时返回位宽
 * 		countl_zero_u32/u64(x)					: 统计前导零个数，x = 0 时返回位宽
 * 		lowbit(x) / clear_lowbit(x) / is_pow2(x)	: 常见位运算封装
 * 		highest_bit_index_u64(x)				: 最高位下标(0-based)，要求 x != 0
 * 		bit_width_u64(x) / bit_floor_u64(x)
 * 		/ bit_ceil_u64(x)						: C++20 <bit> 常用包装
 * 		rotl64(x, k) / rotr64(x, k)			: 64 位循环移位
 * 		foreach_subset(mask, f)				: 枚举 mask 的全部子集(包含 0)
 * Note:
 * 		1. Time: 各包装均为 O(1)；foreach_subset 为 O(子集数)
 * 		2. Space: O(1)
 * 		3. clz/ctz 对 0 的原生 builtin 未定义，本文件包装已统一处理
 * 		4. 若环境支持 C++20 <bit>，bit_width / bit_floor / bit_ceil / rotl / rotr 可直接联想使用
 */

inline int popcount_u32(u32 x) { return std::popcount(x); }
inline int popcount_u64(u64 x) { return std::popcount(x); }

inline int countr_zero_u32(u32 x) { return x ? std::countr_zero(x) : 32; }
inline int countr_zero_u64(u64 x) { return x ? std::countr_zero(x) : 64; }
inline int countl_zero_u32(u32 x) { return x ? std::countl_zero(x) : 32; }
inline int countl_zero_u64(u64 x) { return x ? std::countl_zero(x) : 64; }

template<typename T>
inline T lowbit(T x) {
	static_assert(std::is_integral_v<T>);
	return x & -x;
}

template<typename T>
inline T clear_lowbit(T x) {
	static_assert(std::is_integral_v<T>);
	return x & (x - 1);
}

template<typename T>
inline bool is_pow2(T x) {
	static_assert(std::is_integral_v<T>);
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

template<class F>
inline void foreach_subset(int mask, F&& f) {
	for (int s = mask;; s = (s - 1) & mask) {
		f(s);
		if (!s) break;
	}
}
