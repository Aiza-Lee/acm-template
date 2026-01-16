#include "aizalib.h"
/*
  cpp_lib_helper.cpp (浮点与位运算速查)

  内容包括：
	- <cmath> 常用函数（指数/对数/三角/数值稳定函数/舍入/分类等）
	- 浮点比较与数值稳定技巧（log1p/expm1/hypot/fma/nextafter）
	- 浮点分类与特殊值（isnan/isinf/signbit）
	- 位运算 GCC/Clang 内建与 C++20 <bit> 的常用封装（popcount/clz/ctz/rot/lowbit/is_pow2 等）
*/

// ------------------- 浮点数常用函数（来自 <cmath> / std） -------------------
// 指数/对数
// std::exp(x), std::exp2(x), std::expm1(x)  -- e^x, 2^x, e^x - 1（小 x 更精确）
// std::log(x), std::log2(x), std::log10(x), std::log1p(x) -- ln, log2, log10, ln(1+x)
// 幂与根： std::pow(x,y), std::sqrt(x), std::cbrt(x)

// 三角与反三角
// 三角与反三角（以弧度为单位）
//
// std::sin(x) / std::cos(x) / std::tan(x)
//   - 定义域：sin, cos 在所有实数上定义；tan 在 x = π/2 + kπ 处无定义。
//   - 值域：sin, cos ∈ [-1, 1]；tan 在其定义点上可取任意实数。
//   - 周期：sin、cos 的周期为 2π；tan 的周期为 π。
//   - 性质：cos(x) = sin(x + π/2)。sin 为奇函数（sin(-x) = -sin x），cos 为偶函数（cos(-x) = cos x）。
//   - 数值注意：对于 |x| 非常大时，直接调用三角函数可能导致精度损失；可先用 std::remainder 或 fmod 将 x 规约到 [-π, π] 或 [-π/2, π/2] 区间以提高精度。
//
// std::asin(x) / std::acos(x) / std::atan(x)
//   - 主值范围与定义域：
//       * asin(x): 定义域 x ∈ [-1, 1]，返回值域 [-π/2, π/2]。
//       * acos(x): 定义域 x ∈ [-1, 1]，返回值域 [0, π]。
//       * atan(x): 定义域 x ∈ R，返回值域 (-π/2, π/2)。
//   - 使用建议：反三角函数用于由三角比恢复角度。因为浮点误差，传入 asin/acos 的值可能略超出 [-1,1]，应先用 std::clamp(v, -1.0, 1.0) 做截断以避免 NaN。
//
// std::atan2(y, x)
//   - 说明：返回点 (x, y) 对应的极角 θ，θ 的范围通常为 (-π, π]，能够区分四个象限。
//   - 推荐用法：始终使用 atan2(y, x) 而不是 atan(y/x)，以避免 x = 0 的除零并正确恢复象限信息。
//
// 双曲函数（hyperbolic）：std::sinh(x), std::cosh(x), std::tanh(x)
//   - 定义：sinh(x) = (e^x - e^{-x})/2，cosh(x) = (e^x + e^{-x})/2，tanh(x) = sinh(x)/cosh(x)。
//   - 性质：sinh 为奇函数（值域 R），cosh 为偶函数（值域 [1, +∞)），tanh 值域 (-1, 1)。
//   - 反双曲：asinh, acosh, atanh；注意 acosh 的定义域为 [1, +∞)，atanh 的定义域为 (-1,1)。
//   - 数值注意：当 |x| 很大时，直接计算 sinh/cosh 可能溢出；可用 expm1 或按对数形式处理以改善稳定性。

// 数值稳定函数
// std::hypot(x,y,...)  -- 稳定地计算 sqrt(x^2 + y^2 + ...)
// std::fma(a,b,c)     -- fused multiply-add，计算 a*b + c，减少中间舍入
// std::expm1(x), std::log1p(x) -- 当 x 接近 0 时更精确
// std::nextafter(x,y) -- x 朝 y 的下一个可表示浮点数

// 舍入/截断
// std::ceil, std::floor, std::trunc, std::round
// std::nearbyint, std::rint, std::llround, std::lround, std::lrint

// 分解/组合
// std::frexp(x, &exp)  -- x = frac * 2^exp，frac in [0.5,1)
// std::ldexp(frac, exp) -- 与 frexp 逆运算
// std::modf(x, &intpart) -- 拆分整数和小数部分

// 分类/判定
// std::isnan, std::isinf, std::isfinite, std::signbit, std::fpclassify

// 符号函数
inline int sgn(double x) { if (x > EPS) return 1; if (x < -EPS) return -1; return 0; }
// 浮点数比较（默认误差 EPS，可调整）
inline bool feq(double a, double b, double eps = EPS) { return std::fabs(a - b) <= eps; }

// 常用示例：
inline void float_examples() {
	double a = std::exp(1.5);
	double b = std::log1p(1e-12);        // 小量时更精确
	double d = std::hypot(3.0, 4.0);    // 5.0
	double y = std::fma(1000000000000000000.0, 1e-6, 1.0); // 避免中间舍入，计算 a * b + c

	double x = 0.1 + 0.2;
	bool eq = feq(x, 0.3);

	double frac; double ipart;
	frac = std::modf(3.14, &ipart); // frac = 0.14, ipart = 3.0

	// nextafter 用法：向目标值迈进一个可表示的浮点单位
	double nx = std::nextafter(1.0, 2.0);

	(void)a; (void)b; (void)d; (void)y; (void)eq; (void)frac; (void)ipart; (void)nx;
}

// ------------------- 位运算（GCC/Clang 内建 与 C++20 <bit>） -------------------
// 注意：__builtin_clz/ctz 对 0 未定义，使用前检查 x != 0

// popcount
inline int popcount_u32(u32 x) { return __builtin_popcount(x); }
inline int popcount_u64(u64 x) { return __builtin_popcountll(x); }

// 计数尾零 / 前导零
inline int countr_zero_u32(u32 x) { return x ? __builtin_ctz(x) : 32; }
inline int countr_zero_u64(u64 x) { return x ? __builtin_ctzll(x) : 64; }
inline int countl_zero_u32(u32 x) { return x ? __builtin_clz(x) : 32; }
inline int countl_zero_u64(u64 x) { return x ? __builtin_clzll(x) : 64; }

// lowbit / clear lowbit / is_pow2
template<typename T>
inline T lowbit(T x) { return x & -x; }
template<typename T>
inline T clear_lowbit(T x) { return x & (x - 1); }
template<typename T>
inline bool is_pow2(T x) { return x && ((x & (x - 1)) == 0); }

// highest bit index (0-based), x must be != 0
inline int highest_bit_index_u64(u64 x) { return 63 - __builtin_clzll(x); }

// rotate left / right （对 64-bit 的通用实现）
inline u64 rotl64(u64 x, int k) {
	const int B = 64;
	k &= (B - 1);
	if (k == 0) return x;
	return (x << k) | (x >> (B - k));
}
inline u64 rotr64(u64 x, int k) {
	const int B = 64;
	k &= (B - 1);
	if (k == 0) return x;
	return (x >> k) | (x << (B - k));
}

// 示例：按位子集遍历
inline void iterate_subsets(int mask) {
	for (int s = mask; s; s = s & (s - 1)) {
		int lb = s & -s; // 子集中最低位
		(void)lb;
	}
}

// ------------------- 注意事项 -------------------
// - 使用 std::hypot 计算距离以避免溢出/下溢
// - 使用 std::log1p / std::expm1 处理小量
// - 使用 std::fma 在需要更好数值稳定性与性能时替换 a*b + c
// - 对 clz/ctz 使用前请检查 x != 0，或者用包装函数处理 0
// - 若使用 C++20 的 <bit> 可替换为 std::popcount/std::countr_zero 等，更可移植