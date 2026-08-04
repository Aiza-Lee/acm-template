#pragma once
#include "aizalib.h"

/*
 * 几何工具（三维）
 *
 * Overview:
 * 	提供三维计算几何通用常量、带容差比较和安全的浮点初等函数。
 * 	与计算几何2D 版本行为完全等价,仅命名空间不同,避免跨段依赖。
 *
 * API:
 * 	EPS / PI: 全局浮点容差 / long double 圆周率。
 * 	sgn(x) / cmp(x, y) / is_zero(x): 符号、三路比较和判零。O(1)。
 * 	clamp_unit(x): 将浮点数钳到 [-1, 1]。O(1)。
 * 	safe_sqrt(x): 负输入按 0 处理后开方。O(1)。
 * 	safe_acos(x) / safe_asin(x): 钳制输入后计算反三角函数。O(1)。
 *
 * Notes:
 * 	浮点比较使用绝对误差 EPS；整数比较保持严格。
 * 	本文件刻意与 2D 版本独立维护,修改时请同步两边。
 *
 * Related:
 * 	../计算几何2D/0-utils/GeoUtils·几何工具.hpp: 同名等价文件,namespace Geo2D。
 */
namespace Geo3D {

inline constexpr ld EPS = 1e-10;

inline constexpr ld PI = std::numbers::pi_v<ld>;

template<typename T>
inline int sgn(T x) {
	if constexpr (std::is_floating_point_v<T>) {
		return x < -EPS ? -1 : x > EPS;
	} else {
		return x < 0 ? -1 : x > 0;
	}
}

template<typename T>
inline int cmp(T x, T y) {
	if constexpr (std::is_floating_point_v<T>) {
		return sgn(x - y);
	} else {
		return x < y ? -1 : (x > y ? 1 : 0);
	}
}

template<typename T>
inline bool is_zero(T x) {
	return sgn(x) == 0;
}

template<typename T>
requires std::is_floating_point_v<T>
inline T clamp_unit(T x) {
	return std::clamp(x, (T)-1, (T)1);
}

template<typename T>
requires std::is_floating_point_v<T>
inline T safe_sqrt(T x) {
	return std::sqrt(std::max((T)0, x));
}

template<typename T>
requires std::is_floating_point_v<T>
inline T safe_acos(T x) {
	return std::acos(clamp_unit(x));
}

template<typename T>
requires std::is_floating_point_v<T>
inline T safe_asin(T x) {
	return std::asin(clamp_unit(x));
}

} // namespace Geo3D
