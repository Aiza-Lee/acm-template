#pragma once
#include "aizalib.h"

namespace Geo2D {

/** @brief 浮点容差 (全局默认);用作绝对误差阈值。 */
inline constexpr ld EPS = 1e-10;

/** @brief π (long double),来自 std::numbers。 */
inline constexpr ld PI = std::numbers::pi_v<ld>;

/** @brief 符号函数;浮点路径用 EPS 比较,整数路径严格。O(1)。 */
template<typename T>
inline int sgn(T x) {
	if constexpr (std::is_floating_point_v<T>) {
		return x < -EPS ? -1 : x > EPS;
	} else {
		return x < 0 ? -1 : x > 0;
	}
}

/** @brief 三路比较;浮点经 sgn 走 EPS,整数严格比较。O(1)。 */
template<typename T>
inline int cmp(T x, T y) {
	if constexpr (std::is_floating_point_v<T>) {
		return sgn(x - y);
	} else {
		return x < y ? -1 : (x > y ? 1 : 0);
	}
}

/** @brief 是否在零附近;等价于 sgn(x) == 0。O(1)。 */
template<typename T>
inline bool is_zero(T x) {
	return sgn(x) == 0;
}

/** @brief 把 x 钳到 [-1, 1] (浮点专用,用于 acos/asin 入参)。O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
inline T clamp_unit(T x) {
	return std::clamp(x, (T)-1, (T)1);
}

/** @brief 安全的 sqrt;负输入返回 0 避免 NaN。O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
inline T safe_sqrt(T x) {
	return std::sqrt(std::max((T)0, x));
}

/** @brief 安全的 acos;入参先钳到 [-1, 1] 避免定义域外 NaN。O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
inline T safe_acos(T x) {
	return std::acos(clamp_unit(x));
}

/** @brief 安全的 asin;入参先钳到 [-1, 1] 避免定义域外 NaN。O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
inline T safe_asin(T x) {
	return std::asin(clamp_unit(x));
}

} // namespace Geo2D
