#pragma once
#include "aizalib.h"

namespace Geo2D {

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

} // namespace Geo2D
