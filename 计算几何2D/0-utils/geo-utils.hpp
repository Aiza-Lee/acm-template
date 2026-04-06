#pragma once
#include "aizalib.h"

namespace Geo2D {

const ld PI = acos(-1.0);

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

// 保留兼容性
inline int sign(ld x) { return sgn(x); }
inline int dcmp(ld x, ld y) { return cmp(x, y); }

} // namespace Geo2D