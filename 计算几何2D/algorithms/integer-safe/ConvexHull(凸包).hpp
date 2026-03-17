#pragma once
#include "aizalib.h"
#include "../../2-shapes/Polygon.hpp"

/**
 * [ConvexHull (凸包)]
 * 算法介绍: Andrew 算法求二维点集凸包
 * 模板参数: T (点集坐标类型)
 * Interface: 
 *   - Polygon<T> convex_hull(std::vector<Point<T>> pts)
 * Note:
 * 1. Time: O(N log N)
 * 2. Space: O(N)
 * 3. 结果按逆时针顺序排列，自动去重且不包含多余共线点。若需保留共线点，应将 cross 判断改为 <= 0。
 */

namespace Geo2D {

template<typename T>
Polygon<T> convex_hull(std::vector<Point<T>> pts) {
	std::sort(pts.begin(), pts.end());
	pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
	int n = pts.size(), k = 0;
	if (n <= 1) return pts;
	
	Polygon<T> hull(2 * n);
	// 下凸包
	rep(i, 0, n - 1) {
		while (k > 1 && sgn((hull[k-1] - hull[k-2]).cross(pts[i] - hull[k-1])) <= 0) --k;
		hull[k++] = pts[i];
	}
	// 上凸包
	int t = k;
	per(i, n - 2, 0) {
		while (k > t && sgn((hull[k-1] - hull[k-2]).cross(pts[i] - hull[k-1])) <= 0) --k;
		hull[k++] = pts[i];
	}
	hull.resize(k - 1);
	return hull;
}

} // namespace Geo2D
