#pragma once
#include "../../2-shapes/Polygon·多边形.hpp"

/**
 * [ConvexHull (凸包)]
 * 算法介绍: Andrew 算法求二维点集凸包
 * 模板参数: T (点集坐标类型)
 * Interface: 
 *   - Polygon<T> convex_hull(std::vector<Point<T>> pts)
 * Note:
 * 1. Time: O(N log N)
 * 2. Space: O(N)
 * 3. 结果按逆时针顺序排列，自动去重且不包含多余共线点。若需保留共线点，应将 cross 判断改为 < 0 (而非 <= 0)。
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

/**
 * [ConvexHullCentroid (凸包重心)]
 * 算法介绍: 凸多边形(凸包)的几何重心(面积加权质心 / 形心),由 shoelace 公式 O(N) 求得
 * 模板参数: T (凸包顶点坐标类型)
 * Interface:
 *   - Point<ld> convex_hull_centroid(const Polygon<T>& hull)
 * Note:
 * 1. Time: O(N)
 * 2. 输入必须是逆时针凸多边形 (顶点集为 convex_hull 的返回值);退化情形由 polygon_centroid 处理 (area=0 时返回首顶点)。
 * 3. 返回 Point<ld>,无论输入 T 是整数还是浮点数 (公式含除法)。
 * 4. 底层转调 Polygon::polygon_centroid,与多边形重心共享同一实现。
 */

template<typename T>
Point<ld> convex_hull_centroid(const Polygon<T>& hull) {
	if (hull.empty()) return Point<ld>(0, 0);
	return polygon_centroid(hull);
}

} // namespace Geo2D
