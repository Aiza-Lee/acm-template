#pragma once
#include "../../2-shapes/Polygon·多边形.hpp"

/**
 * [ConvexHull (凸包)]
 * 算法介绍: Andrew 单调链算法求二维点集凸包,以及凸包的几何重心。
 * 模板参数: T (点集/凸包顶点坐标类型)
 * Interface:
 *   - Polygon<T> convex_hull(std::vector<Point<T>> pts)
 *   - Point<ld>  convex_hull_centroid(const Polygon<T>& hull)
 * Note:
 * 1. convex_hull:
 *    a. Time: O(N log N); Space: O(N).
 *    b. 结果按逆时针顺序排列,自动去重且不包含多余共线点。若需保留共线点,应将 cross 判断改为 < 0 (而非 <= 0)。
 * 2. convex_hull_centroid:
 *    a. Time: O(N).
 *    b. 输入必须是逆时针凸多边形 (顶点集为 convex_hull 的返回值);退化情形由 polygon_centroid 处理 (area=0 时返回首顶点)。
 *    c. 返回 Point<ld>,无论输入 T 是整数还是浮点数 (公式含除法)。
 *    d. 底层转调 Polygon·多边形.hpp::polygon_centroid,与多边形重心共享同一实现。
 *    @see Polygon·多边形.hpp::polygon_centroid — 多边形层通用质心;此函数为空凸包兜底 (0,0)。
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

template<typename T>
Point<ld> convex_hull_centroid(const Polygon<T>& hull) {
	if (hull.empty()) return Point<ld>(0, 0);   // 本地兜底:空凸包返回原点
	// 转: Polygon·多边形.hpp::polygon_centroid
	return polygon_centroid(hull);
}

} // namespace Geo2D
