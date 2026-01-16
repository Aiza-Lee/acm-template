#pragma once
#include "../shapes/Polygon.hpp"

namespace Geo2D {
/**
 * @brief 求凸包 (Andrew 算法)
 * 
 * 给定点集，求其凸包。
 * 结果按逆时针顺序排列。
 * 
 * @param pts 输入点集
 * @return Polygon 凸包点集
 * @note 自动去重，结果不包含共线点（除非只有共线点）。
 * @complexity O(N log N) (排序瓶颈)
 */
Polygon convex_hull(std::vector<Point> pts) {
	sort(pts.begin(), pts.end());
	pts.erase(unique(pts.begin(), pts.end()), pts.end());
	int n = pts.size(), k = 0;
	if (n <= 1) return pts;
	Polygon hull(2 * n);
	// 下凸包
	for(int i = 0; i < n; i++) {
		while (k > 1 && sign((hull[k-1] - hull[k-2]).cross(pts[i] - hull[k-1])) < 0) --k;
		hull[k++] = pts[i];
	}
	// 上凸包
	int t = k;
	for(int i = n - 2; i >= 0; i--) {
		while (k > t && sign((hull[k-1] - hull[k-2]).cross(pts[i] - hull[k-1])) < 0) --k;
		hull[k++] = pts[i];
	}
	hull.resize(k - 1); // 若不保留最后一个重复点
	return hull;
}
} // namespace Geo2D
