#pragma once
#include "aizalib.h"
#include "../../2-shapes/Polygon.hpp"

/**
 * [MinkowskiSum (闵可夫斯基和)]
 * 算法介绍: 计算两个凸多边形 P 和 Q 的闵可夫斯基和 (位积)。
 * 模板参数: T (坐标类型)
 * Interface: 
 *   - Polygon<T> minkowski_sum(const Polygon<T>& P, const Polygon<T>& Q)
 * Note:
 * 1. Time: O(N + M)
 * 2. Space: O(N + M)
 * 3. 输入的两个多边形必须是凸多边形并且点集按逆时针有序。
 */

namespace Geo2D {

template<typename T>
Polygon<T> minkowski_sum(const Polygon<T>& P, const Polygon<T>& Q) {
	int n = P.size(), m = Q.size();
	auto reorder = [](Polygon<T>& poly) {
		int pos = 0;
		for(int i = 1; i < (int)poly.size(); i++) {
			if(cmp(poly[i].y, poly[pos].y) < 0 || (cmp(poly[i].y, poly[pos].y) == 0 && cmp(poly[i].x, poly[pos].x) < 0))
				pos = i;
		}
		std::rotate(poly.begin(), poly.begin() + pos, poly.end());
	};
	Polygon<T> p = P, q = Q;
	reorder(p); reorder(q);

	std::vector<Point<T>> e1(n), e2(m);
	rep(i, 0, n - 1) e1[i] = p[(i + 1) % n] - p[i];
	rep(i, 0, m - 1) e2[i] = q[(i + 1) % m] - q[i];

	std::vector<Point<T>> edge(n + m);
	std::merge(
		e1.begin(), e1.end(), 
		e2.begin(), e2.end(), 
		edge.begin(), 
		[](const Point<T>& a, const Point<T>& b) {
			return sgn(a.cross(b)) > 0;
		}
	);

	Polygon<T> res;
	res.reserve(n + m + 1);
	res.push_back(p[0] + q[0]);
	for(const auto& v : edge) res.push_back(res.back() + v);
	res.pop_back();
	return res;
}

} // namespace Geo2D
