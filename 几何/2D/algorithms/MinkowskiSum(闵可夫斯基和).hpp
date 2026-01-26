#pragma once
#include "../2-shapes/Polygon.hpp"
namespace Geo2D {

/**
 * @brief 闵可夫斯基和 (Minkowski Sum)
 * 
 * 计算两个凸多边形 P 和 Q 的闵可夫斯基和。
 * 结果也是一个凸多边形。
 * 
 * @param P 凸多边形 P (需按逆时针排序，且第一个点为 lowest-left) - 代码内部会处理排序
 * @param Q 凸多边形 Q (需按逆时针排序)
 * @return Polygon 闵可夫斯基和的凸包
 * @note 内部会对输入多边形进行重新定序以确保从最低点开始。
 * @complexity O(N + M)
 */
Polygon minkowski_sum(const Polygon& P, const Polygon& Q) {
	int n = P.size(), m = Q.size();
	auto reorder = [](Polygon& poly) {
		int pos = 0;
		for(int i = 1; i < (int)poly.size(); i++) {
			if(poly[i].y < poly[pos].y || (poly[i].y == poly[pos].y && poly[i].x < poly[pos].x))
				pos = i;
		}
		std::rotate(poly.begin(), poly.begin() + pos, poly.end());
	};
	Polygon p = P, q = Q;
	reorder(p); reorder(q);
	
	std::vector<Point> e1(n), e2(m);
	for(int i = 0; i < n; i++) e1[i] = p[(i + 1) % n] - p[i];
	for(int i = 0; i < m; i++) e2[i] = q[(i + 1) % m] - q[i];
	
	std::vector<Point> edge(n + m);
	std::merge(e1.begin(), e1.end(), e2.begin(), e2.end(), edge.begin(), [](const Point& a, const Point& b) {
		return a.cross(b) > 0;
	});
	
	Polygon res;
	res.reserve(n + m + 1);
	res.push_back(p[0] + q[0]);
	for(const auto& v : edge) res.push_back(res.back() + v);
	res.pop_back();
	return res;
}

} // namespace Geo2D
