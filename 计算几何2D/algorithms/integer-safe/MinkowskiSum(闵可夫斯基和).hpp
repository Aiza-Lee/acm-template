#pragma once
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
	if (!n || !m) return {};
	auto shift = [](const Polygon<T>& poly, Point<T> d) {
		Polygon<T> res = poly;
		for (auto& p : res) p = p + d;
		return res;
	};
	if (n == 1) return shift(Q, P[0]);
	if (m == 1) return shift(P, Q[0]);
	auto reorder = [](Polygon<T> poly) {
		int pos = 0;
		rep(i, 1, (int)poly.size() - 1) {
			if (cmp(poly[i].y, poly[pos].y) < 0 || (cmp(poly[i].y, poly[pos].y) == 0 && cmp(poly[i].x, poly[pos].x) < 0))
				pos = i;
		}
		std::rotate(poly.begin(), poly.begin() + pos, poly.end());
		return poly;
	};
	auto normalize = [](Polygon<T> poly) {
		if (!poly.empty() && poly.front() == poly.back()) poly.pop_back();
		if (poly.size() <= 2) return poly;
		Polygon<T> res;
		int sz = poly.size();
		rep(i, 0, sz - 1) {
			Point<T> pre = poly[(i - 1 + sz) % sz];
			Point<T> cur = poly[i];
			Point<T> nxt = poly[(i + 1) % sz];
			if (sgn((cur - pre).cross(nxt - cur)) == 0) continue;
			res.emplace_back(cur);
		}
		return res.empty() ? Polygon<T>{poly[0]} : res;
	};
	Polygon<T> p = reorder(P), q = reorder(Q);

	std::vector<Point<T>> e1(n), e2(m);
	rep(i, 0, n - 1) e1[i] = p[(i + 1) % n] - p[i];
	rep(i, 0, m - 1) e2[i] = q[(i + 1) % m] - q[i];

	Polygon<T> res;
	res.reserve(n + m + 1);
	Point<T> cur = p[0] + q[0];
	res.emplace_back(cur);
	int i = 0, j = 0;
	while (i < n || j < m) {
		if (j == m || (i < n && sgn(e1[i].cross(e2[j])) > 0)) {
			cur = cur + e1[i++];
		} else if (i == n || sgn(e1[i].cross(e2[j])) < 0) {
			cur = cur + e2[j++];
		} else {
			cur = cur + e1[i++] + e2[j++];
		}
		if (i < n || j < m) res.emplace_back(cur);
	}
	return normalize(res);
}

} // namespace Geo2D
