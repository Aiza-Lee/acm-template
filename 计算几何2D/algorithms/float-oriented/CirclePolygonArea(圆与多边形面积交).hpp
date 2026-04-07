#pragma once
#include "../../2-shapes/Circle.hpp"
#include "../../2-shapes/Polygon.hpp"

/**
 * [CirclePolygonArea (圆与多边形面积交)]
 * 算法介绍: 计算圆与简单多边形（不一定是凸多边形）的重叠面积。使用有向面积三角剖分和积分思想。
 * 模板参数: T (浮点数类型，requires std::is_floating_point_v<T>)
 * Interface: 
 *   - T circle_polygon_area(Circle<T> c, const Polygon<T>& poly)
 * Note:
 * 1. Time: O(N)
 * 2. Space: O(1)
 * 3. 多边形点需按逆时针顺序，如果顺时针求出来的结果需要取绝对值。
 */

namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
T circle_polygon_area(Circle<T> c, const Polygon<T>& poly) {
	auto calc = [&](Point<T> a, Point<T> b) {
		if (sgn((a - c.c).cross(b - c.c)) == 0) return (T)0.0;

		T da = dist_to(c.c, a), db = dist_to(c.c, b);
		if (sgn(da - c.r) <= 0 && sgn(db - c.r) <= 0) {
			return (a - c.c).cross(b - c.c) / 2.0;
		}

		Point<T> pa = a - c.c, pb = b - c.c;
		Point<T> v = pb - pa;
		T A = pa.dist2(pb);
		T B = 2 * pa.dot(v);
		T C = pa.len2() - c.r * c.r;
		T delta = B * B - 4 * A * C;

		T t1 = -1, t2 = -1;
		if (sgn(delta) >= 0) {
			delta = safe_sqrt(delta);
			t1 = (-B - delta) / (2 * A);
			t2 = (-B + delta) / (2 * A);
		}

		std::vector<T> t_vals;
		t_vals.push_back(0);
		t_vals.push_back(1);
		if (sgn(t1) > 0 && sgn(t1 - 1) < 0) t_vals.push_back(t1);
		if (sgn(t2) > 0 && sgn(t2 - 1) < 0) t_vals.push_back(t2);
		std::sort(t_vals.begin(), t_vals.end());

		T res = 0;
		for (size_t i = 0; i < t_vals.size() - 1; i++) {
			T mid_t = (t_vals[i] + t_vals[i+1]) / 2;
			Point<T> mid = pa + v * mid_t;
			Point<T> p1 = pa + v * t_vals[i];
			Point<T> p2 = pa + v * t_vals[i+1];

			if (sgn(len(mid) - c.r) <= 0) {
				res += p1.cross(p2) / 2.0;
			} else {
				T angle_val = std::atan2(p1.cross(p2), p1.dot(p2));
				res += c.r * c.r * angle_val / 2.0;
			}
		}
		return res;
	};

	T area_val = 0;
	int n = poly.size();
	rep(i, 0, n - 1) {
		area_val += calc(poly[i], poly[(i + 1) % n]);
	}
	return area_val;
}

} // namespace Geo2D
