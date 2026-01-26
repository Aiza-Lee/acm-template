#pragma once
#include "../2-shapes/Circle.hpp"
#include "../2-shapes/Polygon.hpp"

namespace Geo2D {

/**
 * @brief 圆与多边形面积交
 * 
 * 计算圆与简单多边形（不一定是凸的）的重叠面积。
 * 
 * @param c 圆
 * @param poly 多边形
 * @return ld 重叠面积
 * @note 使用三角剖分和积分思想。
 * @complexity O(N)
 */
ld area_circle_polygon(Circle c, Polygon poly) {
	auto calc = [&](Point a, Point b) {
		if (sign((a - c.c).cross(b - c.c)) == 0) return 0.0L;
		
		ld da = c.c.dist_to(a), db = c.c.dist_to(b);
		if (sign(da - c.r) <= 0 && sign(db - c.r) <= 0) {
			return (a - c.c).cross(b - c.c) / 2.0;
		}
		
		Point pa = a - c.c, pb = b - c.c;
		// 检查线段是否与圆相交
		// 使用 Basic.hpp 中的 distance_to_segment 可能会比较麻烦，因为我们需要交点
		// 这里直接解方程
		
		// 转换坐标系，圆心为原点
		// 线段 pa -> pb
		// |p|^2 = r^2, p = pa + (pb-pa)*t, 0<=t<=1
		
		Point v = pb - pa;
		ld A = v.len2();
		ld B = 2 * pa.dot(v);
		ld C = pa.len2() - c.r * c.r;
		ld delta = B * B - 4 * A * C;
		
		ld t1 = -1, t2 = -1;
		if (sign(delta) >= 0) {
			delta = sqrt(std::max(0.0L, delta));
			t1 = (-B - delta) / (2 * A);
			t2 = (-B + delta) / (2 * A);
		}
		
		std::vector<ld> t_vals;
		t_vals.push_back(0);
		t_vals.push_back(1);
		if (sign(t1) > 0 && sign(t1 - 1) < 0) t_vals.push_back(t1);
		if (sign(t2) > 0 && sign(t2 - 1) < 0) t_vals.push_back(t2);
		std::sort(t_vals.begin(), t_vals.end());
		
		ld res = 0;
		for (size_t i = 0; i < t_vals.size() - 1; i++) {
			ld mid_t = (t_vals[i] + t_vals[i+1]) / 2;
			Point mid = pa + v * mid_t;
			Point p1 = pa + v * t_vals[i];
			Point p2 = pa + v * t_vals[i+1];
			
			if (sign(mid.len() - c.r) <= 0) {
				res += p1.cross(p2) / 2.0;
			} else {
				ld angle = atan2(p1.cross(p2), p1.dot(p2));
				res += c.r * c.r * angle / 2.0;
			}
		}
		return res;
	};
	
	ld area = 0;
	int n = poly.size();
	rep(i, 0, n - 1) {
		area += calc(poly[i], poly[(i + 1) % n]);
	}
	return area;
}

} // namespace Geo2D
