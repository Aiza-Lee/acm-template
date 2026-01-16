#pragma once
#include "../base/Point.hpp"

namespace Geo2D {

struct Segment {
	Point a, b;
	Segment(Point a = Point(), Point b = Point()) : a(a), b(b) {}
};

/**
 * @brief 判断点 p 是否位于线段 s 上（包含端点）。
 * note：
 * 		1. 叉积为 0 (共线)
 * 		2. 点积 <= 0 (点在线段两端点之间)
 */
bool on_segment(Point p, Segment s) {
	return dcmp((s.a - p).cross(s.b - p), 0) == 0 &&
			dcmp((s.a - p).dot(s.b - p), 0) <= 0;
}

/**
 * @brief 点到线段的距离
 */
ld distance_to_segment(Point p, Segment s) {
	if (s.a == s.b) return p.dist_to(s.a);
		
	Point v1 = s.b - s.a, v2 = p - s.a, v3 = p - s.b;
	if (dcmp(v1.dot(v2), 0) < 0) return v2.len();
	if (dcmp(v1.dot(v3), 0) > 0) return v3.len();
		
	return abs(v1.cross(v2)) / v1.len();
}

/**
 * @brief 求两线段的交点
 * 
 * note:
 * 		1. 调用前必须确保两线段确实相交（可先调用 segment_intersect）。
 * 		2. 如果两线段重合（有无穷多交点），该函数行为未定义（通常返回其中一个端点或计算出的某点）。
 */
Point segment_intersection_point(Segment s1, Segment s2) {
	Point a = s1.a, b = s1.b, c = s2.a, d = s2.b;
	ld t = ((c - a).cross(d - c)) / ((b - a).cross(d - c));
	return a + (b - a) * t;
}

/**
 * @brief 判断线段 s1 和 s2 是否相交。
 * 
 * 支持规范相交（交点唯一且不在端点）和非规范相交（端点重合或部分重合）。
 */
bool segment_intersect(Segment s1, Segment s2) {
	Point &a = s1.a, &b = s1.b, &c = s2.a, &d = s2.b;
		
	int d1 = sign((b - a).cross(c - a));
	int d2 = sign((b - a).cross(d - a));
	int d3 = sign((d - c).cross(a - c));
	int d4 = sign((d - c).cross(b - c));
		
	if (d1 * d2 < 0 && d3 * d4 < 0) return true;
	if (d1 == 0 && on_segment(c, s1)) return true;
	if (d2 == 0 && on_segment(d, s1)) return true;
	if (d3 == 0 && on_segment(a, s2)) return true;
	if (d4 == 0 && on_segment(b, s2)) return true;
		
	return false;
}

} // namespace Geo2D