#pragma once
#include "Segment.hpp"
#include "../1-base/PointFP.hpp"

namespace Geo2D {

/**
 * @brief 点到线段的距离
 */
template<typename T>
requires std::is_floating_point_v<T>
T distance_to_segment(Point<T> p, Segment<T> s) {
	if (s.a == s.b) return dist_to(p, s.a);
		
	Point<T> v1 = s.b - s.a, v2 = p - s.a, v3 = p - s.b;
	if (sgn(v1.dot(v2)) < 0) return len(v2);
	if (sgn(v1.dot(v3)) > 0) return len(v3);
		
	return std::abs(v1.cross(v2)) / len(v1);
}

/**
 * @brief 求两线段的交点
 * 
 * note:
 * 		1. 调用前必须确保两线段确实相交（可先调用 segment_intersect）。
 * 		2. 如果两线段重合（有无穷多交点），该函数行为未定义（通常返回其中一个端点或计算出的某点）。
 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> segment_intersection_point(Segment<T> s1, Segment<T> s2) {
	Point<T> a = s1.a, b = s1.b, c = s2.a, d = s2.b;
	T t = ((c - a).cross(d - c)) / ((b - a).cross(d - c));
	return a + (b - a) * t;
}

using SegmentFP = Segment<ld>;

} // namespace Geo2D