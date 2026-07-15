#pragma once
#include "../1-base/Point·点.hpp"

namespace Geo2D {

/** @brief 线段 (有向/无向皆可,由调用方约定);T 为坐标类型。 */
template<typename T>
struct Segment {
	Point<T> a, b;
	Segment() : a(Point<T>()), b(Point<T>()) {}
	Segment(Point<T> a, Point<T> b) : a(a), b(b) {}
};

/**
 * @brief 判断点 p 是否在线段 s 上 (含端点)。
 * @complexity O(1);浮点路径经 sgn 含 EPS。
 */
template<typename T>
bool on_segment(Point<T> p, Segment<T> s) {
	return sgn((s.a - p).cross(s.b - p)) == 0 &&
			sgn((s.a - p).dot(s.b - p)) <= 0;
}

/**
 * @brief 判断两线段是否相交 (含端点接触与共线部分重叠)。
 * @complexity O(1);浮点路径经 sgn 含 EPS。
 */
template<typename T>
bool segment_intersect(Segment<T> s1, Segment<T> s2) {
	Point<T> &a = s1.a, &b = s1.b, &c = s2.a, &d = s2.b;
	
	// 1. 快速排斥实验 (Bounding Box)
	if (std::max(a.x, b.x) < std::min(c.x, d.x) || std::max(c.x, d.x) < std::min(a.x, b.x) ||
	    std::max(a.y, b.y) < std::min(c.y, d.y) || std::max(c.y, d.y) < std::min(a.y, b.y)) {
		return false;
	}

	int d1 = sgn((b - a).cross(c - a));
	int d2 = sgn((b - a).cross(d - a));
	int d3 = sgn((d - c).cross(a - c));
	int d4 = sgn((d - c).cross(b - c));
		
	if (d1 * d2 < 0 && d3 * d4 < 0) return true;
	if (d1 == 0 && on_segment(c, s1)) return true;
	if (d2 == 0 && on_segment(d, s1)) return true;
	if (d3 == 0 && on_segment(a, s2)) return true;
	if (d4 == 0 && on_segment(b, s2)) return true;
		
	return false;
}

} // namespace Geo2D