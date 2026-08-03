#pragma once
#include "../1-base/Point·点.hpp"

/*
 * 线段
 *
 * Overview:
 * 	提供二维线段类型、点在线段上的判定以及两线段相交判定。
 *
 * API:
 * 	Segment<T>(a, b): 由两端点构造。
 * 	on_segment(p, s): 点 p 是否在线段 s 上（含端点）。O(1)。
 * 	segment_intersect(s1, s2): 两线段是否相交（含端点接触与共线部分重叠）。O(1)。
 *
 * Notes:
 * 	浮点路径使用 sgn 和 EPS 进行容差比较。
 */
namespace Geo2D {

template<typename T>
struct Segment {
	Point<T> a, b;
	Segment() : a(Point<T>()), b(Point<T>()) {}
	Segment(Point<T> a, Point<T> b) : a(a), b(b) {}
};

template<typename T>
bool on_segment(Point<T> p, Segment<T> s) {
	return sgn((s.a - p).cross(s.b - p)) == 0 &&
			sgn((s.a - p).dot(s.b - p)) <= 0;
}

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
