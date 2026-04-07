#pragma once
#include "PointFP.hpp"

namespace Geo2D {

template<typename T>
struct Line {
	T a, b, c;
	
	Line() : a(0), b(0), c(0) {}
	Line(T a, T b, T c) : a(a), b(b), c(c) {}

	// 两点构造有向直线 (p1 -> p2), 左侧为正, 右侧为负
	Line(Point<T> p1, Point<T> p2) {
		a = p1.y - p2.y;
		b = p2.x - p1.x;
		c = -a * p1.x - b * p1.y;
	}
		
	static Line from_point_and_dir(Point<T> p, Point<T> dir) {
		return Line(p, p + dir);
	}
	
	static Line from_point_and_normal(Point<T> p, Point<T> normal) {
		return Line(normal.x, normal.y, -(normal.x * p.x + normal.y * p.y));
	}
	
	// 代入点坐标计算值: ax + by + c
	T eval(Point<T> p) const {
		return a * p.x + b * p.y + c;
	}

	// 1: 左侧, -1: 右侧, 0: 直线上
	int side(Point<T> p) const {
		return sgn(eval(p));
	}
};

template<typename T>
requires std::is_floating_point_v<T>
Point<T> direction(const Line<T>& l) { return normalize(Point<T>(l.b, -l.a)); }

template<typename T>
requires std::is_floating_point_v<T>
Point<T> normal(const Line<T>& l) { return normalize(Point<T>(l.a, l.b)); }

template<typename T>
bool parallel(const Line<T>& l1, const Line<T>& l2) {
	return sgn(l1.a * l2.b - l1.b * l2.a) == 0;
}

template<typename T>
bool same_dir(const Line<T>& l1, const Line<T>& l2) {
	return parallel(l1, l2) && sgn(l1.a * l2.a + l1.b * l2.b) > 0;
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> intersection(const Line<T>& l1, const Line<T>& l2) {
	AST(!parallel(l1, l2));
	T det = l1.a * l2.b - l1.b * l2.a;
	return Point<T>((l1.b * l2.c - l1.c * l2.b) / det, (l1.c * l2.a - l1.a * l2.c) / det);
}

template<typename T>
requires std::is_floating_point_v<T>
T distance_to_point(const Line<T>& l, Point<T> p) {
	return std::abs(l.eval(p)) / std::sqrt(l.a * l.a + l.b * l.b);
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> projection(const Line<T>& l, Point<T> p) {
	Point<T> n(l.a, l.b);
	return p - n * (l.eval(p) / n.len2());
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> reflection(const Line<T>& l, Point<T> p) {
	return p + (projection(l, p) - p) * 2;
}

using LineFP = Line<ld>;

} // namespace Geo2D
