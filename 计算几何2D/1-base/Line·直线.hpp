#pragma once
#include "PointFP·浮点二维点.hpp"

/*
 * 有向直线
 *
 * Overview:
 * 	提供一般式有向直线、构造方法、位置判定及常用浮点直线运算。
 *
 * API:
 * 	Line<T>(a, b, c) / Line<T>(p1, p2): 由一般式 / 两点构造 ax + by + c = 0。
 * 	Line<T>::from_point_and_dir(p, dir): 由一点和方向向量构造。
 * 	Line<T>::from_point_and_normal(p, normal): 由一点和法向量构造。
 * 	eval(p) / side(p): 计算 ax + by + c / 返回点在左侧、右侧、直线上的 1 / -1 / 0。O(1)。
 * 	direction(l) / normal(l): 单位方向向量 (b, -a) / 单位法向量 (a, b)。O(1)。
 * 	parallel(l1, l2) / same_dir(l1, l2): 判断平行 / 平行且同向。O(1)。
 * 	intersection(l1, l2): 两条非平行直线的唯一交点。O(1)。
 * 	distance_to_point(l, p): 点到直线的无符号距离。O(1)。
 * 	projection(l, p) / reflection(l, p): 点在直线上的投影 / 关于直线的对称点。O(1)。
 * 	LineFP: Line<ld> 的常用别名。
 *
 * Notes:
 * 	有向直线左侧（法向量一侧）为正；浮点 side / 平行判定使用 EPS。
 * 	direction、normal 和 distance_to_point 要求 (a, b) != (0, 0)；intersection 要求两线不平行。
 */
namespace Geo2D {

template<typename T>
struct Line {
	T a, b, c;

	Line() : a(0), b(0), c(0) {}
	Line(T a, T b, T c) : a(a), b(b), c(c) {}

	// 两点构造有向直线 (p1 -> p2),左侧为正,右侧为负
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

	T eval(Point<T> p) const {
		return a * p.x + b * p.y + c;
	}

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
