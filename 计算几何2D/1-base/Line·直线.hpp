#pragma once
#include "PointFP·浮点二维点.hpp"

namespace Geo2D {

/** @brief 有向直线 ax + by + c = 0;(a,b) 是法向量,左侧 (法向一侧) 为正。T 为坐标类型。 */
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

	/** @brief 由一点 + 方向向量构造。 */
	static Line from_point_and_dir(Point<T> p, Point<T> dir) {
		return Line(p, p + dir);
	}

	/** @brief 由一点 + 法向量构造。 */
	static Line from_point_and_normal(Point<T> p, Point<T> normal) {
		return Line(normal.x, normal.y, -(normal.x * p.x + normal.y * p.y));
	}

	/** @brief 代入点坐标: ax + by + c。O(1)。 */
	T eval(Point<T> p) const {
		return a * p.x + b * p.y + c;
	}

	/** @brief 点与有向直线的位置关系;1 左侧, -1 右侧, 0 在直线上(浮点经 sgn 含 EPS)。O(1)。 */
	int side(Point<T> p) const {
		return sgn(eval(p));
	}
};

/** @brief 直线的单位方向向量 (b, -a)。浮点专用,O(1)。退化:直线为点 (a=b=0) 会触发 normalize 的 assert。 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> direction(const Line<T>& l) { return normalize(Point<T>(l.b, -l.a)); }

/** @brief 直线的单位法向量 (a, b)。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> normal(const Line<T>& l) { return normalize(Point<T>(l.a, l.b)); }

/** @brief 两直线是否平行(浮点经 sgn 含 EPS,整数严格)。O(1)。 */
template<typename T>
bool parallel(const Line<T>& l1, const Line<T>& l2) {
	return sgn(l1.a * l2.b - l1.b * l2.a) == 0;
}

/** @brief 两直线是否平行且同向(法向量内积 > 0)。O(1)。 */
template<typename T>
bool same_dir(const Line<T>& l1, const Line<T>& l2) {
	return parallel(l1, l2) && sgn(l1.a * l2.a + l1.b * l2.b) > 0;
}

/** @brief 两直线唯一交点;调用前必须保证不平行 (否则 assert)。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> intersection(const Line<T>& l1, const Line<T>& l2) {
	AST(!parallel(l1, l2));
	T det = l1.a * l2.b - l1.b * l2.a;
	return Point<T>((l1.b * l2.c - l1.c * l2.b) / det, (l1.c * l2.a - l1.a * l2.c) / det);
}

/** @brief 点到直线的距离 (无符号)。浮点专用,O(1)。退化:直线为点 (a=b=0) 会除零返回 NaN/inf。 */
template<typename T>
requires std::is_floating_point_v<T>
T distance_to_point(const Line<T>& l, Point<T> p) {
	return std::abs(l.eval(p)) / std::sqrt(l.a * l.a + l.b * l.b);
}

/** @brief 点在直线上的正交投影。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> projection(const Line<T>& l, Point<T> p) {
	Point<T> n(l.a, l.b);
	return p - n * (l.eval(p) / n.len2());
}

/** @brief 点关于直线的对称像。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> reflection(const Line<T>& l, Point<T> p) {
	return p + (projection(l, p) - p) * 2;
}

/** @brief 常用浮点实例别名。 */
using LineFP = Line<ld>;

} // namespace Geo2D
