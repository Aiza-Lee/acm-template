#pragma once
#include "../1-base/Line·直线.hpp"

namespace Geo2D {

/**
 * @brief 三角形外心 — 外接圆圆心,到三顶点距离相等。
 * @complexity O(1)
 * @note 退化:三顶点共线返回 (NaN, NaN);若三角形退化为线段亦同。
 * @see MinEnclosingCircle·最小圆覆盖.hpp::smallest_enclosing_circle — 当前函数的已知消费者之一
 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> circum_center(Point<T> a, Point<T> b, Point<T> c) {
	Point<T> ab = b - a, ac = c - a;
	if (sgn(ab.cross(ac)) == 0) return Point<T>(NAN, NAN);
	Point<T> perpAB = Point<T>(-ab.y, ab.x);
	Point<T> perpAC = Point<T>(-ac.y, ac.x);

	Line<T> l1(a + ab * 0.5, a + ab * 0.5 + perpAB);
	Line<T> l2(a + ac * 0.5, a + ac * 0.5 + perpAC);

	return intersection(l1, l2);
}

/**
 * @brief 三角形内心 — 角平分线交点,到三边距离相等。
 * @complexity O(1)
 * @note 退化:任两边退化 (三角形面积为 0) 时公式仍返回某点,但无几何意义;勿用于非三角形输入。
 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> in_center(Point<T> a, Point<T> b, Point<T> c) {
	T ab = len(b - a);
	T bc = len(c - b);
	T ca = len(a - c);

	return (a * bc + b * ca + c * ab) / (ab + bc + ca);
}

} // namespace Geo2D
