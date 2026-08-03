#pragma once
#include "../1-base/Line·直线.hpp"

/*
 * 三角形特征点
 *
 * Overview:
 * 	提供三角形外心与内心的直接计算。
 *
 * API:
 * 	circum_center(a, b, c) -> Point<T>: 外接圆圆心。O(1)。
 * 	in_center(a, b, c) -> Point<T>: 角平分线交点。O(1)。
 *
 * Notes:
 * 	仅支持浮点坐标。
 * 	circum_center: 三顶点共线或三角形退化为线段时返回 (NaN, NaN)。
 * 	in_center: 任两边退化（面积为 0）时公式仍返回某点，但无几何意义，勿用于非三角形输入。
 *
 * Related:
 * 	MinEnclosingCircle·最小圆覆盖.hpp::smallest_enclosing_circle: 当前 circum_center 的已知消费者之一。
 */
namespace Geo2D {

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

template<typename T>
requires std::is_floating_point_v<T>
Point<T> in_center(Point<T> a, Point<T> b, Point<T> c) {
	T ab = len(b - a);
	T bc = len(c - b);
	T ca = len(a - c);

	return (a * bc + b * ca + c * ab) / (ab + bc + ca);
}

} // namespace Geo2D
