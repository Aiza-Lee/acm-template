#pragma once
#include "../1-base/Line.hpp"

namespace Geo2D {

/**
 * @brief 三角形外心
 * 
 * 外接圆圆心：到三角形三个顶点距离相等。
 * 
 * @param a 顶点 A
 * @param b 顶点 B
 * @param c 顶点 C
 * @return Point<T> 外心坐标
 * @complexity O(1)
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
 * @brief 三角形内心
 * 
 * 内切圆圆心：角平分线交点，到三边距离相等。
 * 
 * @param a 顶点 A
 * @param b 顶点 B
 * @param c 顶点 C
 * @return Point<T> 内心坐标
 * @complexity O(1)
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
