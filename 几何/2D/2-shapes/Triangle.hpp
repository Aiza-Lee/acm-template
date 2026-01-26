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
 * @return Point 外心坐标
 * @complexity O(1)
 */
Point circum_center(Point a, Point b, Point c) {
	Point ab = b - a, ac = c - a;
	if (sign(ab.cross(ac)) == 0) return Point(NAN, NAN);
	Point perpAB = Point(-ab.y, ab.x);
	Point perpAC = Point(-ac.y, ac.x);
		
	Line l1(a + ab * 0.5, a + ab * 0.5 + perpAB);
	Line l2(a + ac * 0.5, a + ac * 0.5 + perpAC);
		
	return l1.intersection(l2);
}

/**
 * @brief 三角形内心
 * 
 * 内切圆圆心：角平分线交点，到三边距离相等。
 * 
 * @param a 顶点 A
 * @param b 顶点 B
 * @param c 顶点 C
 * @return Point 内心坐标
 * @complexity O(1)
 */
Point in_center(Point a, Point b, Point c) {
	ld ab = (b - a).len();
	ld bc = (c - b).len();
	ld ca = (a - c).len();
		
	return (a * bc + b * ca + c * ab) / (ab + bc + ca);
}

} // namespace Geo2D
