#pragma once
#include "Point.hpp"
namespace Geo2D {

// 直线（一般式 ax + by + c = 0）
struct Line {
	ld a, b, c;
	Line(Point p1, Point p2) {
		a = p2.y - p1.y;
		b = p1.x - p2.x;
		c = p2.x * p1.y - p1.x * p2.y;
	}
	Line(Point p, Point normal) {
		a = normal.x;
		b = normal.y;
		c = -a * p.x - b * p.y;
	}
		
	// 点和方向向量
	static Line from_point_and_dir(Point p, Point dir) {
		return Line(p, p + dir);
	}
		
	// 获取直线的方向向量（单位化）
	Point direction() const {  return Point(-b, a).normalize(); }
		
	Point intersection(const Line& l) const {
		ld det = a * l.b - b * l.a;
		return Point((b * l.c - c * l.b) / det, (c * l.a - a * l.c) / det);
	}
	ld distance_to_point(Point p) const {
		return abs(a * p.x + b * p.y + c) / sqrt(a * a + b * b);
	}
	
	// 点在直线上的投影
	Point projection(Point p) const {
		Point normal(a, b);
		return p - normal * (a * p.x + b * p.y + c) / normal.len2();
	}
};

/**
 * @brief 计算点 p 在直线 l 上的正交投影点（垂足）。
 */
Point projection(Point p, Line l) {
	Point normal(l.a, l.b);
	return p - normal * (normal.dot(p) + l.c) / normal.len2();
}

/**
 * @brief 计算点 p 关于直线 l 的轴对称点。
 */
Point reflection(Point p, Line l) {
	return p + (projection(p, l) - p) * 2;
}

} // namespace Geo2D