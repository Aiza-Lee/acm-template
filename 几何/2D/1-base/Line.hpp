#pragma once
#include "Point.hpp"

namespace Geo2D {

// 直线（一般式 ax + by + c = 0）
struct Line {
	ld a, b, c;
	
	Line() : a(0), b(0), c(0) {}
	Line(ld a, ld b, ld c) : a(a), b(b), c(c) {}

	// 两点构造直线 (p1 -> p2)
	// 如果需要用点和法向量构造，请使用 static 方法 from_point_and_normal
	Line(Point p1, Point p2) {
		a = p2.y - p1.y;
		b = p1.x - p2.x;
		c = -a * p1.x - b * p1.y;
	}
		
	// 点和方向向量
	static Line from_point_and_dir(Point p, Point dir) {
		return Line(p, p + dir);
	}
	
	// 点和法向量
	static Line from_point_and_normal(Point p, Point normal) {
		return Line(normal.x, normal.y, -(normal.x * p.x + normal.y * p.y));
	}
		
	// 获取直线的方向向量（单位化）
	Point direction() const {  return Point(-b, a).normalize(); }
	
	// 获取直线的法向量（单位化）
	Point normal() const { return Point(a, b).normalize(); }
		
	Point intersection(const Line& l) const {
		ld det = a * l.b - b * l.a;
		// 需注意 det 接近 0 的情况（平行）
		return Point((b * l.c - c * l.b) / det, (c * l.a - a * l.c) / det);
	}

	ld distance_to_point(Point p) const {
		return std::abs(a * p.x + b * p.y + c) / std::sqrt(a * a + b * b);
	}
	
	// 代入点坐标计算值: ax + by + c
	// 几何意义:
	// > 0: 点在直线法向量指向的一侧 (通常是 "右侧", 因为 (dy, -dx) 是方向向量顺时针旋转 90 度)
	// < 0: 点在直线法向量反向的一侧 (通常是 "左侧")
	// = 0: 点在直线上
	ld eval(Point p) const {
		return a * p.x + b * p.y + c;
	}

	// 判断点在直线的哪一侧
	// 1: "右侧" / 法向量侧
	// -1: "左侧" / 法向量反侧
	// 0: 直线上
	int side(Point p) const {
		ld val = eval(p);
		if (std::abs(val) < EPS) return 0;
		return val > 0 ? 1 : -1;
	}

	// 点在直线上的投影
	Point projection(Point p) const {
		Point n(a, b);
		return p - n * (eval(p) / n.len2());
	}

	// 点关于直线的对称点
	Point reflection(Point p) const {
		return p + (projection(p) - p) * 2;
	}
};

} // namespace Geo2D