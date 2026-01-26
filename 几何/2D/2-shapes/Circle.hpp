#pragma once
#include "../1-base/Point.hpp"
#include "../1-base/Line.hpp"
namespace Geo2D {

struct Circle {
	Point c; ld r;
	Circle() : c(Point()), r(0) {}
	Circle(Point c, ld r) : c(c), r(r) {}
	bool contains(Point p) const { return dcmp(c.dist_to(p), r) <= 0; }
	// 圆周上的点，angle为弧度
	Point point(ld angle) const {
		return Point(c.x + r * cos(angle), c.y + r * sin(angle));
	}
};

/**
 * @brief 圆与直线交点
 */
std::vector<Point> circle_line_intersection(Circle c, Line l) {
	std::vector<Point> res;
	ld d = l.distance_to_point(c.c);
		
	if (sign(d - c.r) > 0) return res;  // 无交点
		
	// 直线垂足
	Point normal(l.a, l.b);
	Point p0 = c.c - normal * (normal.dot(c.c) + l.c) / normal.len2();
		
	if (sign(d - c.r) == 0) {  // 一个交点（切线）
		res.push_back(p0);
		return res;
	}
		
	// 两个交点
	ld len = sqrt(c.r * c.r - d * d);
	Point dir = l.direction() * len;
	res.push_back(p0 + dir);
	res.push_back(p0 - dir);
	return res;
}

/**
 * @brief 圆与圆交点
 * note:
 * 		1.不处理两圆完全重合的情况（返回空）。
 */
std::vector<Point> circle_circle_intersection(Circle c1, Circle c2) {
	std::vector<Point> res;
	ld d = c1.c.dist_to(c2.c);
		
	if (sign(d - c1.r - c2.r) > 0 || sign(d - abs(c1.r - c2.r)) < 0) 
		return res; // 相离或内含
		
	if (sign(d) == 0 && sign(c1.r - c2.r) == 0) 
		return res; // 完全重合
		
	ld a = (c1.r * c1.r - c2.r * c2.r + d * d) / (2 * d);
	ld h = sqrt(std::max(0.0L, c1.r * c1.r - a * a));
		
	Point v = (c2.c - c1.c).normalize();
	Point mid = c1.c + v * a;
	Point normal = v.rotate(PI/2) * h;
		
	res.push_back(mid + normal);
	if (sign(h) != 0) res.push_back(mid - normal);
	return res;
}

/**
 * @brief 过点 p 做圆 c 的切线，返回切点。
 * 
 * note: 
 * 		1.点在圆内无切线；点在圆上有一个切点（即点本身）；点在圆外有两个切点。
 */
std::vector<Point> tangents_point_circle(Point p, Circle c) {
	std::vector<Point> res;
	ld d = c.c.dist_to(p);
	int s = sign(d - c.r);
	if (s < 0) return res; // 点在圆内
	if (s == 0) { // 点在圆上
		res.push_back(p);
		return res;
	}
	
	ld angle = asin(c.r / d);
	Point v = (p - c.c).normalize() * c.r;
	res.push_back(c.c + v.rotate(angle));
	res.push_back(c.c + v.rotate(-angle));
	return res;
}

/**
 * @brief 计算两个圆的公切线，返回切点对。
 * 
 * @return std::vector<std::pair<Point, Point>> 切点对列表
 *         每个 pair 包含两个点 {p1, p2}，分别位于 c1 和 c2 上。
 *         列表顺序通常为：外公切线1, 外公切线2, 内公切线1, 内公切线2。
 */
std::vector<std::pair<Point, Point>> tangents_circle_circle(Circle c1, Circle c2) {
	std::vector<std::pair<Point, Point>> res;
	if (c1.r < c2.r) {
		auto tmp = tangents_circle_circle(c2, c1);
		for (auto p : tmp) res.push_back({p.second, p.first});
		return res;
	}
	
	ld d2 = c1.c.dist_to(c2.c);
	ld rdiff = c1.r - c2.r;
	ld rsum = c1.r + c2.r;
	
	if (dcmp(d2, rdiff * rdiff) < 0) return res; // 内含
	
	Point v = (c2.c - c1.c).normalize();
	
	// 外公切线
	ld angle_out = acos(rdiff / sqrt(d2)); // 修正：d2是距离平方还是距离？dist_to返回距离
	// 修正：d2应该是距离
	ld d = sqrt(d2); // dist_to返回的是距离，不是平方
	// 重新写一下
	d = c1.c.dist_to(c2.c);
	if (sign(d) == 0 && sign(c1.r - c2.r) == 0) return res; // 重合，无数条切线

	// 外公切线
	if (sign(d - rdiff) >= 0) {
		ld angle = acos(rdiff / d);
		Point v1 = v.rotate(angle);
		Point v2 = v.rotate(-angle);
		res.push_back({c1.c + v1 * c1.r, c2.c + v1 * c2.r});
		if (sign(angle) != 0)
			res.push_back({c1.c + v2 * c1.r, c2.c + v2 * c2.r});
	}

	// 内公切线
	if (sign(d - rsum) >= 0) {
		ld angle = acos(rsum / d);
		Point v1 = v.rotate(angle);
		Point v2 = v.rotate(-angle);
		res.push_back({c1.c + v1 * c1.r, c2.c + v1 * (-c2.r)}); // 注意内切线方向相反
		if (sign(angle) != 0)
			res.push_back({c1.c + v2 * c1.r, c2.c + v2 * (-c2.r)});
	}
	
	return res;
}

} // namespace Geo2D