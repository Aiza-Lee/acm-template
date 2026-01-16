#pragma once
#include "../utils/geo-utils.hpp"
namespace Geo2D {

struct Point {
	ld x, y;
	Point(ld x = 0, ld y = 0) : x(x), y(y) {}
		
	Point operator+(const Point& p) const { return Point(x + p.x, y + p.y); }
	Point operator-(const Point& p) const { return Point(x - p.x, y - p.y); }
	Point operator*(ld k) const { return Point(x * k, y * k); }
	Point operator/(ld k) const { return Point(x / k, y / k); }
	bool operator<(const Point& p) const { return dcmp(x, p.x) == 0 ? dcmp(y, p.y) < 0 : dcmp(x, p.x) < 0; }
	bool operator==(const Point& p) const { return dcmp(x, p.x) == 0 && dcmp(y, p.y) == 0; }

	ld len() const { return hypot(x, y); }
	ld len2() const { return x * x + y * y; }
	ld dist_to(const Point& p) const { return (*this - p).len(); }
	Point normalize() const { ld l = len(); return Point(x / l, y / l); }

	// 向量旋转，rad为弧度，逆时针为正
	Point rotate(ld rad) const {
		return Point(x * cos(rad) - y * sin(rad), x * sin(rad) + y * cos(rad));
	}
		
	// 向量的极角 (-PI, PI]
	ld angle() const { return atan2(y, x); }

	// 计算从this向量到p向量的有向夹角，范围为(-PI, PI]
	ld directed_angle(const Point& p) const {
		ld angle = atan2(cross(p), dot(p));
		return angle;
	}
		
	ld dot(const Point& p) const { return x * p.x + y * p.y; }
	ld cross(const Point& p) const { return x * p.y - y * p.x; }
};

} // namespace Geo2D