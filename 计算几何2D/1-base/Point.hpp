#pragma once
#include "../0-utils/geo-utils.hpp"
namespace Geo2D {

template<typename T>
struct Point {
	T x, y;
	Point() : x(0), y(0) {}
	Point(T x, T y) : x(x), y(y) {}
		
	Point operator+(const Point& p) const { return Point(x + p.x, y + p.y); }
	Point operator-(const Point& p) const { return Point(x - p.x, y - p.y); }
	Point operator*(T k) const { return Point(x * k, y * k); }
	Point operator/(T k) const { return Point(x / k, y / k); }
	
	bool operator<(const Point& p) const { 
		int cx = cmp(x, p.x);
		if (cx != 0) return cx < 0;
		return cmp(y, p.y) < 0; 
	}
	bool operator==(const Point& p) const { 
		return cmp(x, p.x) == 0 && cmp(y, p.y) == 0; 
	}

	T len2() const { return x * x + y * y; }
	T dist2(const Point& p) const { return (*this - p).len2(); }
		
	T dot(const Point& p) const { return x * p.x + y * p.y; }
	T cross(const Point& p) const { return x * p.y - y * p.x; }

	// 极角所在半平面: 0 for [0, pi), 1 for [pi, 2pi)
	int half() const {
		int sy = sgn(y), sx = sgn(x);
		return (sy > 0 || (sy == 0 && sx > 0)) ? 0 : 1;
	}
};

// 极角排序
template<typename T>
bool polar_cmp(const Point<T>& a, const Point<T>& b) {
	int ha = a.half(), hb = b.half();
	if (ha != hb) return ha < hb;
	T c = a.cross(b);
	if (sgn(c) != 0) return sgn(c) > 0;
	return cmp(a.len2(), b.len2()) < 0;
}

// 旋转方向判定: 1(逆时针, 左转), -1(顺时针, 右转), 0(共线)
template<typename T>
int ccw(const Point<T>& p0, const Point<T>& p1, const Point<T>& p2) {
	return sgn((p1 - p0).cross(p2 - p0));
}

} // namespace Geo2D
