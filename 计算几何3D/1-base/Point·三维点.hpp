#pragma once
#include "../0-utils/GeoUtils·几何工具.hpp"

/*
 * 三维点
 *
 * Overview:
 * 	提供整数 / 浮点通用的三维点（向量）类型及基础运算。
 *
 * API:
 * 	Point<T>(x, y, z): 构造三维点 / 向量。
 * 	operator + / -: 向量加减。O(1)。
 * 	operator * k / / k: 数乘和数除（数除仅在浮点下数值稳定）。O(1)。
 * 	operator < / ==: 按 x、y、z 字典序比较 / 判等。O(1)。
 * 	len2() / dist2(p): 模长平方 / 两点距离平方。O(1)。
 * 	dot(p): 点积,返回 T。O(1)。
 * 	cross(p): 叉积,返回 Point<T>。O(1)。
 * 	triple(p, q): 三重标量积 (this × p) · q; 用于共面判定 / 有向体积。O(1)。
 *
 * Notes:
 * 	浮点比较经 cmp / sgn 使用 EPS。
 * 	整数坐标下,叉积与标量三重积需调用方自行控制 T 溢出（|x|,|y|,|z| ≤ 1e6 时 i64 安全）。
 *
 * Related:
 * 	PointFP·浮点三维点.hpp::normalize: 浮点归一化前需 AST(!is_zero(p.len2()))。
 */
namespace Geo3D {

template<typename T>
struct Point {
	T x, y, z;

	Point() : x(0), y(0), z(0) {}
	Point(T x, T y, T z) : x(x), y(y), z(z) {}

	Point operator+(const Point& p) const { return Point(x + p.x, y + p.y, z + p.z); }  // 向量加
	Point operator-(const Point& p) const { return Point(x - p.x, y - p.y, z - p.z); }  // 向量减
	Point operator*(T k)        const { return Point(x * k, y * k, z * k); }           // 数乘
	Point operator/(T k)        const { return Point(x / k, y / k, z / k); }           // 数除

	bool operator<(const Point& p) const {  // 字典序 (先 x 后 y 后 z),经 cmp 兼容整数严格 / 浮点 EPS
		int c = cmp(x, p.x);
		if (c != 0) return c < 0;
		c = cmp(y, p.y);
		if (c != 0) return c < 0;
		return cmp(z, p.z) < 0;
	}
	bool operator==(const Point& p) const {  // 经 cmp 容差后判等
		return cmp(x, p.x) == 0 && cmp(y, p.y) == 0 && cmp(z, p.z) == 0;
	}

	T len2()                const { return x * x + y * y + z * z; }                    // |p|²
	T dist2(const Point& p) const { return (*this - p).len2(); }                       // |p - q|²

	T dot(const Point& p)   const { return x * p.x + y * p.y + z * p.z; }             // 点积 p·q
	Point cross(const Point& p) const {                                              // 叉积 p × q (向量)
		return Point(
			y * p.z - z * p.y,
			z * p.x - x * p.z,
			x * p.y - y * p.x
		);
	}
	T triple(const Point& p, const Point& q) const {                                   // (this × p) · q
		return cross(p).dot(q);
	}
};

} // namespace Geo3D
