#pragma once
#include "Point·三维点.hpp"

/*
 * 平面
 *
 * Overview:
 * 	提供一般式有向平面 ax + by + cz + d = 0 及位置 / 距离 / 投影运算。
 * 	半空间 {p | ax + by + cz + d ≥ 0} 与法向量 (a, b, c) 指向一致。
 *
 * API:
 * 	Plane<T>(a, b, c, d): 一般式构造。
 * 	Plane<T>::from_point_and_normal(p, n): 由一点和法向量构造。
 * 	Plane<T>::from_points(p1, p2, p3): 三点构造;(p2-p1)×(p3-p1) 作为法向量。
 * 	eval(p) -> T: 代入 ax + by + cz + d。O(1)。
 * 	side(p) -> int: 返回 +1 / 0 / -1（正 / 平面 / 负），使用 sgn。O(1)。
 * 	normal() -> Point<T>: 法向量 (a, b, c)；浮点要求归一化请在外部用 normalize。
 * 	distance_to_point(plane, p) -> T: 有符号距离（正值在法向量一侧）。浮点。O(1)。
 * 	projection(plane, p) -> Point<T>: 点在平面上的正交投影。浮点。O(1)。
 * 	reflection(plane, p) -> Point<T>: 点关于平面的对称点。浮点。O(1)。
 * 	PlaneFP: Plane<ld> 的常用别名。
 *
 * Notes:
 * 	模板参数 T: eval / side 对整数保持严格;浮点路径走 sgn / EPS。
 * 	from_points 要求三点不共线,否则法向量零向量,投影 / 距离公式无意义,勿用。
 * 	distance_to_point 为有符号距离:正值在法向量 (a, b, c) 一侧。
 *
 * Related:
 * 	PointFP·浮点三维点.hpp::normalize: 法向量归一化时使用。
 * 	Tetrahedron·四面体.hpp::point_in_tet: 以四个平面 face 符号判定点是否在四面体内。
 */
namespace Geo3D {

template<typename T>
struct Plane {
	T a, b, c, d;

	Plane() : a(0), b(0), c(0), d(0) {}
	Plane(T a, T b, T c, T d) : a(a), b(b), c(c), d(d) {}

	static Plane from_point_and_normal(Point<T> p, Point<T> n) {
		return Plane(n.x, n.y, n.z, -(n.x * p.x + n.y * p.y + n.z * p.z));
	}

	// 三点构造;法向量方向 = (p2 - p1) × (p3 - p1)
	Plane(Point<T> p1, Point<T> p2, Point<T> p3) {
		Point<T> n = (p2 - p1).cross(p3 - p1);
		a = n.x; b = n.y; c = n.z; d = -(a * p1.x + b * p1.y + c * p1.z);
	}

	T eval(Point<T> p) const {                    // ax + by + cz + d
		return a * p.x + b * p.y + c * p.z + d;
	}

	int side(Point<T> p) const {                  // +1 / 0 / -1
		return sgn(eval(p));
	}

	Point<T> normal() const {                     // 法向量 (a, b, c),未归一
		return Point<T>(a, b, c);
	}
};

template<typename T>
requires std::is_floating_point_v<T>
T distance_to_point(const Plane<T>& pl, Point<T> p) {
	T nl = std::hypot(pl.a, pl.b, pl.c);
	return pl.eval(p) / nl;
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> projection(const Plane<T>& pl, Point<T> p) {
	Point<T> n(pl.a, pl.b, pl.c);
	return p - n * (pl.eval(p) / n.len2());
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> reflection(const Plane<T>& pl, Point<T> p) {
	return p + (projection(pl, p) - p) * 2;
}

using PlaneFP = Plane<ld>;

} // namespace Geo3D
