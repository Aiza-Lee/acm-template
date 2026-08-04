#pragma once
#include "Point·三维点.hpp"
#include "Plane·平面.hpp"

/*
 * 空间直线与线段
 *
 * Overview:
 * 	参数式空间直线 p0 + t·d,以及三维线段类型;含点 / 线段 / 异面直线相关运算。
 *
 * API:
 * 	Line<T>(p0, dir): 由一点和方向向量构造（dir 非零）。
 * 	Line<T>::from_points(p1, p2): 由两点构造,方向 = p2 - p1。
 * 	Segment<T>(a, b): 三维线段类型。
 * 	closest_point_on_line(line, p) -> Point<T>: p 到直线的最近点(无 t 范围约束)。浮点。O(1)。
 * 	distance_point_to_line(line, p) -> T: 点到直线的距离。浮点。O(1)。
 * 	distance_line_to_line(l1, l2) -> T: 两线最短距离;平行时降级为点到线距离。浮点。O(1)。
 * 	on_segment(p, seg) -> bool: 点是否在线段上（含端点）。O(1)。
 * 	closest_point_on_segment(p, seg) -> Point<T>: 点到线段的最近点（钳至端点）。浮点。O(1)。
 * 	distance_point_to_segment(p, seg) -> T: 点到线段的最短距离。浮点。O(1)。
 * 	line_plane_intersection(line, plane) -> optional<Point<T>>: 直线与平面交点;平行时 nullopt。浮点。O(1)。
 * 	LineFP / SegmentFP: 浮点常用别名。
 *
 * Notes:
 * 	closest_point_on_line 返回的最近点对应 t 可能为负（落在直线反向延长线）。
 * 	distance_line_to_line:两线平行 (cross = 0) 降级为 distance_point_to_line(l1, l2.p0)。
 * 	line_plane_intersection:|dir · normal| < EPS 时直线与平面平行,返回 nullopt。
 *
 * Related:
 * 	Plane·平面.hpp::eval / side: line_plane_intersection 的辅助判定。
 * 	PointFP·浮点三维点.hpp: 浮点长度 / 归一化的复用入口。
 */
namespace Geo3D {

template<typename T>
struct Line {
	Point<T> p0;   // 直线经过的一点
	Point<T> dir;  // 方向向量 (非零)

	Line() : p0(Point<T>()), dir(Point<T>(1, 0, 0)) {}
	Line(Point<T> p0, Point<T> dir) : p0(p0), dir(dir) {}

	static Line from_points(Point<T> p1, Point<T> p2) { return Line(p1, p2 - p1); }
};

template<typename T>
struct Segment {
	Point<T> a, b;
	Segment() : a(Point<T>()), b(Point<T>()) {}
	Segment(Point<T> a, Point<T> b) : a(a), b(b) {}
};

// p 在直线上的最近点(无最近 t 范围约束): t = (p - p0) · dir / |dir|^2
template<typename T>
requires std::is_floating_point_v<T>
Point<T> closest_point_on_line(const Line<T>& l, Point<T> p) {
	T t = (p - l.p0).dot(l.dir) / l.dir.len2();
	return l.p0 + l.dir * t;
}

// 点到直线距离:| (p - p0) × dir | / |dir|
template<typename T>
requires std::is_floating_point_v<T>
T distance_point_to_line(const Line<T>& l, Point<T> p) {
	Point<T> cross = (p - l.p0).cross(l.dir);
	return std::sqrt(cross.len2()) / std::sqrt(l.dir.len2());
}

// 两线最短距离:异面情形 = | (q - p) · (d1 × d2) | / |d1 × d2|;平行时降级为点到线。
template<typename T>
requires std::is_floating_point_v<T>
T distance_line_to_line(const Line<T>& l1, const Line<T>& l2) {
	Point<T> cross = l1.dir.cross(l2.dir);
	T den2 = cross.len2();
	if (is_zero(den2)) {                                 // 平行或重合
		return distance_point_to_line(l1, l2.p0);
	}
	Point<T> w = l1.p0 - l2.p0;
	return std::abs(w.dot(cross)) / std::sqrt(den2);
}

// 共线 + 投影落在 [a, b]:用 cross 与 dot 判定,浮点路径走 sgn / cmp。
template<typename T>
bool on_segment(Point<T> p, Segment<T> s) {
	Point<T> c = (s.a - p).cross(s.b - p);                      // 3D cross 是向量
	if constexpr (std::is_floating_point_v<T>) {
		if (sgn(c.x) != 0 || sgn(c.y) != 0 || sgn(c.z) != 0) return false;
	} else {
		if (c.x != 0 || c.y != 0 || c.z != 0) return false;
	}
	return sgn((s.a - p).dot(s.b - p)) <= 0;
}

// 点到线段最近点:先投影到参数 t ∈ [0,1]
template<typename T>
requires std::is_floating_point_v<T>
Point<T> closest_point_on_segment(Point<T> p, Segment<T> s) {
	Point<T> ab = s.b - s.a;
	T t = (p - s.a).dot(ab) / ab.len2();
	t = std::clamp(t, (T)0, (T)1);
	return s.a + ab * t;
}

template<typename T>
requires std::is_floating_point_v<T>
T distance_point_to_segment(Point<T> p, Segment<T> s) {
	Point<T> cp = closest_point_on_segment(p, s);
	return std::hypot(p.x - cp.x, p.y - cp.y, p.z - cp.z);
}

// 直线与平面交点:dir 与 normal 不平行时 t = -eval(p0) / (dir · normal)
template<typename T>
requires std::is_floating_point_v<T>
std::optional<Point<T>> line_plane_intersection(const Line<T>& l, const Plane<T>& pl) {
	T den = l.dir.x * pl.a + l.dir.y * pl.b + l.dir.z * pl.c;
	if (is_zero(den)) return std::nullopt;             // 平行
	T num = -(pl.a * l.p0.x + pl.b * l.p0.y + pl.c * l.p0.z + pl.d);
	T t = num / den;
	return l.p0 + l.dir * t;
}

using LineFP = Line<ld>;
using SegmentFP = Segment<ld>;

} // namespace Geo3D
