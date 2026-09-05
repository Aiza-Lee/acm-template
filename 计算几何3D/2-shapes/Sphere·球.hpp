#pragma once
#include "../1-base/Point·三维点.hpp"
#include "../1-base/Plane·平面.hpp"
#include "../1-base/PointFP·浮点三维点.hpp"

/*
 * 球
 *
 * Overview:
 *  提供球类型、点包含判定、与平面 / 球的求交,以及点到球外切点。
 *
 * API:
 *  Sphere<T>(c, r)                                      — 由圆心 c 与半径 r 构造;要求浮点类型。
 *  contains(s, p) -> bool                               — 点 p 是否在球内(含边界)。O(1)。
 *  on_sphere(s, p) -> bool                              — 点 p 是否在球面上(浮点容差)。O(1)。
 *  Circle<T>(c, r)                                      — 球求交返回的圆类型;r < 0 表示无效/无交。
 *  sphere_plane_intersection(s, pl) -> Circle<T>        — 球与平面交圆;球心到平面 |d| > r 时 r < 0。O(1)。
 *  sphere_sphere_intersection(s1, s2) -> Circle<T>      — 球与球交圆;返回 Circle.r 的判定规则:
 *      r = 0 外切, r > 0 相交, r < 0 包含 / 无意义(同心不等径 r=-1;同心等径 r=-1)。
 *  tangent_points_to_sphere(p, s) -> optional<Point<T>> — 从点 p 到球的外切点（两点中较近的一个）;点位于球内或球面上返回 nullopt。O(1)。
 *  SphereFP                                             — Sphere<ld> 的常用别名。
 *
 * Notes:
 *  sphere_sphere_intersection 返回的圆心位于两球心连线上,半径 = sqrt(r1² - a²),a = (r1² - r2² + d²) / (2d)。
 *  tangent_points_to_sphere 仅返回一个切点(较近 P 的那条);另一条需绕 PC 轴旋转 180°(当前接口不重复提供)。
 *
 * Related:
 *  Line·空间直线.hpp::Line: 由切点 T 与点 p 即可拼出切线 Line(p, T - p)。
 */
namespace Geo3D {

template<typename T>
requires std::is_floating_point_v<T>
struct Circle {
    Point<T> c;
    T r;
    Circle() : c(Point<T>()), r(0) {}
    Circle(Point<T> c, T r) : c(c), r(r) {}
    bool valid() const { return r > 0; }                         // r > 0 表有意义的相交圆;r == 0 表外切;r < 0 表无交
};

template<typename T>
requires std::is_floating_point_v<T>
struct Sphere {
    Point<T> c;
    T r;
    Sphere() : c(Point<T>()), r(0) {}
    Sphere(Point<T> c, T r) : c(c), r(r) {}
};

// 点 p 是否在球内(含边界):|p - c|² ≤ r²
template<typename T>
requires std::is_floating_point_v<T>
bool contains(const Sphere<T>& s, Point<T> p) {
    return sgn((p - s.c).len2() - s.r * s.r) <= 0;
}

// 点是否在球面上:容差比较
template<typename T>
requires std::is_floating_point_v<T>
bool on_sphere(const Sphere<T>& s, Point<T> p) {
    return is_zero((p - s.c).len2() - s.r * s.r);
}

// 球与平面交圆:半径 r = sqrt(r_sphere² - d²),d = 球心到平面有符号距离
template<typename T>
requires std::is_floating_point_v<T>
Circle<T> sphere_plane_intersection(const Sphere<T>& s, const Plane<T>& pl) {
    T d = distance_to_point(pl, s.c);
    T rd = s.r * s.r - d * d;
    if (rd < 0) return Circle<T>(s.c, -1);                       // 无交(r < 0)
    Point<T> n(pl.a, pl.b, pl.c);
    T nlen = std::hypot(n.x, n.y, n.z);
    Point<T> c = s.c - n * (d / nlen);                            // 球心在平面的投影
    T r = safe_sqrt(rd) / nlen;
    return Circle<T>(c, r);
}

// 球与球交圆:r = sqrt(r1² - a²);外切时 r = 0;同心 / 内含时 r < 0
template<typename T>
requires std::is_floating_point_v<T>
Circle<T> sphere_sphere_intersection(const Sphere<T>& s1, const Sphere<T>& s2) {
    Point<T> dd = s2.c - s1.c;
    T d2 = dd.len2();
    if (is_zero(d2)) {                                            // 同心
        return Circle<T>(s1.c, -1);                               // 同心无交 / 等径重合,统一 r = -1
    }
    T d = std::sqrt(d2);
    T a = (s1.r * s1.r - s2.r * s2.r + d2) / (2 * d);            // 交圆心在 d 方向距 s1.c 的距离
    T h2 = s1.r * s1.r - a * a;
    if (h2 < 0) return Circle<T>(s1.c, -1);                       // 内含无圆交
    T h = safe_sqrt(h2);
    Point<T> mid = s1.c + dd * (a / d);
    return Circle<T>(mid, h);
}

// 从外部点 p 到球的外切点（较近 p 的那个）;点位于球内或球面上时返回 nullopt
template<typename T>
requires std::is_floating_point_v<T>
std::optional<Point<T>> tangent_points_to_sphere(Point<T> p, const Sphere<T>& s) {
    Point<T> v = s.c - p;
    T v2 = v.len2();
    if (cmp(v2, s.r * s.r) <= 0) return std::nullopt;            // 球内或球上
    T len_v = std::sqrt(v2);
    return s.c - v * (s.r / len_v);                              // 沿 PC 方向投影到球面
}

using SphereFP = Sphere<ld>;

} // namespace Geo3D
