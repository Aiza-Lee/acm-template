#pragma once
#include "1-1-Geo3D·三维浮点基础.hpp"

/*
 * 三维空间球与线/面/球求交
 *
 * Overview:
 *     提供三维空间中的球 (Sphere) 与空间交圆 (Circle3) 结构，以及空间直线与平面、
 *     两平面交线、三平面交点、球与直线交点、球与平面交圆、两球相交圆的完整求交函数。
 *
 * API:
 *     Sphere<T>(c, r)                — 三维球（球心 c，半径 r）。
 *     Sphere::contains(p)            — 点 p 是否在球内（含球面），O(1)。
 *     Sphere::on_sphere(p)           — 点 p 是否在球面上（容差 EPS），O(1)。
 *     Circle3<T>(c, r)               — 三维空间交圆（圆心 c，半径 r；valid() 判 r
 *                                       > 0）。
 *     closest_on_sphere(s, p)        — 球外点 p 到球面 s 的最近点；若 p
 *                                       在球内或球面上返回 std::nullopt，O(1)。
 *     line_plane_inter(l, pl)        — 空间直线 l 与平面 pl 的唯一交点；
 *                                       平行或直线在平面内返回 std::nullopt，O(1)。
 *     two_planes_inter(p1, p2)       — 两平面交线（方向向量已单位化）；
 *                                       平行或重合返回 std::nullopt，O(1)。
 *     three_planes_inter(p1, p2, p3) — 三平面唯一交点（Cramer 法则）；
 *                                       法向量共面退化返回 std::nullopt，O(1)。
 *     sphere_line_inter(s, l)        — 球 s 与直线 l 的交点列表（0/1/2 个，
 *                                       按直线参数 t 从小到大排列），O(1)。
 *     sphere_plane_inter(s, pl)      — 球 s 与平面 pl 的截面圆；相离时返回 r
 *                                       = -1（valid() == false），相切时 r = 0，
 *                                       O(1)。
 *     sphere_sphere_inter(s1, s2)    — 两球相交圆；相离、内含或同心退化时返回 r
 *                                       = -1，内外切时 r = 0，O(1)。
 *
 * Notes:
 *     1. Circle3::valid() 定义为 r > 0（严格交圆）；若需包含相切退化为单点的情形，
 *        请判断 sgn(cir.r) >= 0，此时 cir.c 即为切点。
 *
 * Related:
 *     1-1-Geo3D·三维浮点基础.hpp — 提供 Point、Line、Plane 及距离/投影运算。
 */

namespace Geo3D {

template<typename T = ld>
struct Circle3 {
    Point<T> c;
    T r;
    Circle3() : c(0, 0, 0), r(0) {}
    Circle3(Point<T> c, T r) : c(c), r(r) {}
    bool valid() const { return r > 0; }
};

template<typename T = ld>
struct Sphere {
    Point<T> c;
    T r;
    Sphere() : c(0, 0, 0), r(0) {}
    Sphere(Point<T> c, T r) : c(c), r(r) {}

    bool contains(Point<T> p) const { return cmp(c.dist2(p), r * r) <= 0; }
    bool on_sphere(Point<T> p) const { return is_zero(c.dist2(p) - r * r); }
};

// 球外点 p 到球面 s 的最近点；点在球内或球面上时返回 std::nullopt
template<typename T = ld>
std::optional<Point<T>> closest_on_sphere(const Sphere<T>& s, Point<T> p) {
    Point<T> v = s.c - p;
    T d2 = v.len2();
    if (cmp(d2, s.r * s.r) <= 0) return std::nullopt;
    return s.c - v * (s.r / std::sqrt(d2));
}

// 空间直线 l 与平面 pl 求唯一交点；平行或共面返回 std::nullopt
template<typename T = ld>
std::optional<Point<T>> line_plane_inter(const Line<T>& l, const Plane<T>& pl) {
    T den = l.dir.dot(pl.normal());
    if (is_zero(den)) return std::nullopt;
    return l.p0 - l.dir * (pl.eval(l.p0) / den);
}

// 两平面 p1, p2 求交线（选取最大非零 2x2 子行列式固定一维坐标为 0）；平行返回
// nullopt
template<typename T = ld>
std::optional<Line<T>> two_planes_inter(const Plane<T>& p1, const Plane<T>& p2) {
    Point<T> n1 = p1.normal(), n2 = p2.normal();
    Point<T> dir = n1.cross(n2);
    if (is_zero(dir.len2())) return std::nullopt;
    T dxy = n1.x * n2.y - n1.y * n2.x;
    T dyz = n1.y * n2.z - n1.z * n2.y;
    T dzx = n1.z * n2.x - n1.x * n2.z;
    Point<T> base;
    if (!is_zero(dxy)) {
        base = {(-p1.d * n2.y + p2.d * n1.y) / dxy,
                (-n1.x * p2.d + n2.x * p1.d) / dxy, 0};
    } else if (!is_zero(dyz)) {
        base = {0, (-p1.d * n2.z + p2.d * n1.z) / dyz,
                (-n1.y * p2.d + n2.y * p1.d) / dyz};
    } else if (!is_zero(dzx)) {
        base = {(-n1.z * p2.d + n2.z * p1.d) / dzx, 0,
                (-p1.d * n2.x + p2.d * n1.x) / dzx};
    } else {
        return std::nullopt;
    }
    return Line<T>(base, dir.unit());
}

// 三平面求唯一交点：利用标量三重积 D = (n1 x n2) . n3 求解；退化返回 nullopt
template<typename T = ld>
std::optional<Point<T>> three_planes_inter(
    const Plane<T>& p1, const Plane<T>& p2, const Plane<T>& p3
) {
    Point<T> n1 = p1.normal(), n2 = p2.normal(), n3 = p3.normal();
    T D = n1.triple(n2, n3);
    if (is_zero(D)) return std::nullopt;
    Point<T> c23 = n2.cross(n3), c31 = n3.cross(n1), c12 = n1.cross(n2);
    return (c23 * (-p1.d) + c31 * (-p2.d) + c12 * (-p3.d)) / D;
}

// 球 s 与空间直线 l 求交点（0/1/2 个，按直线参数 t 升序排列）
template<typename T = ld>
std::vector<Point<T>> sphere_line_inter(const Sphere<T>& s, const Line<T>& l) {
    Point<T> oc = l.p0 - s.c;
    T a = l.dir.len2(), b = 2 * oc.dot(l.dir), c = oc.len2() - s.r * s.r;
    T disc = b * b - 4 * a * c;
    if (sgn(disc) < 0) return {};
    if (is_zero(disc)) return {l.p0 + l.dir * (-b / (2 * a))};
    T sq = safe_sqrt(disc);
    T t1 = (-b - sq) / (2 * a), t2 = (-b + sq) / (2 * a);
    if (t1 > t2) std::swap(t1, t2);
    return {l.p0 + l.dir * t1, l.p0 + l.dir * t2};
}

// 球 s 与平面 pl 求交圆；无交点时返回 r = -1
template<typename T = ld>
Circle3<T> sphere_plane_inter(const Sphere<T>& s, const Plane<T>& pl) {
    T d = dist_plane(pl, s.c);
    T r2 = s.r * s.r - d * d;
    if (sgn(r2) < 0) return {s.c, -1};
    return {proj(pl, s.c), safe_sqrt(r2)};
}

// 两球 s1, s2 求交圆；同心、相离或内含时返回 r = -1
template<typename T = ld>
Circle3<T> sphere_sphere_inter(const Sphere<T>& s1, const Sphere<T>& s2) {
    Point<T> dd = s2.c - s1.c;
    T d2 = dd.len2();
    if (is_zero(d2)) return {s1.c, -1};
    T d = std::sqrt(d2);
    T a = (s1.r * s1.r - s2.r * s2.r + d2) / (2 * d);
    T h2 = s1.r * s1.r - a * a;
    if (sgn(h2) < 0) return {s1.c, -1};
    return {s1.c + dd * (a / d), safe_sqrt(h2)};
}

using SphereFP = Sphere<ld>;

} // namespace Geo3D
