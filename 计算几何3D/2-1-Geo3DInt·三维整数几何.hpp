#pragma once
#include "aizalib.h"

/*
 * 三维整数几何
 *
 * Overview:
 *     提供纯整数精确三维计算几何：三维整数点/向量、点积/叉积/标量三重积、
 *     空间线段包含判定、三角形 4 倍面积平方、四面体 6 倍有向体积与严格体内判定。
 *
 * API:
 *     sgn(x)                   — 严格整数符号函数，返回 -1、0、1，O(1)。
 *     Point<T>(x, y, z)        — 三维整数点/向量，默认 T = i64。
 *     Point::dot(p) / cross(p) — 精确点积 / 叉积向量，O(1)。
 *     Point::triple(p, q)      — 标量三重积 (this x p) . q，O(1)。
 *     Point::len2() / dist2(p) — 模长平方 / 距离平方，O(1)。
 *     Segment<T>(a, b)         — 三维整数线段。
 *     on_segment(s, p)         — 点 p 是否在线段 s 上（含两端点），O(1)。
 *     triangle_area4(a, b, c)  — 空间三角形 4 倍面积平方 |(b-a)x(c-a)|^2，O(1)。
 *     tet_volume6(a, b, c, d)  — 四面体 6 倍有向体积 ((b-a)x(c-a)).(d-a)，O(1)。
 *     in_tet(a, b, c, d, p)    — 点在四面体内精确判定：严格内 1、面/退化 0、外 -1，
 *                                 O(1)。
 *     closest_pair_3d(pts)     — 求三维点集最近点对的距离平方与下标对，O(N^2)。
 *
 * Notes:
 *     1. 全程零浮点误差，默认 T = i64（|坐标| <= 1e6 时 triple 不溢出；更大坐标可传
 *        i128）。
 *     2. in_tet 通过 4 个子四面体的有向体积符号与总有向体积符号比对实现，
 *        与顶点朝向无关。
 *
 * Related:
 *     1-1-Geo3D·三维浮点基础.hpp — 三维浮点几何版本。
 */

namespace Geo3DInt {

template<typename T = i64>
inline int sgn(T x) {
    return (x > 0) - (x < 0);
}

template<typename T = i64>
struct Point {
    T x, y, z;
    Point() : x(0), y(0), z(0) {}
    Point(T x, T y, T z) : x(x), y(y), z(z) {}

    Point operator+(const Point& p) const { return {x + p.x, y + p.y, z + p.z}; }
    Point operator-(const Point& p) const { return {x - p.x, y - p.y, z - p.z}; }
    Point operator-() const { return {-x, -y, -z}; }
    Point operator*(T k) const { return {x * k, y * k, z * k}; }
    Point operator/(T k) const { return {x / k, y / k, z / k}; }

    bool operator<(const Point& p) const {
        if (x != p.x) return x < p.x;
        return y != p.y ? y < p.y : z < p.z;
    }
    bool operator==(const Point& p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!=(const Point& p) const { return !(*this == p); }

    T dot(const Point& p) const { return x * p.x + y * p.y + z * p.z; }
    Point cross(const Point& p) const {
        return {y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x};
    }
    T triple(const Point& p, const Point& q) const { return cross(p).dot(q); }
    T len2() const { return x * x + y * y + z * z; }
    T dist2(const Point& p) const { return (*this - p).len2(); }
};

template<typename T = i64>
struct Segment {
    Point<T> a, b;
    Segment() = default;
    Segment(Point<T> a, Point<T> b) : a(a), b(b) {}
};

// 纯整数判定空间点 p 是否在线段 s 上（叉积为零向量且点积 <= 0）
template<typename T = i64>
bool on_segment(const Segment<T>& s, Point<T> p) {
    Point<T> c = (s.a - p).cross(s.b - p);
    if (c.x != 0 || c.y != 0 || c.z != 0) return false;
    return (s.a - p).dot(s.b - p) <= 0;
}

// 空间三角形面积平方的 4 倍 |(b - a) x (c - a)|^2
template<typename T = i64>
T triangle_area4(Point<T> a, Point<T> b, Point<T> c) {
    return (b - a).cross(c - a).len2();
}

// 四面体有向体积的 6 倍 ((b - a) x (c - a)) . (d - a)
template<typename T = i64>
T tet_volume6(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    return (b - a).cross(c - a).dot(d - a);
}

// 纯整数判定点 p 是否在四面体 abcd 内：1 严格内部，0 边界面或退化，-1 严格外部
template<typename T = i64>
int in_tet(Point<T> a, Point<T> b, Point<T> c, Point<T> d, Point<T> p) {
    T v = tet_volume6(a, b, c, d);
    if (v == 0) return 0;
    int sv = sgn(v);
    int sa = sgn(tet_volume6(p, b, c, d)) * sv;
    int sb = sgn(tet_volume6(a, p, c, d)) * sv;
    int sc = sgn(tet_volume6(a, b, p, d)) * sv;
    int sd = sgn(tet_volume6(a, b, c, p)) * sv;
    if (sa < 0 || sb < 0 || sc < 0 || sd < 0) return -1;
    if (sa == 0 || sb == 0 || sc == 0 || sd == 0) return 0;
    return 1;
}

// 求三维点集最近点对，返回 (最小距离平方, 下标 i, 下标 j)，O(N^2)
template<typename T = i64>
std::tuple<T, int, int> closest_pair_3d(const std::vector<Point<T>>& pts) {
    int n = pts.size();
    if (n <= 1) return {-1, -1, -1};
    T best = -1;
    int ia = -1, ib = -1;
    rep(i, 0, n - 1) rep(j, i + 1, n - 1) {
        T d2 = pts[i].dist2(pts[j]);
        if (best < 0 || d2 < best) {
            best = d2;
            ia = i;
            ib = j;
        }
    }
    return {best, ia, ib};
}

} // namespace Geo3DInt
