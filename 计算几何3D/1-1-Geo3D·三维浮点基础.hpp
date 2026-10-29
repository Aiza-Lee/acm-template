#pragma once
#include "aizalib.h"

/*
 * 三维浮点几何基础
 *
 * Overview:
 *     提供三维浮点计算几何的核心图元与度量底座：三维点/向量、点向式空间直线与线段、
 *     一般式有向平面（含点到平面距离、正交投影与镜像反射），
 *     以及空间三角形与四面体的面积、体积、重心坐标与包含判定。
 *
 * API:
 *     EPS / PI                             — 浮点容差 (1e-10) / 圆周率常数。
 *     sgn(x) / cmp(x, y) / is_zero(x)      — 容差符号 (-1/0/1)、三路比较与判零，
 *                                             O(1)。
 *     clamp_unit(x) / safe_sqrt(x)         — 截断到 [-1, 1] / 非负安全开方，O(1)。
 *     safe_acos(x) / safe_asin(x)          — 防越界反余弦 / 反正弦，O(1)。
 *     Point<T>(x, y, z)                    — 三维点/向量，默认 T = ld。
 *     Point::dot(p) / cross(p)             — 点积 / 叉积向量，O(1)。
 *     Point::triple(p, q)                  — 标量三重积 (this x p) . q，O(1)。
 *     Point::len2() / len() / unit()       — 模长平方 / 模长 (hypot) / 单位向量，
 *                                             O(1)。
 *     Point::dist2(p) / dist(p) / ang_xy() — 距离平方 / 距离 / XY 投影极角，O(1)。
 *     angle_between(a, b)                  — 两三维向量无符号夹角 [0, PI]，O(1)。
 *     Line<T>(p0, dir) / Segment<T>(a, b)  — 点向式空间直线 p0 + t*dir 与三维线段。
 *     Line::from_points(a, b)              — 由两点构造空间直线 Line(a, b - a)，
 *                                             O(1)。
 *     closest_on_line(l, p) / dist_line    — 点 p 在直线 l 上的最近点 / 最短距离，
 *                                             O(1)。
 *     dist_line_line(l1, l2)               — 两空间直线最短距离（含异面与平行），
 *                                             O(1)。
 *     on_segment(s, p)                     — 点 p 是否在线段 s 上（含两端点），
 *                                             O(1)。
 *     closest_on_seg(s, p) / dist_seg      — 点 p 到线段 s 的最近点 / 最短距离，
 *                                             O(1)。
 *     Plane<T>(a, b, c, d)                 — 一般式有向平面 ax + by + cz + d = 0。
 *     Plane::from_point_normal(p, n)       — 由点 p 与法向量 n 构造有向平面，O(1)。
 *     Plane::from_points(p1, p2, p3)       — 由三点构造平面，法向 (p2-p1)x(p3-p1)，
 *                                             O(1)。
 *     Plane::eval(p) / side(p) / normal()  — 代入求值 / 符号 (-1/0/1) / 法向量，
 *                                             O(1)。
 *     dist_plane(pl, p)                    — 点 p 到平面的有符号距离（法向侧为
 *                                             正），O(1)。
 *     proj(pl, p) / refl(pl, p)            — 点 p 在平面上的正交投影 / 镜像对称点，
 *                                             O(1)。
 *     triangle_area(a, b, c)               — 空间三角形面积，O(1)。
 *     triangle_normal(a, b, c)             — 空间三角形单位法向量（共线退化返回零向
 *                                             量），O(1)。
 *     barycentric(a, b, c, p, al, be, ga)  — 求点 p 关于三角形 abc 的重心坐标，
 *                                             O(1)。
 *     in_triangle(a, b, c, p)              — 点在空间三角形内：严格内 1、边/顶点
 *                                             0、外或不共面 -1，O(1)。
 *     tet_volume6(a, b, c, d)              — 四面体 6 倍有向体积 ((b-a)
 *                                             x(c-a)).(d-a)，O(1)。
 *     tet_volume(a, b, c, d)               — 四面体无符号体积 |tet_volume6| / 6，
 *                                             O(1)。
 *     tet_centroid(a, b, c, d)             — 四面体几何质心 (a+b+c+d)/4，O(1)。
 *     in_tet(a, b, c, d, p)                — 点在四面体内：严格内 1、边界面/退化
 *                                             0、外 -1，O(1)。
 *
 * Notes:
 *     1. 本文件提供三维浮点底座（默认 T = ld）；球与线/面/球求交请见
 *        1-2-Intersections3D·空间球与求交.hpp，纯整数三维几何请见 Geo3DInt。
 *     2. 参数顺序统一遵循 (形状/线/面, 查询点)，如 dist_line(l, p)、proj(pl, p)。
 *
 * Related:
 *     1-2-Intersections3D·空间球与求交.hpp — 球、空间圆及线/面/球求交。
 *     1-3-Rotation3D·三维旋转.hpp          — Rodrigues 绕轴旋转与欧拉角变换。
 *     1-4-ConvexHull3D·三维凸包.hpp        — O(N^2) 增量法三维凸包。
 */

namespace Geo3D {

inline constexpr ld EPS = 1e-10;
inline constexpr ld PI = std::numbers::pi_v<ld>;

template<typename T = ld>
inline int sgn(T x) {
    return x < -EPS ? -1 : (x > EPS ? 1 : 0);
}

template<typename T = ld>
inline int cmp(T x, T y) {
    return sgn(x - y);
}

template<typename T = ld>
inline bool is_zero(T x) {
    return sgn(x) == 0;
}

template<typename T = ld>
inline T clamp_unit(T x) {
    return std::clamp(x, (T)-1, (T)1);
}

template<typename T = ld>
inline T safe_sqrt(T x) {
    return std::sqrt(std::max((T)0, x));
}

template<typename T = ld>
inline T safe_acos(T x) {
    return std::acos(clamp_unit(x));
}

template<typename T = ld>
inline T safe_asin(T x) {
    return std::asin(clamp_unit(x));
}

template<typename T = ld>
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
        int c = cmp(x, p.x);
        if (c != 0) return c < 0;
        c = cmp(y, p.y);
        return c != 0 ? c < 0 : cmp(z, p.z) < 0;
    }
    bool operator==(const Point& p) const {
        return cmp(x, p.x) == 0 && cmp(y, p.y) == 0 && cmp(z, p.z) == 0;
    }
    bool operator!=(const Point& p) const { return !(*this == p); }

    T dot(const Point& p) const { return x * p.x + y * p.y + z * p.z; }
    Point cross(const Point& p) const {
        return {y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x};
    }
    T triple(const Point& p, const Point& q) const { return cross(p).dot(q); }

    T len2() const { return x * x + y * y + z * z; }
    T len() const { return std::hypot(x, y, z); }
    T dist2(const Point& p) const { return (*this - p).len2(); }
    T dist(const Point& p) const { return (*this - p).len(); }

    Point unit() const {
        T l = len();
        AST(!is_zero(l));
        return {x / l, y / l, z / l};
    }
    T ang_xy() const { return std::atan2(y, x); }
};

template<typename T = ld>
T angle_between(const Point<T>& a, const Point<T>& b) {
    return safe_acos(a.dot(b) / (a.len() * b.len()));
}

template<typename T = ld>
struct Line {
    Point<T> p0, dir;
    Line() : p0(0, 0, 0), dir(1, 0, 0) {}
    Line(Point<T> p0, Point<T> dir) : p0(p0), dir(dir) {}
    static Line from_points(Point<T> a, Point<T> b) { return {a, b - a}; }
};

template<typename T = ld>
struct Segment {
    Point<T> a, b;
    Segment() = default;
    Segment(Point<T> a, Point<T> b) : a(a), b(b) {}
};

// 点 p 在空间直线 l 上的正交投影点
template<typename T = ld>
Point<T> closest_on_line(const Line<T>& l, Point<T> p) {
    T t = (p - l.p0).dot(l.dir) / l.dir.len2();
    return l.p0 + l.dir * t;
}

template<typename T = ld>
T dist_line(const Line<T>& l, Point<T> p) {
    return (p - l.p0).cross(l.dir).len() / l.dir.len();
}

// 两空间直线最短距离（平行或重合时退化为点到直线距离）
template<typename T = ld>
T dist_line_line(const Line<T>& l1, const Line<T>& l2) {
    Point<T> cr = l1.dir.cross(l2.dir);
    T den2 = cr.len2();
    if (is_zero(den2)) return dist_line(l1, l2.p0);
    return std::abs((l1.p0 - l2.p0).dot(cr)) / std::sqrt(den2);
}

template<typename T = ld>
bool on_segment(const Segment<T>& s, Point<T> p) {
    Point<T> c = (s.a - p).cross(s.b - p);
    if (!is_zero(c.x) || !is_zero(c.y) || !is_zero(c.z)) return false;
    return sgn((s.a - p).dot(s.b - p)) <= 0;
}

template<typename T = ld>
Point<T> closest_on_seg(const Segment<T>& s, Point<T> p) {
    Point<T> ab = s.b - s.a;
    if (is_zero(ab.len2())) return s.a;
    T t = std::clamp((p - s.a).dot(ab) / ab.len2(), (T)0, (T)1);
    return s.a + ab * t;
}

template<typename T = ld>
T dist_seg(const Segment<T>& s, Point<T> p) {
    return p.dist(closest_on_seg(s, p));
}

template<typename T = ld>
struct Plane {
    T a, b, c, d;
    Plane() : a(0), b(0), c(0), d(0) {}
    Plane(T a, T b, T c, T d) : a(a), b(b), c(c), d(d) {}
    Plane(Point<T> p1, Point<T> p2, Point<T> p3) {
        Point<T> n = (p2 - p1).cross(p3 - p1);
        a = n.x; b = n.y; c = n.z; d = -n.dot(p1);
    }

    static Plane from_point_normal(Point<T> p, Point<T> n) {
        return {n.x, n.y, n.z, -n.dot(p)};
    }
    static Plane from_points(Point<T> p1, Point<T> p2, Point<T> p3) {
        return Plane(p1, p2, p3);
    }

    Point<T> normal() const { return {a, b, c}; }
    T eval(Point<T> p) const { return a * p.x + b * p.y + c * p.z + d; }
    int side(Point<T> p) const { return sgn(eval(p)); }
};

// 点 p 到平面 pl 的有符号距离（法向量指向侧为正）
template<typename T = ld>
T dist_plane(const Plane<T>& pl, Point<T> p) {
    return pl.eval(p) / pl.normal().len();
}

template<typename T = ld>
Point<T> proj(const Plane<T>& pl, Point<T> p) {
    Point<T> n = pl.normal();
    return p - n * (pl.eval(p) / n.len2());
}

template<typename T = ld>
Point<T> refl(const Plane<T>& pl, Point<T> p) {
    return proj(pl, p) * 2 - p;
}

template<typename T = ld>
T triangle_area(Point<T> a, Point<T> b, Point<T> c) {
    return (b - a).cross(c - a).len() / 2;
}

template<typename T = ld>
Point<T> triangle_normal(Point<T> a, Point<T> b, Point<T> c) {
    Point<T> n = (b - a).cross(c - a);
    return is_zero(n.len2()) ? Point<T>(0, 0, 0) : n.unit();
}

// 求点 p 关于空间三角形 abc 的重心坐标 (alpha, beta, gamma)，退化时置 0
template<typename T = ld>
void barycentric(Point<T> a, Point<T> b, Point<T> c, Point<T> p,
                 T& alpha, T& beta, T& gamma) {
    Point<T> v0 = b - a, v1 = c - a, v2 = p - a;
    T d00 = v0.dot(v0), d01 = v0.dot(v1), d11 = v1.dot(v1);
    T d20 = v2.dot(v0), d21 = v2.dot(v1);
    T den = d00 * d11 - d01 * d01;
    if (is_zero(den)) { alpha = beta = gamma = 0; return; }
    beta = (d11 * d20 - d01 * d21) / den;
    gamma = (d00 * d21 - d01 * d20) / den;
    alpha = 1 - beta - gamma;
}

// 四面体 abcd 的 6 倍有向体积 ((b - a) x (c - a)) . (d - a)
template<typename T = ld>
T tet_volume6(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    return (b - a).cross(c - a).dot(d - a);
}

template<typename T = ld>
T tet_volume(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    return std::abs(tet_volume6(a, b, c, d)) / 6;
}

template<typename T = ld>
Point<T> tet_centroid(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    return (a + b + c + d) / 4;
}

// 判定空间点 p 是否在三角形 abc 内：1 严格内部，0 边/顶点，-1 外部或不共面
template<typename T = ld>
int in_triangle(Point<T> a, Point<T> b, Point<T> c, Point<T> p) {
    if (!is_zero(tet_volume6(a, b, c, p))) return -1;
    T al, be, ga;
    barycentric(a, b, c, p, al, be, ga);
    if (sgn(al) < 0 || sgn(be) < 0 || sgn(ga) < 0) return -1;
    if (is_zero(al) || is_zero(be) || is_zero(ga)) return 0;
    return 1;
}

// 判定空间点 p 是否在四面体 abcd 内：1 严格内部，0 边界面/退化，-1 外部
template<typename T = ld>
int in_tet(Point<T> a, Point<T> b, Point<T> c, Point<T> d, Point<T> p) {
    T v = tet_volume6(a, b, c, d);
    if (is_zero(v)) return 0;
    int sv = sgn(v);
    int sa = sgn(tet_volume6(p, b, c, d)) * sv;
    int sb = sgn(tet_volume6(a, p, c, d)) * sv;
    int sc = sgn(tet_volume6(a, b, p, d)) * sv;
    int sd = sgn(tet_volume6(a, b, c, p)) * sv;
    if (sa < 0 || sb < 0 || sc < 0 || sd < 0) return -1;
    if (sa == 0 || sb == 0 || sc == 0 || sd == 0) return 0;
    return 1;
}

using PointFP = Point<ld>;
using LineFP = Line<ld>;
using SegmentFP = Segment<ld>;
using PlaneFP = Plane<ld>;

} // namespace Geo3D
