#pragma once
#include "aizalib.h"

/*
 * 二维浮点几何基础 (Point / Line / Segment)
 *
 * Overview:
 *     提供二维浮点计算几何的基础图元：容差比较、安全初等函数、二维向量/点、
 *     以及统一采用两点式 {a, b} 表示的有向直线与线段（投影、对称、相交、距离）。
 *
 * API:
 *     EPS / PI                        — 全局浮点容差 (1e-10) / 圆周率常数。
 *     sgn(x) / cmp(x, y) / is_zero(x) — 容差符号 (-1/0/1)、三路比较与判零，O(1)。
 *     clamp_unit(x) / safe_sqrt(x)    — 截断至 [-1, 1] / 负数按 0 开方，O(1)。
 *     safe_acos(x) / safe_asin(x)     — 钳制输入至 [-1, 1] 后求反三角函数，O(1)。
 *     Point<T>(x, y)                  — 二维浮点点/向量，默认 T = ld。
 *     Point::dot(p) / cross(p)        — 点积 / 叉积，O(1)。
 *     Point::len2() / len()           — 模长平方 / 模长 (std::hypot 防溢出)，O(1)。
 *     Point::dist2(p) / dist(p)       — 两点距离平方 / 欧氏距离，O(1)。
 *     Point::unit()                   — 单位化向量（零向量触发 AST），O(1)。
 *     Point::perp() / rot(rad)        — 逆时针转 90 度 (-y, x) / 逆时针旋转 rad
 *                                        弧度。
 *     Point::ang() / half()           — 极角 (-PI, PI] / 上半平面 [0,PI) 返 0 否则
 *                                        1。
 *     ccw(p0, p1, p2)                 — p0->p1->p2 逆时针返 1，顺时针 -1，共线 0。
 *     polar_less(a, b)                — 极角排序严格弱序谓词（先半平面再叉积再模
 *                                        长）。
 *     directed_angle(a, b)            — 从向量 a 到 b 的有向夹角，范围 (-PI, PI]。
 *     Line<T>(a, b) / Segment<T>      — 两点式有向直线 (a -> b) 与线段类型别名。
 *     Line::from_point_dir(p, dir)    — 由直线上一点 p 与方向向量 dir
 *                                        构造有向直线。
 *     Line::from_point_normal(p, n)   — 由一点 p 与指向左侧的法向量 n
 *                                        构造有向直线。
 *     Line::from_general(A, B, C)     — 由一般式 Ax + By + C = 0 构造（法向 (A,B)
 *                                        朝左）。
 *     Line::v() / dir() / normal()    — 方向向量 b-a / 单位方向 /
 *                                        指向左侧的单位法向。
 *     Line::eval(p) / side(p)         — 叉积 (b-a)x(p-a) / 左侧 1、右侧 -1、线上
 *                                        0。
 *     parallel(l1, l2) / same_dir     — 判定两直线平行 / 平行且同向，O(1)。
 *     line_inter(l1, l2)              — 两非平行直线的唯一交点（平行触发 AST），
 *                                        O(1)。
 *     proj(l, p) / refl(l, p)         — 点 p 在直线 l 上的正交投影 / 对称点，O(1)。
 *     dist_line(l, p)                 — 点 p 到直线 l 的无符号距离，O(1)。
 *     on_segment(s, p)                — 点 p 是否在线段 s 上（含两端点），O(1)。
 *     seg_inter(s1, s2)               — 两线段是否相交（含端点接触与共线重叠），
 *                                        O(1)。
 *     dist_seg(s, p)                  — 点 p 到线段 s 的最短距离，O(1)。
 *     seg_dist(s1, s2)                — 两线段 s1 与 s2 之间的最短距离，O(1)。
 *
 * Notes:
 *     1. 本文件专为浮点几何设计（默认 T = ld）；纯整数几何请使用 Geo2DInt。
 *     2. 有向直线 Line(a, b) 方向为 a -> b，左手侧为正半平面 (side(p) == 1)。
 *     3. 两相交且不共线的线段求交点可直接调用 line_inter(s1, s2)。
 *     4. 参数顺序统一遵循 (线/线段, 查询点)，如 proj(l, p)、on_segment(s, p)。
 *
 * Related:
 *     1-2-PolygonFP·浮点多边形与三角形.hpp — 多边形度量、包含、切割与三角形四心。
 *     2-1-Geo2DInt·二维整数基础.hpp — 纯整数精确几何底座。
 */

namespace Geo2D {

inline constexpr ld EPS = 1e-10;
inline constexpr ld PI = std::numbers::pi_v<ld>;

template<typename T = ld>
inline int sgn(T x) { return x < -EPS ? -1 : (x > EPS ? 1 : 0); }

template<typename T = ld>
inline int cmp(T x, T y) { return sgn(x - y); }

template<typename T = ld>
inline bool is_zero(T x) { return sgn(x) == 0; }

template<typename T = ld>
inline T clamp_unit(T x) { return std::clamp(x, (T)-1, (T)1); }

template<typename T = ld>
inline T safe_sqrt(T x) { return std::sqrt(std::max((T)0, x)); }

template<typename T = ld>
inline T safe_acos(T x) { return std::acos(clamp_unit(x)); }

template<typename T = ld>
inline T safe_asin(T x) { return std::asin(clamp_unit(x)); }

template<typename T = ld>
struct Point {
    T x, y;
    Point() : x(0), y(0) {}
    Point(T x, T y) : x(x), y(y) {}

    Point operator+(const Point& p) const { return {x + p.x, y + p.y}; }
    Point operator-(const Point& p) const { return {x - p.x, y - p.y}; }
    Point operator-() const { return {-x, -y}; }
    Point operator*(T k) const { return {x * k, y * k}; }
    Point operator/(T k) const { return {x / k, y / k}; }

    bool operator<(const Point& p) const {
        int cx = cmp(x, p.x);
        return cx != 0 ? cx < 0 : cmp(y, p.y) < 0;
    }
    bool operator==(const Point& p) const {
        return cmp(x, p.x) == 0 && cmp(y, p.y) == 0;
    }
    bool operator!=(const Point& p) const { return !(*this == p); }

    T dot(const Point& p) const { return x * p.x + y * p.y; }
    T cross(const Point& p) const { return x * p.y - y * p.x; }
    T len2() const { return x * x + y * y; }
    T len() const { return std::hypot(x, y); }
    T dist2(const Point& p) const { return (*this - p).len2(); }
    T dist(const Point& p) const { return (*this - p).len(); }

    Point unit() const {
        T l = len();
        AST(!is_zero(l));
        return {x / l, y / l};
    }
    Point perp() const { return {-y, x}; }
    Point rot(T rad) const {
        T c = std::cos(rad), s = std::sin(rad);
        return {x * c - y * s, x * s + y * c};
    }
    T ang() const { return std::atan2(y, x); }
    // 极角所在半平面：[0, PI) 返回 0，[PI, 2PI) 及原点返回 1
    int half() const {
        int sy = sgn(y), sx = sgn(x);
        return (sy > 0 || (sy == 0 && sx > 0)) ? 0 : 1;
    }
};

// 旋转方向：逆时针(左转) 1，顺时针(右转) -1，共线 0
template<typename T = ld>
int ccw(const Point<T>& p0, const Point<T>& p1, const Point<T>& p2) {
    return sgn((p1 - p0).cross(p2 - p0));
}

// 极角排序谓词（可直接传给 std::sort）
template<typename T = ld>
bool polar_less(const Point<T>& a, const Point<T>& b) {
    int ha = a.half(), hb = b.half();
    if (ha != hb) return ha < hb;
    int c = sgn(a.cross(b));
    return c != 0 ? c > 0 : cmp(a.len2(), b.len2()) < 0;
}

template<typename T = ld>
T directed_angle(const Point<T>& a, const Point<T>& b) {
    return std::atan2(a.cross(b), a.dot(b));
}

// 两点式有向直线 a -> b（左侧为正半平面）
template<typename T = ld>
struct Line {
    Point<T> a, b;
    Line() = default;
    Line(Point<T> a, Point<T> b) : a(a), b(b) {}

    static Line from_point_dir(Point<T> p, Point<T> dir) { return {p, p + dir}; }
    static Line from_point_normal(Point<T> p, Point<T> n) {
        return {p, p + Point<T>(n.y, -n.x)};
    }
    static Line from_general(T A, T B, T C) {
        Point<T> p0 = std::abs(A) > std::abs(B) ?
                      Point<T>(-C / A, 0) : Point<T>(0, -C / B);
        return {p0, p0 + Point<T>(B, -A)};
    }

    Point<T> v() const { return b - a; }
    T len2() const { return a.dist2(b); }
    T len() const { return a.dist(b); }
    Point<T> dir() const { return v().unit(); }
    Point<T> normal() const { return v().perp().unit(); }
    T eval(Point<T> p) const { return v().cross(p - a); }
    int side(Point<T> p) const { return sgn(eval(p)); } // 1 左, -1 右, 0 线上
};

template<typename T = ld>
using Segment = Line<T>;

template<typename T = ld>
bool parallel(const Line<T>& l1, const Line<T>& l2) {
    return sgn(l1.v().cross(l2.v())) == 0;
}

template<typename T = ld>
bool same_dir(const Line<T>& l1, const Line<T>& l2) {
    return parallel(l1, l2) && sgn(l1.v().dot(l2.v())) > 0;
}

template<typename T = ld>
Point<T> line_inter(const Line<T>& l1, const Line<T>& l2) {
    Point<T> v1 = l1.v(), v2 = l2.v();
    T det = v1.cross(v2);
    AST(!is_zero(det));
    return l1.a + v1 * ((l2.a - l1.a).cross(v2) / det);
}

template<typename T = ld>
Point<T> proj(const Line<T>& l, Point<T> p) {
    Point<T> v = l.v();
    return l.a + v * ((p - l.a).dot(v) / v.len2());
}

template<typename T = ld>
Point<T> refl(const Line<T>& l, Point<T> p) {
    return proj(l, p) * 2 - p;
}

template<typename T = ld>
T dist_line(const Line<T>& l, Point<T> p) {
    return std::abs(l.eval(p)) / l.len();
}

template<typename T = ld>
bool on_segment(const Segment<T>& s, Point<T> p) {
    return sgn((s.a - p).cross(s.b - p)) == 0 &&
           sgn((s.a - p).dot(s.b - p)) <= 0;
}

template<typename T = ld>
bool seg_inter(const Segment<T>& s1, const Segment<T>& s2) {
    const auto &a = s1.a, &b = s1.b, &c = s2.a, &d = s2.b;
    if (cmp(std::max(a.x, b.x), std::min(c.x, d.x)) < 0 ||
        cmp(std::max(c.x, d.x), std::min(a.x, b.x)) < 0 ||
        cmp(std::max(a.y, b.y), std::min(c.y, d.y)) < 0 ||
        cmp(std::max(c.y, d.y), std::min(a.y, b.y)) < 0) {
        return false;
    }
    int d1 = ccw(a, b, c), d2 = ccw(a, b, d);
    int d3 = ccw(c, d, a), d4 = ccw(c, d, b);
    if (d1 * d2 < 0 && d3 * d4 < 0) return true;
    if (d1 == 0 && on_segment(s1, c)) return true;
    if (d2 == 0 && on_segment(s1, d)) return true;
    if (d3 == 0 && on_segment(s2, a)) return true;
    if (d4 == 0 && on_segment(s2, b)) return true;
    return false;
}

template<typename T = ld>
T dist_seg(const Segment<T>& s, Point<T> p) {
    if (s.a == s.b) return p.dist(s.a);
    Point<T> v1 = s.v(), v2 = p - s.a, v3 = p - s.b;
    if (sgn(v1.dot(v2)) < 0) return v2.len();
    if (sgn(v1.dot(v3)) > 0) return v3.len();
    return std::abs(v1.cross(v2)) / v1.len();
}

template<typename T = ld>
T seg_dist(const Segment<T>& s1, const Segment<T>& s2) {
    if (seg_inter(s1, s2)) return 0;
    return std::min({
        dist_seg(s2, s1.a), dist_seg(s2, s1.b),
        dist_seg(s1, s2.a), dist_seg(s1, s2.b)
    });
}

using PointFP = Point<ld>;
using LineFP = Line<ld>;
using SegmentFP = Segment<ld>;

} // namespace Geo2D
