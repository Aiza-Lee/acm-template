#pragma once
#include "aizalib.h"

/*
 * 二维整数几何基础
 *
 * Overview:
 *     提供纯整数精确计算几何的底座：整数点/向量、无 atan2 精确极角排序、
 *     线段相交（严格与含端点）、多边形两倍有向面积、纯整数绕数法点包含判定，
 *     以及基于 Pick 定理的边界与内部格点计数。
 *
 * API:
 *     sgn(x)                   — 严格整数符号函数：x < 0 返回 -1，x == 0 返回 0，x
 *                                 > 0 返回 1，O(1)。
 *     Point<T>(x, y)           — 二维整数点/向量，默认 T = i64（大坐标可传 i128）。
 *     Point::dot(p) / cross(p) — 精确点积 x1*x2 + y1*y2 / 叉积 x1*y2 - y1*x2，
 *                                 O(1)。
 *     Point::len2() / dist2(p) — 模长平方 / 两点欧氏距离平方，O(1)。
 *     Point::perp() / half()   — 逆时针转 90 度 (-y, x) / 上半平面 [0, pi) 判别
 *                                 (上半或 x 正半轴返回 0，否则 1)，O(1)。
 *     ccw(p0, p1, p2)          — 转向判定：逆时针 1，顺时针 -1，三点共线 0，O(1)。
 *     polar_less(a, b)         — 纯整数精确极角严格弱序谓词（先按 [0, 2pi) 极角，
 *                                 共线同向按 len2 升序），O(1)。
 *     Segment<T>(a, b)         — 二维整数线段。
 *     on_segment(s, p)         — 点 p 是否在线段 s 上（含两端点），O(1)。
 *     seg_inter_strict(s1, s2) — 两线段是否严格规范相交（交点严格在两线段内部），
 *                                 O(1)。
 *     seg_inter(s1, s2)        — 两线段是否相交（含端点接触与共线重叠），O(1)。
 *     Polygon<T>               — 多边形顶点序列 std::vector<Point<T>>。
 *     poly_area2(poly)         — 多边形两倍有向面积（逆时针为正，恒为整数），O(N)。
 *     is_convex(poly)          — 判定多边形是否严格凸（支持 CCW/CW，无共线角），
 *                                 O(N)。
 *     in_polygon(poly, p)      — 纯整数绕数法点包含：边/顶点上 2，内部 1，外部 0，
 *                                 O(N)。
 *     boundary_points(poly)    — 多边形边界上的整点数量（各边 gcd(|dx|, |dy|)
 *                                 求和），O(N log C)。
 *     interior_points(poly)    — Pick 定理求多边形严格内部整点数量，O(N log C)。
 *
 * Notes:
 *     1. 本文件全程零浮点、零 EPS，所有运算均为精确整数运算。若坐标绝对值达 10^9，
 *        点积/叉积/面积在 i64 下不会溢出；若超过 2*10^9，请实例化 Point<i128>。
 *     2. Pick 定理：多边形面积 S = I + B/2 - 1，故严格内部整点数 I =
 *        (|poly_area2| - B + 2) / 2。
 *
 * Related:
 *     2-2-ConvexInt·整数凸包与二分查询.hpp — 纯整数凸包、旋转卡尺直径平方与 O(log
 *     N) 查询。
 *     2-4-ClosestPair·平面最近点对.hpp     — 分治法 O(N log N) 求平面最近点对。
 *     2-5-PolygonTriangulate·多边形三角剖分.hpp — 简单多边形耳剪法三角剖分。
 */

namespace Geo2DInt {

template<typename T = i64>
inline int sgn(T x) {
    return (x > 0) - (x < 0);
}

template<typename T = i64>
struct Point {
    T x, y;
    Point() : x(0), y(0) {}
    Point(T x, T y) : x(x), y(y) {}

    Point operator+(const Point& p) const { return {x + p.x, y + p.y}; }
    Point operator-(const Point& p) const { return {x - p.x, y - p.y}; }
    Point operator-() const { return {-x, -y}; }
    Point operator*(T k) const { return {x * k, y * k}; }
    Point operator/(T k) const { return {x / k, y / k}; }

    bool operator<(const Point& p) const { return x != p.x ? x < p.x : y < p.y; }
    bool operator==(const Point& p) const { return x == p.x && y == p.y; }
    bool operator!=(const Point& p) const { return !(*this == p); }

    T dot(const Point& p) const { return x * p.x + y * p.y; }
    T cross(const Point& p) const { return x * p.y - y * p.x; }
    T len2() const { return x * x + y * y; }
    T dist2(const Point& p) const { return (*this - p).len2(); }
    Point perp() const { return {-y, x}; }
    // 上半平面 (y > 0) 或 x 正半轴 (y == 0 && x > 0) 返回 0，其余返回 1
    int half() const { return (y > 0 || (y == 0 && x > 0)) ? 0 : 1; }
};

// 转向判定：p0->p1->p2 逆时针返回 1，顺时针返回 -1，共线返回 0
template<typename T = i64>
int ccw(const Point<T>& p0, const Point<T>& p1, const Point<T>& p2) {
    return sgn((p1 - p0).cross(p2 - p0));
}

// 纯整数精确极角排序谓词：按 [0, 2pi) 极角升序，同向按模长平方升序
template<typename T = i64>
bool polar_less(const Point<T>& a, const Point<T>& b) {
    int ha = a.half(), hb = b.half();
    if (ha != hb) return ha < hb;
    T c = a.cross(b);
    return c != 0 ? c > 0 : a.len2() < b.len2();
}

template<typename T = i64>
struct Segment {
    Point<T> a, b;
    Segment() = default;
    Segment(Point<T> a, Point<T> b) : a(a), b(b) {}
};

// 判定点 p 是否在线段 s 上（含两端点）
template<typename T = i64>
bool on_segment(const Segment<T>& s, Point<T> p) {
    return (s.a - p).cross(s.b - p) == 0 && (s.a - p).dot(s.b - p) <= 0;
}

// 线段严格规范相交（交点严格在两线段内部，不含端点）
template<typename T = i64>
bool seg_inter_strict(const Segment<T>& s1, const Segment<T>& s2) {
    return ccw(s1.a, s1.b, s2.a) * ccw(s1.a, s1.b, s2.b) < 0 &&
           ccw(s2.a, s2.b, s1.a) * ccw(s2.a, s2.b, s1.b) < 0;
}

// 线段相交判定（含端点接触与共线重叠）
template<typename T = i64>
bool seg_inter(const Segment<T>& s1, const Segment<T>& s2) {
    const auto &a = s1.a, &b = s1.b, &c = s2.a, &d = s2.b;
    if (std::max(a.x, b.x) < std::min(c.x, d.x) ||
        std::max(c.x, d.x) < std::min(a.x, b.x) ||
        std::max(a.y, b.y) < std::min(c.y, d.y) ||
        std::max(c.y, d.y) < std::min(a.y, b.y)) {
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

template<typename T = i64>
using Polygon = std::vector<Point<T>>;

// 多边形两倍有向面积（逆时针为正，恒为整数）
template<typename T = i64>
T poly_area2(const Polygon<T>& poly) {
    T area = 0;
    int n = poly.size();
    rep(i, 0, n - 1) area += poly[i].cross(poly[(i + 1) % n]);
    return area;
}

// 严格凸多边形判定（相邻边叉积同号且全非零）
template<typename T = i64>
bool is_convex(const Polygon<T>& poly) {
    int n = poly.size();
    if (n < 3) return false;
    bool pos = false, neg = false;
    rep(i, 0, n - 1) {
        int j = (i + 1) % n, k = (j + 1) % n;
        int s = sgn((poly[j] - poly[i]).cross(poly[k] - poly[j]));
        if (s > 0) pos = true;
        if (s < 0) neg = true;
        if (pos && neg) return false;
    }
    return pos ^ neg;
}

// 纯整数绕数法判定点在多边形内：2 边界上，1 严格内部，0 严格外部
template<typename T = i64>
int in_polygon(const Polygon<T>& poly, Point<T> p) {
    int n = poly.size(), wn = 0;
    rep(i, 0, n - 1) {
        int j = (i + 1) % n;
        if (on_segment(Segment<T>(poly[i], poly[j]), p)) return 2;
        int k = ccw(poly[i], poly[j], p);
        T d1 = poly[i].y - p.y, d2 = poly[j].y - p.y;
        if (k > 0 && d1 <= 0 && d2 > 0) ++wn;
        if (k < 0 && d2 <= 0 && d1 > 0) --wn;
    }
    return wn ? 1 : 0;
}

// 整点多边形边界上的格点总数 B = sum gcd(|dx|, |dy|)
template<typename T = i64>
i64 boundary_points(const Polygon<T>& poly) {
    int n = poly.size();
    if (n == 0) return 0;
    if (n == 1) return 1;
    i64 b = 0;
    rep(i, 0, n - 1) {
        int j = (i + 1) % n;
        i64 dx = std::abs((i64)(poly[j].x - poly[i].x));
        i64 dy = std::abs((i64)(poly[j].y - poly[i].y));
        b += std::gcd(dx, dy);
    }
    return n == 2 ? b / 2 + 1 : b;
}

// Pick 定理求整点多边形严格内部格点数 I = (|2S| - B + 2) / 2
template<typename T = i64>
i64 interior_points(const Polygon<T>& poly) {
    if (poly.size() < 3) return 0;
    T a2 = poly_area2(poly);
    if (a2 < 0) a2 = -a2;
    if (a2 == 0) return 0;
    return ((i64)a2 - boundary_points(poly) + 2) / 2;
}

} // namespace Geo2DInt
