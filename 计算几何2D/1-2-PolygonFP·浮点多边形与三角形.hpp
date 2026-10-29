#pragma once
#include "1-1-Geo2D·二维浮点基础.hpp"

/*
 * 浮点多边形与三角形特征点
 *
 * Overview:
 *     提供二维浮点多边形的基础度量（有向面积、周长、面积加权质心）、凸性判定、
 *     绕数法点包含判定、有向直线半平面切割多边形，以及三角形外心与内心求解。
 *
 * API:
 *     Polygon<T>             — 多边形顶点序列别名 std::vector<Point<T>>。
 *     poly_area2(poly)       — 两倍有向面积（逆时针为正，顺时针为负，空返 0），
 *                               O(N)。
 *     poly_area(poly)        — 有向面积（逆时针为正，顺时针为负），O(N)。
 *     poly_perimeter(poly)   — 闭合多边形周长（n<2 返 0，n=2 返两倍距离），O(N)。
 *     poly_centroid(poly)    — 面积加权几何质心（空返原点，面积为 0 返首顶点），
 *                               O(N)。
 *     is_convex(poly)        — 判定多边形是否严格凸（支持 CCW/CW，共线或 n<3 返
 *                               false），O(N)。
 *     in_polygon(poly, p)    — 绕数法点包含：边或顶点上返 2，严格内部 1，外部 0，
 *                               O(N)。
 *     polygon_cut(poly, l)   — 用有向直线 l 的左半平面 (l.side(p)>=0) 切割多边形，
 *                               O(N)。
 *     circum_center(a, b, c) — 三角形外接圆圆心（三点共线退化时返回 (NaN, NaN)），
 *                               O(1)。
 *     in_center(a, b, c)     — 三角形内心（角平分线交点，要求三点不共线），O(1)。
 *
 * Notes:
 *     1. in_polygon 对凹凸简单多边形均适用，参数顺序为 (多边形 poly, 查询点 p)。
 *     2. polygon_cut 基于 Sutherland-Hodgman 裁剪，
 *        要求被切部分连通（通常用于凸多边形）。
 *
 * Related:
 *     1-3-RotatingCalipers·凸包与旋转卡尺FP.hpp — 浮点凸包、直径、
 *     最小宽度与最小外接矩形。
 *     1-4-HalfPlane·半平面交.hpp::polygon_kernel — 一般简单多边形的核（O(N log
 *     N)）。
 */

namespace Geo2D {

template<typename T = ld>
using Polygon = std::vector<Point<T>>;
using PolygonFP = Polygon<ld>;

template<typename T = ld>
T poly_area2(const Polygon<T>& poly) {
    T s = 0;
    int n = poly.size();
    rep(i, 0, n - 1) s += poly[i].cross(poly[(i + 1) % n]);
    return s;
}

template<typename T = ld>
T poly_area(const Polygon<T>& poly) {
    return poly_area2(poly) / 2;
}

template<typename T = ld>
T poly_perimeter(const Polygon<T>& poly) {
    int n = poly.size();
    if (n < 2) return 0;
    T s = 0;
    rep(i, 0, n - 1) s += poly[i].dist(poly[(i + 1) % n]);
    return s;
}

template<typename T = ld>
Point<T> poly_centroid(const Polygon<T>& poly) {
    if (poly.empty()) return {0, 0};
    Point<T> c(0, 0);
    T area2 = 0;
    int n = poly.size();
    rep(i, 0, n - 1) {
        int j = (i + 1) % n;
        T cr = poly[i].cross(poly[j]);
        area2 += cr;
        c = c + (poly[i] + poly[j]) * cr;
    }
    if (is_zero(area2)) return poly[0];
    return c / (3 * area2);
}

template<typename T = ld>
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

// 绕数法判定点 p 与多边形关系：2 在边/顶点上，1 严格内部，0 严格外部
template<typename T = ld>
int in_polygon(const Polygon<T>& poly, Point<T> p) {
    int n = poly.size(), wn = 0;
    rep(i, 0, n - 1) {
        int j = (i + 1) % n;
        if (on_segment(Segment<T>(poly[i], poly[j]), p)) return 2;
        int k = ccw(poly[i], poly[j], p);
        int d1 = sgn(poly[i].y - p.y), d2 = sgn(poly[j].y - p.y);
        if (k > 0 && d1 <= 0 && d2 > 0) ++wn;
        if (k < 0 && d2 <= 0 && d1 > 0) --wn;
    }
    return wn ? 1 : 0;
}

// 用有向直线 l 的左侧半平面 (l.side(p) >= 0) 切割凸多边形
template<typename T = ld>
Polygon<T> polygon_cut(const Polygon<T>& poly, const Line<T>& l) {
    Polygon<T> res;
    int n = poly.size();
    rep(i, 0, n - 1) {
        Point<T> cur = poly[i], nxt = poly[(i + 1) % n];
        int s1 = l.side(cur), s2 = l.side(nxt);
        if (s1 >= 0) res.push_back(cur);
        if (s1 * s2 < 0) res.push_back(line_inter(l, Line<T>(cur, nxt)));
    }
    return res;
}

// 三角形外心（三点共线返回 (NaN, NaN)）
template<typename T = ld>
Point<T> circum_center(Point<T> a, Point<T> b, Point<T> c) {
    Point<T> ab = b - a, ac = c - a;
    T d = 2 * ab.cross(ac);
    if (is_zero(d)) {
        T nan = std::numeric_limits<T>::quiet_NaN();
        return {nan, nan};
    }
    T l2_ab = ab.len2(), l2_ac = ac.len2();
    return a + Point<T>((ac.y * l2_ab - ab.y * l2_ac) / d,
                        (ab.x * l2_ac - ac.x * l2_ab) / d);
}

// 三角形内心（角平分线交点）
template<typename T = ld>
Point<T> in_center(Point<T> a, Point<T> b, Point<T> c) {
    T ab = a.dist(b), bc = b.dist(c), ca = c.dist(a);
    return (a * bc + b * ca + c * ab) / (ab + bc + ca);
}

} // namespace Geo2D
