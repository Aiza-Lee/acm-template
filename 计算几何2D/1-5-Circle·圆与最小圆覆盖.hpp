#pragma once
#include "1-2-PolygonFP·浮点多边形与三角形.hpp"

/*
 * 圆、切线与最小圆覆盖
 *
 * Overview:
 *     提供二维圆的基础表示、圆与直线求交、两圆求交、过点对圆做切线、
 *     两圆外/内公切线，以及求覆盖点集最小圆的 Welzl 随机增量法。
 *
 * API:
 *     Circle<T>(c, r)                — 构造圆心为 c、半径为 r 的圆。
 *     Circle::contains(p)            — 判定点 p 是否在圆内（含圆周边界），O(1)。
 *     Circle::point(ang)             — 返回圆周上极角 ang 对应的点坐标，O(1)。
 *     circle_line_inter(c, l)        — 圆 c 与直线 l 的交点列表（0/1/2 个，沿 l
 *                                       方向），O(1)。
 *     circle_circle_inter(c1, c2)    — 两圆交点列表（0/1/2 个；相离/内含/
 *                                       完全重合返空），O(1)。
 *     tangents_point_circle(c, p)    — 过点 p 对圆 c 的切点列表（圆内 0 个，圆上返
 *                                       {p}，圆外返 2 个切点），O(1)。
 *     tangents_circle_circle(c1, c2) — 两圆公切线切点对 (p_on_c1, p_on_c2)
 *                                       列表（内含或重合返空；外离 4 条、外切 3 条、
 *                                       相交 2 条、内切 1 条），O(1)。
 *     min_enclosing_circle(pts)      — Welzl 算法求覆盖点集的最小圆（空集返原点零
 *                                       圆），期望时间 O(N)。
 *
 * Notes:
 *     1. 参数顺序统一为 (圆 c, 直线 l / 点 p)，如 circle_line_inter(c, l)、
 *        tangents_point_circle(c, p)。
 *     2. tangents_circle_circle 返回的每对切点 {q1, q2} 分别位于 c1 和 c2 上，
 *        直接构造 Line(q1, q2) 即可得到公切线。
 *
 * Related:
 *     1-6-CircleArea·圆面积交并.hpp — 圆与多边形面积交、多圆面积并。
 *     1-7-Inversion·反演变换.hpp — 点、圆、直线的反演变换。
 */

namespace Geo2D {

template<typename T = ld>
struct Circle {
    Point<T> c;
    T r;
    Circle() : c(0, 0), r(0) {}
    Circle(Point<T> c, T r) : c(c), r(r) {}

    bool contains(Point<T> p) const { return cmp(c.dist2(p), r * r) <= 0; }
    Point<T> point(T ang) const {
        return {c.x + r * std::cos(ang), c.y + r * std::sin(ang)};
    }
};

using CircleFP = Circle<ld>;

// 圆与直线交点：0 个（相离）、1 个（相切）、2 个（相交，按直线方向排列）
template<typename T = ld>
std::vector<Point<T>> circle_line_inter(Circle<T> c, const Line<T>& l) {
    Point<T> p0 = proj(l, c.c);
    T d = c.c.dist(p0);
    int s = sgn(d - c.r);
    if (s > 0) return {};
    if (s == 0) return {p0};
    Point<T> off = l.dir() * safe_sqrt(c.r * c.r - d * d);
    return {p0 - off, p0 + off};
}

// 两圆交点：相离、内含或完全重合时返回空列表
template<typename T = ld>
std::vector<Point<T>> circle_circle_inter(Circle<T> c1, Circle<T> c2) {
    T d = c1.c.dist(c2.c);
    if (sgn(d - c1.r - c2.r) > 0 || sgn(d - std::abs(c1.r - c2.r)) < 0) return {};
    if (is_zero(d) && is_zero(c1.r - c2.r)) return {};
    T a = (c1.r * c1.r - c2.r * c2.r + d * d) / (2 * d);
    T h = safe_sqrt(c1.r * c1.r - a * a);
    Point<T> v = (c2.c - c1.c).unit();
    Point<T> mid = c1.c + v * a, off = v.perp() * h;
    if (is_zero(h)) return {mid};
    return {mid + off, mid - off};
}

// 过点 p 对圆 c 做切线，返回圆 c 上的切点列表（0 / 1 / 2 个）
template<typename T = ld>
std::vector<Point<T>> tangents_point_circle(Circle<T> c, Point<T> p) {
    T d = c.c.dist(p);
    int s = sgn(d - c.r);
    if (s < 0) return {};
    if (s == 0) return {p};
    T ang = safe_acos(c.r / d);
    Point<T> v = (p - c.c).unit() * c.r;
    return {c.c + v.rot(ang), c.c + v.rot(-ang)};
}

// 两圆公切线：返回切点对 {c1 上的切点, c2 上的切点} 列表（先外公切线，后内公切线）
template<typename T = ld>
std::vector<std::pair<Point<T>, Point<T>>> tangents_circle_circle(
    Circle<T> c1, Circle<T> c2
) {
    std::vector<std::pair<Point<T>, Point<T>>> res;
    if (c1.r < c2.r) {
        auto tmp = tangents_circle_circle(c2, c1);
        for (auto& [p1, p2] : tmp) res.push_back({p2, p1});
        return res;
    }
    T d = c1.c.dist(c2.c), rdiff = c1.r - c2.r, rsum = c1.r + c2.r;
    if (sgn(d - rdiff) < 0 || (is_zero(d) && is_zero(rdiff))) return res;
    Point<T> v = (c2.c - c1.c).unit();
    if (sgn(d - rdiff) >= 0) {
        T ang = safe_acos(rdiff / d);
        Point<T> v1 = v.rot(ang), v2 = v.rot(-ang);
        res.push_back({c1.c + v1 * c1.r, c2.c + v1 * c2.r});
        if (!is_zero(ang)) res.push_back({c1.c + v2 * c1.r, c2.c + v2 * c2.r});
    }
    if (sgn(d - rsum) >= 0) {
        T ang = safe_acos(rsum / d);
        Point<T> v1 = v.rot(ang), v2 = v.rot(-ang);
        res.push_back({c1.c + v1 * c1.r, c2.c - v1 * c2.r});
        if (!is_zero(ang)) res.push_back({c1.c + v2 * c1.r, c2.c - v2 * c2.r});
    }
    return res;
}

// Welzl 随机增量法求覆盖 pts 所有点的最小圆，期望时间 O(N)
template<typename T = ld>
Circle<T> min_enclosing_circle(std::vector<Point<T>> pts) {
    if (pts.empty()) return {{0, 0}, 0};
    if (pts.size() == 1) return {pts[0], 0};
    std::mt19937 rng(std::random_device{}());
    std::shuffle(pts.begin(), pts.end(), rng);
    Circle<T> c(pts[0], 0);
    int n = pts.size();
    rep(i, 1, n - 1) if (!c.contains(pts[i])) {
        c = {pts[i], 0};
        rep(j, 0, i - 1) if (!c.contains(pts[j])) {
            c = {(pts[i] + pts[j]) / 2, pts[i].dist(pts[j]) / 2};
            rep(k, 0, j - 1) if (!c.contains(pts[k])) {
                Point<T> o = circum_center(pts[i], pts[j], pts[k]);
                c = {o, o.dist(pts[i])};
            }
        }
    }
    return c;
}

} // namespace Geo2D
