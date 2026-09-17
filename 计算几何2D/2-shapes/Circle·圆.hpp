#pragma once
#include "../1-base/Line·直线.hpp"

/*
 * 圆 (Circle)
 *
 * Overview:
 *     提供圆的几何表达、点与圆位置关系判定，以及圆与直线、圆与圆相交和公切线求解。
 *     - 内部结构：包含圆心点 c 与半径 r。
 *     - 几何计算：
 *       1. 点在圆内：dist2(c, p) <= r^2。
 *       2. 圆线相交：通过圆心到直线垂足与弦长勾股定理求交点。
 *       3. 圆圆相交：余弦定理求两圆心连线与公弦夹角，旋转偏移得交点。
 *       4. 公切线：包括外公切线与内公切线，向量旋转求得切点。
 *     - 工具：Circle、contains、point、circle_line_intersection、
 *       circle_circle_intersection、tangents_point_circle、tangents_circle_circle。
 *
 * API:
 *     Circle<T>(c, r)                    — 构造指定圆心与半径的圆。
 *     contains(p)                        — 判定点 p 是否在圆内（含边界），O(1)。
 *     point(angle)                       — 返回圆周上对应极角 angle 的点坐标，
 *                                           O(1)。
 *     circle_line_intersection(c, l)     — 求圆与直线的交点列表（0/1/2 个），O(1)。
 *     circle_circle_intersection(c1, c2) — 求两圆的交点列表（0/1/2 个），O(1)。
 *     tangents_point_circle(p, c)        — 求过点 p 的切点列表（0/1/2 个），O(1)。
 *     tangents_circle_circle(c1, c2)     — 求两圆的公切线切点对列表，O(1)。
 *     CircleFP                           — Circle<ld> 的常用别名。
 *
 * Notes:
 *     1. 仅支持浮点类型 T（要求 std::is_floating_point_v<T>）。
 *     2. circle_circle_intersection 不处理两圆完全重合情形（返回空集）。
 *     3. tangents_point_circle: 点在圆内返回空；点在圆上返回 1 个切点（p 自身）；
 *        点在圆外返回 2 个切点。
 *     4. tangents_circle_circle: 内含或完全重合时返回空；
 *        否则分别给出外公切线与内公切线对应的切点对。
 *
 * Related:
 *     MinEnclosingCircle·最小圆覆盖.hpp::smallest_enclosing_circle:
 *     给定点集的最小覆盖圆（Welzl）。
 *     CirclePolygonArea·圆与多边形面积交.hpp::circle_polygon_area:
 *     圆与多边形有向面积交。
 */
namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
struct Circle {
    Point<T> c; T r;
    Circle() : c(Point<T>()), r(0) {}
    Circle(Point<T> c, T r) : c(c), r(r) {}
    bool contains(Point<T> p) const { return cmp(c.dist2(p), r * r) <= 0; }
    Point<T> point(T angle) const {
        return Point<T>(c.x + r * std::cos(angle), c.y + r * std::sin(angle));
    }
};

template<typename T>
requires std::is_floating_point_v<T>
std::vector<Point<T>> circle_line_intersection(Circle<T> c, Line<T> l) {
    std::vector<Point<T>> res;
    T d = distance_to_point(l, c.c);

    if (sgn(d - c.r) > 0) return res;  // 无交点

    // 直线垂足
    Point<T> n(l.a, l.b);
    Point<T> p0 = c.c - n * (l.eval(c.c) / n.len2());

    if (sgn(d - c.r) == 0) {  // 一个交点（切线）
        res.push_back(p0);
        return res;
    }

    // 两个交点
    T len2 = safe_sqrt(c.r * c.r - d * d);
    Point<T> dir = direction(l) * len2;
    res.push_back(p0 + dir);
    res.push_back(p0 - dir);
    return res;
}

template<typename T>
requires std::is_floating_point_v<T>
std::vector<Point<T>> circle_circle_intersection(Circle<T> c1, Circle<T> c2) {
    std::vector<Point<T>> res;
    T d = dist_to(c1.c, c2.c);

    if (sgn(d - c1.r - c2.r) > 0 || sgn(d - std::abs(c1.r - c2.r)) < 0)
        return res; // 相离或内含

    if (sgn(d) == 0 && sgn(c1.r - c2.r) == 0)
        return res; // 完全重合

    T a = (c1.r * c1.r - c2.r * c2.r + d * d) / (2 * d);
    T h = safe_sqrt(c1.r * c1.r - a * a);

    Point<T> v = normalize(c2.c - c1.c);
    Point<T> mid = c1.c + v * a;
    Point<T> normal = rotate(v, (T)(PI/2)) * h;

    res.push_back(mid + normal);
    if (sgn(h) != 0) res.push_back(mid - normal);
    return res;
}

template<typename T>
requires std::is_floating_point_v<T>
std::vector<Point<T>> tangents_point_circle(Point<T> p, Circle<T> c) {
    std::vector<Point<T>> res;
    T d = dist_to(c.c, p);
    int s = sgn(d - c.r);
    if (s < 0) return res; // 点在圆内
    if (s == 0) { // 点在圆上
        res.push_back(p);
        return res;
    }

    T angle_val = safe_acos(c.r / d);
    Point<T> v = normalize(p - c.c) * c.r;
    res.push_back(c.c + rotate(v, angle_val));
    res.push_back(c.c + rotate(v, -angle_val));
    return res;
}

template<typename T>
requires std::is_floating_point_v<T>
std::vector<std::pair<Point<T>, Point<T>>> tangents_circle_circle(
    Circle<T> c1, Circle<T> c2
) {
    std::vector<std::pair<Point<T>, Point<T>>> res;
    if (c1.r < c2.r) {
        auto tmp = tangents_circle_circle(c2, c1);
        for (auto p : tmp) res.push_back({p.second, p.first});
        return res;
    }

    T d = dist_to(c1.c, c2.c);
    T rdiff = c1.r - c2.r;
    T rsum = c1.r + c2.r;

    if (sgn(d - rdiff) < 0) return res; // 内含

    if (sgn(d) == 0 && sgn(c1.r - c2.r) == 0) return res; // 重合，无数条切线

    Point<T> v = normalize(c2.c - c1.c);

    // 外公切线
    if (sgn(d - rdiff) >= 0) {
        T angle_val = safe_acos(rdiff / d);
        Point<T> v1 = rotate(v, angle_val);
        Point<T> v2 = rotate(v, -angle_val);
        res.push_back({c1.c + v1 * c1.r, c2.c + v1 * c2.r});
        if (sgn(angle_val) != 0)
            res.push_back({c1.c + v2 * c1.r, c2.c + v2 * c2.r});
    }

    // 内公切线
    if (sgn(d - rsum) >= 0) {
        T angle_val = safe_acos(rsum / d);
        Point<T> v1 = rotate(v, angle_val);
        Point<T> v2 = rotate(v, -angle_val);
        res.push_back({c1.c + v1 * c1.r, c2.c + v1 * (-c2.r)}); // 注意内切线方向相反
        if (sgn(angle_val) != 0)
            res.push_back({c1.c + v2 * c1.r, c2.c + v2 * (-c2.r)});
    }

    return res;
}

using CircleFP = Circle<ld>;

} // namespace Geo2D
