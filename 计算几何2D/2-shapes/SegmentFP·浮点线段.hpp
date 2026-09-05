#pragma once
#include "Segment·线段.hpp"
#include "../1-base/PointFP·浮点二维点.hpp"

/*
 * 浮点线段运算
 *
 * Overview:
 *  提供点与线段、线段与线段之间的距离及线段交点运算。
 *
 * API:
 *  distance_to_segment(p, s) -> T                 — 点到线段的有符号最短距离。O(1)。
 *  segment_intersection_point(s1, s2) -> Point<T> — 两线段交点，要求二者确实相交。O(1)。
 *  segment_segment_distance(s1, s2) -> T          — 两线段间的最短距离。O(1)。
 *  SegmentFP                                      — Segment<ld> 的常用别名。
 *
 * Notes:
 *  仅支持浮点坐标；调用 segment_intersection_point 前应先用 segment_intersect 确认相交。
 *  两线段重合（有无穷多交点）时 segment_intersection_point 的行为未定义。
 *
 * Related:
 *  Segment·线段.hpp: 整数 / 浮点通用的线段类型与相交判定。
 */
namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
T distance_to_segment(Point<T> p, Segment<T> s) {
    if (s.a == s.b) return dist_to(p, s.a);

    Point<T> v1 = s.b - s.a, v2 = p - s.a, v3 = p - s.b;
    if (sgn(v1.dot(v2)) < 0) return len(v2);
    if (sgn(v1.dot(v3)) > 0) return len(v3);

    return std::abs(v1.cross(v2)) / len(v1);
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> segment_intersection_point(Segment<T> s1, Segment<T> s2) {
    Point<T> a = s1.a, b = s1.b, c = s2.a, d = s2.b;
    T t = ((c - a).cross(d - c)) / ((b - a).cross(d - c));
    return a + (b - a) * t;
}

template<typename T>
requires std::is_floating_point_v<T>
T segment_segment_distance(Segment<T> s1, Segment<T> s2) {
    if (segment_intersect(s1, s2)) return (T)0;
    return std::min({
        distance_to_segment(s1.a, s2),
        distance_to_segment(s1.b, s2),
        distance_to_segment(s2.a, s1),
        distance_to_segment(s2.b, s1)
    });
}

using SegmentFP = Segment<ld>;

} // namespace Geo2D
