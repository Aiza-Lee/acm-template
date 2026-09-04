#pragma once
#include "Point·点.hpp"

/*
 * 浮点二维点
 *
 * Overview:
 *  提供浮点二维点的长度、距离、单位化、旋转和夹角运算。
 *
 * API:
 *  len(p) -> T                           向量长度，使用 std::hypot 避免中间溢出。O(1)。
 *  dist_to(a, b) -> T                    两点距离。O(1)。
 *  normalize(p) -> Point<T>              单位化；零向量触发 AST。O(1)。
 *  rotate(p, rad) -> Point<T>            逆时针旋转 rad 弧度。O(1)。
 *  angle(p) -> T                         向量极角，范围 (-PI, PI]。O(1)。
 *  directed_angle(p1, p2) -> T           从 p1 到 p2 的有向夹角，范围 (-PI, PI]。O(1)。
 *  PointFP                               Point<ld> 的常用别名。
 *
 * Notes:
 *  本文件 API 仅支持浮点类型。
 *
 * Related:
 *  Point·点.hpp: 通用二维点类型与整数安全的基础运算。
 */
namespace Geo2D {

template<typename T>
requires std::is_floating_point_v<T>
T len(const Point<T>& p) {
    return std::hypot(p.x, p.y);
}

template<typename T>
requires std::is_floating_point_v<T>
T dist_to(const Point<T>& a, const Point<T>& b) {
    return len(a - b);
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> normalize(const Point<T>& p) {
    AST(!is_zero(p.len2()));
    T l = len(p);
    return Point<T>(p.x / l, p.y / l);
}

template<typename T>
requires std::is_floating_point_v<T>
Point<T> rotate(const Point<T>& p, T rad) {
    return Point<T>(p.x * std::cos(rad) - p.y * std::sin(rad), p.x * std::sin(rad) + p.y * std::cos(rad));
}

template<typename T>
requires std::is_floating_point_v<T>
T angle(const Point<T>& p) {
    return std::atan2(p.y, p.x);
}

template<typename T>
requires std::is_floating_point_v<T>
T directed_angle(const Point<T>& p1, const Point<T>& p2) {
    return std::atan2(p1.cross(p2), p1.dot(p2));
}

using PointFP = Point<ld>;

} // namespace Geo2D
