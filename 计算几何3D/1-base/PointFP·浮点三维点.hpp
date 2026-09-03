#pragma once
#include "Point·三维点.hpp"

/*
 * 浮点三维点
 *
 * Overview:
 *  提供浮点三维点的长度、距离、单位化与角度运算。
 *
 * API:
 *  len(p) -> T: 向量长度,使用 std::hypot 避免中间溢出。O(1)。
 *  dist_to(a, b) -> T: 两点距离。O(1)。
 *  normalize(p) -> Point<T>: 单位化;零向量触发 AST。O(1)。
 *  angle_xy(p) -> T: 向量 XY 平面投影的极角,范围 (-PI, PI]。O(1)。
 *  angle_between(a, b) -> T: 两向量夹角（无符号, [0, PI]）。O(1)。
 *  PointFP: Point<ld> 的常用别名。
 *
 * Notes:
 *  本文件 API 仅支持浮点类型。
 *  normalize 对零向量 AST(!is_zero(p.len2())): 调用方应先用 is_zero 守门。
 *  三维空间的有向夹角依赖旋转轴;若需要带符号版本,调用方结合 cross 自行 atan2。
 *
 * Related:
 *  Point·三维点.hpp: 通用三维点类型与整数安全的基础运算。
 */
namespace Geo3D {

template<typename T>
requires std::is_floating_point_v<T>
T len(const Point<T>& p) {
    return std::hypot(p.x, p.y, p.z);  // std::hypot 三参重载避免平方溢出
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
    return Point<T>(p.x / l, p.y / l, p.z / l);
}

template<typename T>
requires std::is_floating_point_v<T>
T angle_xy(const Point<T>& p) {
    return std::atan2(p.y, p.x);   // XY 平面投影的极角; z 仅作为向量分量
}

template<typename T>
requires std::is_floating_point_v<T>
T angle_between(const Point<T>& a, const Point<T>& b) {
    T cos_theta = clamp_unit(a.dot(b) / (len(a) * len(b)));
    return std::acos(cos_theta);   // 无符号, 范围 [0, PI]
}

using PointFP = Point<ld>;

} // namespace Geo3D
