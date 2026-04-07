#pragma once
#include "Point.hpp"

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

// 向量旋转，rad为弧度，逆时针为正
template<typename T>
requires std::is_floating_point_v<T>
Point<T> rotate(const Point<T>& p, T rad) {
    return Point<T>(p.x * std::cos(rad) - p.y * std::sin(rad), p.x * std::sin(rad) + p.y * std::cos(rad));
}

// 向量的极角 (-PI, PI]
template<typename T>
requires std::is_floating_point_v<T>
T angle(const Point<T>& p) {
    return std::atan2(p.y, p.x);
}

// 计算从p1向量到p2向量的有向夹角，范围为(-PI, PI]
template<typename T>
requires std::is_floating_point_v<T>
T directed_angle(const Point<T>& p1, const Point<T>& p2) {
    return std::atan2(p1.cross(p2), p1.dot(p2));
}

using PointFP = Point<ld>;

} // namespace Geo2D
