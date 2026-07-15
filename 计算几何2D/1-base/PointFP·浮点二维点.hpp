#pragma once
#include "Point·点.hpp"

namespace Geo2D {

/** @brief |p|;走 std::hypot 避免中间溢出。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
T len(const Point<T>& p) {
    return std::hypot(p.x, p.y);
}

/** @brief |a - b|;浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
T dist_to(const Point<T>& a, const Point<T>& b) {
    return len(a - b);
}

/** @brief 单位化向量。浮点专用,O(1)。退化:零向量会触发 assert。 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> normalize(const Point<T>& p) {
	AST(!is_zero(p.len2()));
	T l = len(p);
	return Point<T>(p.x / l, p.y / l);
}

/** @brief 向量旋转 rad 弧度,逆时针为正。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
Point<T> rotate(const Point<T>& p, T rad) {
    return Point<T>(p.x * std::cos(rad) - p.y * std::sin(rad), p.x * std::sin(rad) + p.y * std::cos(rad));
}

/** @brief 向量的极角,范围 (-π, π]。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
T angle(const Point<T>& p) {
    return std::atan2(p.y, p.x);
}

/** @brief 从 p1 到 p2 的有向夹角,范围 (-π, π]。浮点专用,O(1)。 */
template<typename T>
requires std::is_floating_point_v<T>
T directed_angle(const Point<T>& p1, const Point<T>& p2) {
    return std::atan2(p1.cross(p2), p1.dot(p2));
}

/** @brief 常用浮点实例别名。 */
using PointFP = Point<ld>;

} // namespace Geo2D
