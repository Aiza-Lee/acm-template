#pragma once
#include "../../2-shapes/Polygon.hpp"

/**
 * [RotatingCalipers (旋转卡尺)]
 * 算法介绍: 凸多边形上的旋转卡尺技术（整数安全部分）
 * 模板参数: T (坐标类型)
 * Interface:
 *   - std::pair<Point<T>, Point<T>> diameter(const Polygon<T>& poly)
 * Note:
 * 1. Time: O(N)
 * 2. 输入必须是逆时针严格凸多边形。
 * 3. diameter 内部转调 Polygon::farthest_point_pair，返回距离最远的一对顶点。
 * 4. width 与 min_bounding_rect 因结果一般为实数，放在 float-oriented/RotatingCalipersFP.hpp 中。
 */

namespace Geo2D {

template<typename T>
std::pair<Point<T>, Point<T>> diameter(const Polygon<T>& poly) {
	return farthest_point_pair(poly);
}

} // namespace Geo2D
