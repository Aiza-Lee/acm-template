#pragma once
#include "../../2-shapes/Polygon·多边形.hpp"

/*
 * 旋转卡尺（整数安全）
 *
 * Overview:
 *  凸多边形上的旋转卡尺技术：求直径（最远点对）。
 *
 * API:
 *  diameter(poly) -> pair<Point<T>, Point<T>> — 凸多边形直径，返回距离最远的一对顶点。Time O(N)。
 *
 * Notes:
 *  模板参数 T: 坐标类型。
 *  输入必须是逆时针严格凸多边形。
 *  width / min_bounding_rect 因结果一般为实数，参见浮点版 RotatingCalipersFP。
 *
 * Related:
 *  Polygon·多边形.hpp::farthest_point_pair: diameter 的实际实现。
 *  RotatingCalipersFP·旋转卡尺FP.hpp::width: 凸多边形最小宽度（浮点专用）。
 *  RotatingCalipersFP·旋转卡尺FP.hpp::min_bounding_rect: 凸多边形最小面积包围矩形（浮点专用）。
 */

namespace Geo2D {

template<typename T>
std::pair<Point<T>, Point<T>> diameter(const Polygon<T>& poly) {
    // 转: Polygon·多边形.hpp::farthest_point_pair
    return farthest_point_pair(poly);
}

} // namespace Geo2D
