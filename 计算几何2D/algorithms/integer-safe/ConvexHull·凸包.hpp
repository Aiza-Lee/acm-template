#pragma once
#include "../../2-shapes/Polygon·多边形.hpp"

/*
 * 凸包
 *
 * Overview:
 *  Andrew 单调链算法求二维点集凸包，以及凸包的几何重心。
 *
 * API:
 *  convex_hull(pts) -> Polygon<T>          — 点集凸包；结果按逆时针排列，自动去重且不包含多余共线点。Time O(N log N), Space O(N)。
 *  convex_hull_centroid(hull) -> Point<ld> — 凸包（逆时针）几何重心；空凸包返回 (0, 0)。Time O(N)。
 *
 * Notes:
 *  模板参数 T: 点集 / 凸包顶点坐标类型。
 *  若需保留共线点，应将 cross 判断改为 < 0（而非 <= 0）。
 *  convex_hull_centroid 返回 Point<ld>，与输入 T 无关；底层转调 polygon_centroid；area=0 时由其兜底返回首顶点。
 *
 * Related:
 *  Polygon·多边形.hpp::polygon_centroid: 多边形层通用质心，area=0 时返回首顶点。
 */

namespace Geo2D {

template<typename T>
Polygon<T> convex_hull(std::vector<Point<T>> pts) {
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    int n = pts.size(), k = 0;
    if (n <= 1) return pts;

    Polygon<T> hull(2 * n);
    // 下凸包
    rep(i, 0, n - 1) {
        while (k > 1 && sgn((hull[k-1] - hull[k-2]).cross(pts[i] - hull[k-1])) <= 0) --k;
        hull[k++] = pts[i];
    }
    // 上凸包
    int t = k;
    per(i, n - 2, 0) {
        while (k > t && sgn((hull[k-1] - hull[k-2]).cross(pts[i] - hull[k-1])) <= 0) --k;
        hull[k++] = pts[i];
    }
    hull.resize(k - 1);
    return hull;
}

template<typename T>
Point<ld> convex_hull_centroid(const Polygon<T>& hull) {
    if (hull.empty()) return Point<ld>(0, 0);   // 本地兜底:空凸包返回原点
    // 转: Polygon·多边形.hpp::polygon_centroid
    return polygon_centroid(hull);
}

} // namespace Geo2D
