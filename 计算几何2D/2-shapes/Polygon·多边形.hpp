#pragma once
#include "Segment·线段.hpp"
#include "../1-base/Line·直线.hpp"
#include "../algorithms/float-oriented/HalfPlane·半平面交.hpp"
namespace Geo2D {

template<typename T>
using Polygon = std::vector<Point<T>>;

/*
 * 多边形
 *
 * Overview:
 *  提供多边形基本度量、点包含判定、凸性、重心、多边形核和最远点对。
 *
 * API:
 *  Polygon<T>                                              — 点序列别名，等价于 std::vector<Point<T>>。
 *  point_in_polygon(p, poly) -> int                        — 射线法 / winding number；边上 2，内部 1，外部 0。O(N)。
 *  polygon_area_2(poly) -> T                               — 两倍有向面积；逆时针为正，顺时针为负；空多边形返回 0。O(N)。
 *  polygon_area(poly) -> ld                                — 有向面积，逆时针为正，顺时针为负。O(N)。
 *  polygon_perimeter(poly) -> ld                           — 闭合路径周长；空 / 单点返回 0；两点返回 2×距离。O(N)。
 *  is_convex(poly) -> bool                                 — 顶点数 < 3 或全共线返回 false；CCW 严格凸返回 true。O(N)。
 *  polygon_centroid(poly) -> Point<ld>                     — 面积加权质心；area=0 时返回首顶点。O(N)。
 *  polygon_kernel(poly) -> Polygon<T>                      — 多边形核（仅支持浮点坐标）。O(N log N)。
 *  farthest_point_pair(convex) -> pair<Point<T>, Point<T>> — 凸多边形直径；n<=1 返回 ({0,0}, {0,0})，n=2 返回两端点。O(N)。
 *
 * Notes:
 *  面积符号约定：逆时针为正、顺时针为负。
 *  polygon_kernel 要求输入为逆时针；CW 输入会得到错误方向的半平面。退化为空 / 单点 / 两点 / 全共线多边形时核为空。
 *
 * Related:
 *  HalfPlane·半平面交.hpp::half_plane_intersection: polygon_kernel 的底层实现。
 *  ConvexHull·凸包.hpp::convex_hull_centroid: 凸包质心入口，空凸包兜底 (0,0)。
 *  RotatingCalipers·旋转卡尺.hpp::diameter: farthest_point_pair 的算法层 wrapper。
 */

template<typename T>
int point_in_polygon(Point<T> p, const Polygon<T>& poly) {
    int n = poly.size(), wn = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        Segment<T> s(poly[i], poly[j]);
        if (on_segment(p, s)) return 2;

        int k = sgn((poly[j] - poly[i]).cross(p - poly[i]));
        int d1 = sgn(poly[i].y - p.y);
        int d2 = sgn(poly[j].y - p.y);

        if (k > 0 && d1 <= 0 && d2 > 0) wn++;
        if (k < 0 && d2 <= 0 && d1 > 0) wn--;
    }
    return wn ? 1 : 0;
}

template<typename T>
T polygon_area_2(const Polygon<T>& poly) {
    T area = 0;
    int n = poly.size();
    for (int i = 0; i < n; i++) {
        area += poly[i].cross(poly[(i + 1) % n]);
    }
    return area;
}

template<typename T>
ld polygon_area(const Polygon<T>& poly) {
    return static_cast<ld>(polygon_area_2(poly)) / 2.0;
}

template<typename T>
ld polygon_perimeter(const Polygon<T>& poly) {
    int n = poly.size();
    if (n < 2) return (ld)0;
    ld s = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        s += std::hypot((ld)(poly[j].x - poly[i].x), (ld)(poly[j].y - poly[i].y));
    }
    return s;
}

template<typename T>
bool is_convex(const Polygon<T>& poly) {
    int n = poly.size();
    if (n < 3) return false;
    bool has_pos = false, has_neg = false;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        int k = (j + 1) % n;
        int cross = sgn((poly[j] - poly[i]).cross(poly[k] - poly[j]));
        if (cross > 0) has_pos = true;
        if (cross < 0) has_neg = true;
        if (has_pos && has_neg) return false;
    }
    return has_pos ^ has_neg;
}

template<typename T>
Point<ld> polygon_centroid(const Polygon<T>& poly) {
    if (poly.empty()) return Point<ld>(0, 0);
    Point<ld> c(0, 0);
    ld area = 0;
    int n = poly.size();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        ld temp = static_cast<ld>(poly[i].cross(poly[j]));
        area += temp;
        c = c + Point<ld>(poly[i].x + poly[j].x, poly[i].y + poly[j].y) * temp;
    }
    if (sgn(area) == 0) return Point<ld>(poly[0].x, poly[0].y); // 退化情况
    return c / (3.0 * area);
}

template<typename T>
requires std::is_floating_point_v<T>
Polygon<T> polygon_kernel(const Polygon<T>& poly) {
    int n = poly.size();
    if (n < 3) return {};
    std::vector<Line<T>> lines;
    lines.reserve(n);
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        // 边 poly[i] → poly[j] 在 CCW 多边形下,左侧为多边形内部
        lines.push_back(Line<T>(poly[i], poly[j]));
    }
    // 转: HalfPlane·半平面交.hpp::half_plane_intersection
    return half_plane_intersection(lines);
}

template<typename T>
std::pair<Point<T>, Point<T>> farthest_point_pair(const Polygon<T>& convex) {
    int n = convex.size();
    if (n <= 1) return {Point<T>(), Point<T>()};
    if (n == 2) return {convex[0], convex[1]};

    // 找到x最小和最大的点
    int left = 0, right = 0;
    for(int i = 1; i < n; i++) {
        if (convex[i].x < convex[left].x) left = i;
        if (convex[i].x > convex[right].x) right = i;
    }

    T maxDist2 = 0;
    std::pair<Point<T>, Point<T>> res;

    // 旋转卡尺
    for (int i = 0; i < n; i++) {
        int curr = left, next = (left + 1) % n;
        while (sgn((convex[next] - convex[curr]).cross(convex[(right + 1) % n] - convex[right])) > 0) {
            right = (right + 1) % n;
        }

        T d2 = convex[curr].dist2(convex[right]);
        if (d2 > maxDist2) {
            maxDist2 = d2;
            res = {convex[curr], convex[right]};
        }

        left = next;
    }

    return res;
}

}
