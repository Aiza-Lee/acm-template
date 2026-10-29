#pragma once
#include "1-2-PolygonFP·浮点多边形与三角形.hpp"

/*
 * 半平面交与多边形核
 *
 * Overview:
 *     使用排序增量双端队列法（S&I 算法）求多个有向直线左侧半平面的交集，
 *     并提供简单多边形核（能直接看到多边形所有内壁的可见点集区域）的求解。
 *
 * API:
 *     half_plane_inter(lines) — 求有向直线集合的左半平面交（返回逆时针凸多边形顶
 *                                点），空交集或退化为点/线段时返回 {}，O(N log N)。
 *     polygon_kernel(poly)    — 求逆时针简单多边形的核，空核返回 {}，O(N log N)。
 *
 * Notes:
 *     1. 每个半平面定义为有向直线 Line(a, b) 的左侧区域（满足 l.side(p) >= 0）。
 *     2. 极角排序后同向平行直线仅保留最靠左（保留区域最小）的一条；
 *        反向平行半平面矛盾时在队列闭合检查中自动返回空集 {}。
 *     3. 当去重后的方向角出现 >= PI 的空隙（交集无界）时，自动追加 1e6
 *        倍输入坐标量级的包围盒将其裁剪为有限凸多边形；若题目已给定显式边界，
 *        建议直接将 4 条边界加入。
 *     4. polygon_kernel 要求输入多边形顶点按逆时针顺序排列（顺时针输入会导致半平面反
 *        向）。
 *
 * Related:
 *     1-2-PolygonFP·浮点多边形与三角形.hpp::polygon_cut — 单条有向直线 O(N)
 *     切割凸多边形。
 */

namespace Geo2D {

// S&I 双端队列求有向直线集合的左半平面交，返回逆时针凸多边形顶点，O(N log N)
template<typename T = ld>
Polygon<T> half_plane_inter(std::vector<Line<T>> lines) {
    if (lines.size() < 3) return {};

    // 按方向向量极角排序，同向平行直线去重并保留最靠左的一条
    auto sort_dedup = [](std::vector<Line<T>>& in) {
        std::sort(in.begin(), in.end(), [](const Line<T>& x, const Line<T>& y) {
            Point<T> vx = x.v(), vy = y.v();
            if (polar_less(vx, vy)) return true;
            if (polar_less(vy, vx)) return false;
            return x.side(y.a) > 0;
        });
        std::vector<Line<T>> ls;
        ls.reserve(in.size());
        for (const auto& l : in) {
            if (ls.empty() || !same_dir(ls.back(), l)) {
                ls.push_back(l);
            } else if (ls.back().side(l.a) > 0) {
                ls.back() = l;
            }
        }
        return ls;
    };

    auto dir_ang = [](const Line<T>& l) -> T {
        T t = l.v().ang();
        return t < 0 ? t + 2 * PI : t;
    };

    auto ls = sort_dedup(lines);
    if (ls.size() < 3) return {};

    // 无界判定：相邻方向角最大空隙 >= PI 时自动追加大包围盒
    bool unbounded = (dir_ang(ls[0]) + 2 * PI - dir_ang(ls.back()) >= PI - EPS);
    for (size_t i = 1; i < ls.size() && !unbounded; ++i) {
        unbounded = (dir_ang(ls[i]) - dir_ang(ls[i - 1]) >= PI - EPS);
    }
    if (unbounded) {
        T mx = 1;
        for (const auto& l : ls) {
            mx = std::max({mx, std::abs(l.a.x), std::abs(l.a.y),
                               std::abs(l.b.x), std::abs(l.b.y)});
        }
        T B = mx * (T)1e6;
        ls.push_back({{-B, -B}, {B, -B}});
        ls.push_back({{B, -B}, {B, B}});
        ls.push_back({{B, B}, {-B, B}});
        ls.push_back({{-B, B}, {-B, -B}});
        ls = sort_dedup(ls);
        if (ls.size() < 3) return {};
    }

    auto outside = [](const Line<T>& l, Point<T> p) { return l.side(p) < 0; };
    std::deque<Line<T>> q;
    std::deque<Point<T>> p;
    for (const auto& l : ls) {
        while (!p.empty() && outside(l, p.back())) { p.pop_back(); q.pop_back(); }
        while (!p.empty() && outside(l, p.front())) { p.pop_front(); q.pop_front(); }
        if (!q.empty()) {
            if (parallel(q.back(), l)) return {};
            p.push_back(line_inter(q.back(), l));
        }
        q.push_back(l);
    }
    while (!p.empty() && outside(q.front(), p.back())) {
        p.pop_back(); q.pop_back();
    }
    while (!p.empty() && outside(q.back(), p.front())) {
        p.pop_front(); q.pop_front();
    }
    if (q.size() < 3 || parallel(q.back(), q.front())) return {};
    p.push_back(line_inter(q.back(), q.front()));
    return {p.begin(), p.end()};
}

// 求逆时针简单多边形的核，O(N log N)
template<typename T = ld>
Polygon<T> polygon_kernel(const Polygon<T>& poly) {
    int n = poly.size();
    if (n < 3) return {};
    std::vector<Line<T>> lines;
    lines.reserve(n);
    rep(i, 0, n - 1) lines.emplace_back(poly[i], poly[(i + 1) % n]);
    return half_plane_inter(lines);
}

} // namespace Geo2D
